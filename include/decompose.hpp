#ifndef _MGARD_DECOMPOSE_HPP
#define _MGARD_DECOMPOSE_HPP

#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include "utils.hpp"

namespace MGARD{

using namespace std;

template <class T>
class Decomposer{
public:
	Decomposer(){};
	~Decomposer(){
		if(data_buffer) free(data_buffer);
		if(correction_buffer) free(correction_buffer);	
		if(load_v_buffer) free(load_v_buffer);
	};
	void decompose(T * data_, const vector<size_t>& dims, size_t target_level){
		data = data_;
		size_t num_elements = 1;
		for(const auto& d:dims){
			num_elements *= d;
		}
		data_buffer_size = num_elements * sizeof(T);
		init(dims);
		if(dims.size() == 1){
			size_t h = 1;
			size_t n = dims[0];
			for(int i=0; i<target_level; i++){
				decompose_level_1D(data, n, h);
				n = (n >> 1) + 1;
				h <<= 1;
			}
		}
		else if(dims.size() == 2){
			size_t h = 1;
			size_t n1 = dims[0];
			size_t n2 = dims[1];
			for(int i=0; i<target_level; i++){
				decompose_level_2D(data, n1, n2, (T)h, dims[1]);
				n1 = (n1 >> 1) + 1;
				n2 = (n2 >> 1) + 1;
				h <<= 1;
			}
		}
	}

private:
	unsigned int default_batch_size = 32;
	size_t data_buffer_size = 0;
	T * data = NULL;			// pointer to the original data
	T * data_buffer = NULL;		// buffer for reordered data
	T * load_v_buffer = NULL;
	T * correction_buffer = NULL;

	void init(const vector<size_t>& dims){
		size_t buffer_size = default_batch_size * (*max_element(dims.begin(), dims.end())) * sizeof(T);
		cerr << "buffer_size = " << buffer_size << endl;
		if(data_buffer) free(data_buffer);
		if(correction_buffer) free(correction_buffer);
		if(load_v_buffer) free(load_v_buffer);
		data_buffer = (T *) malloc(data_buffer_size);
		correction_buffer = (T *) malloc(buffer_size);
		load_v_buffer = (T *)malloc(buffer_size);
	}
	// reorder the data to put all the coefficient to the back
	void data_reorder_1D(const T * data_pos, size_t n_nodal, size_t n_coeff, T * nodal_buffer, T * coeff_buffer){
		T * nodal_pos = nodal_buffer;
		T * coeff_pos = coeff_buffer;
		T const * cur_data_pos = data_pos;
		for(int i=0; i<n_coeff; i++){
			*(nodal_pos++) = *(cur_data_pos++);
			*(coeff_pos++) = *(cur_data_pos++);
		}
		*(nodal_pos++) = *(cur_data_pos++);
		if(n_nodal == n_coeff + 2){
			// if even, add a nodal value such that the interpolant
			// of the last two nodal values equal to the last coefficient
			*nodal_pos = 2*cur_data_pos[0] - nodal_pos[-1];
		}
	}
	// reorder the data to put all the coefficient to the back
	/*
		oxoxo		oooxx		oooxx
		xxxxx	(1)	xxxxx	(2)	oooxx
		oxoxo	=>	oooxx	=>	oooxx
		xxxxx		xxxxx		xxxxx
		oxoxo		oooxx		xxxxx
	*/
	void data_reorder_2D(T * data_pos, size_t n1_nodal, size_t n1_coeff, size_t n2_nodal, size_t n2_coeff, size_t stride){
		size_t n1 = n1_nodal + n1_coeff;
		size_t n2 = n2_nodal + n2_coeff;
		T * cur_data_pos = data_pos;
		T * nodal_pos = data_buffer;
		T * coeff_pos = data_buffer + n2_nodal;
		// do reorder (1)
		for(int i=0; i<n1; i++){
			data_reorder_1D(cur_data_pos, n2_nodal, n2_coeff, nodal_pos, coeff_pos);
			memcpy(cur_data_pos, data_buffer, n2 * sizeof(T));
			cur_data_pos += stride;
		}
		if(!(n1 & 1)){
			// n1 is even, change the last coeff row into nodal row
			cur_data_pos -= stride;
			for(int j=0; j<n2; j++){
				cur_data_pos[j] = 2 * cur_data_pos[j] - cur_data_pos[-stride + j];
			}
		}
		// do reorder (2)
		// TODO: change to online processing for memory saving
		switch_rows_2D_by_buffer(data_pos, data_buffer, n1_nodal + n1_coeff, n2_nodal + n2_coeff, stride);
	}
	// compute the difference between original value 
	// and interpolant (I - PI_l)Q_l
	// overwrite the data in N_l \ N_(l-1) in place
	void compute_interpolant_difference_1D(size_t n_coeff, const T * nodal_buffer, T * coeff_buffer){
		for(int i=0; i<n_coeff; i++){
			coeff_buffer[i] -= (nodal_buffer[i] + nodal_buffer[i+1]) / 2; 
		}
	}
	void add_correction(size_t n_nodal, T * nodal_buffer){
		for(int i=0; i<n_nodal; i++){
			nodal_buffer[i] += correction_buffer[i];
		}
	}
	// decompose a level with n element and the given stride
	// to a level with n/2 element
	void decompose_level_1D(T * data_pos, size_t n, T h){
		size_t n_nodal = (n >> 1) + 1;
		size_t n_coeff = n - n_nodal;
		T * nodal_buffer = data_buffer;
		T * coeff_buffer = data_buffer + n_nodal;
		data_reorder_1D(data_pos, n_nodal, n_coeff, nodal_buffer, coeff_buffer);
		compute_interpolant_difference_1D(n_coeff, nodal_buffer, coeff_buffer);
		compute_load_vector(load_v_buffer, n_nodal, n_coeff, h, coeff_buffer);
		compute_correction(correction_buffer, n_nodal, h, load_v_buffer);
		add_correction(n_nodal, nodal_buffer);
		memcpy(data_pos, data_buffer, n*sizeof(T));
	}
	// compute and add the corrections
	void compute_and_add_correction_2D(T * data_pos, size_t n1_nodal, size_t n1_coeff, size_t n2_nodal, size_t n2_coeff, T h, size_t stride){
		size_t n1 = n1_nodal + n1_coeff;
		size_t n2 = n2_nodal + n2_coeff;
		// compute horizontal correction
		T * nodal_pos = data_pos;
		const T * coeff_pos = data_pos + n2_nodal;
		// store horizontal corrections in the data_buffer
		T * correction_pos = data_buffer;
		for(int i=0; i<n1; i++){
			compute_load_vector(load_v_buffer, n2_nodal, n2_coeff, h, coeff_pos);
			compute_correction(correction_pos, n2_nodal, h, load_v_buffer);
			// print(correction_pos, 1, n2_nodal, "horizontal_correction");
			// add_correction(n2_nodal, nodal_pos);
			nodal_pos += stride, coeff_pos += stride;
			correction_pos += n2_nodal;
		}
		// compute vertical correction
		compute_and_apply_correction_2D_vertical(data_pos, n1, n2, h, stride, data_buffer, load_v_buffer, correction_buffer, default_batch_size, true);
	}
	// compute the difference between original value 
	// and interpolant (I - PI_l)Q_l for the coefficient rows in 2D
	// overwrite the data in N_l \ N_(l-1) in place
	// Note: interpolant difference in the nodal rows have already been computed
	void compute_interpolant_difference_2D_vertical(T * data_pos, size_t n1, size_t n2, size_t stride){
		size_t n1_nodal = (n1 >> 1) + 1;
		size_t n1_coeff = n1 - n1_nodal;
		size_t n2_nodal = (n2 >> 1) + 1;
		size_t n2_coeff = n2 - n2_nodal;
		bool even_n2 = !(n2 & 1);
		T * n1_nodal_data = data_pos;
		T * n1_coeff_data = data_pos + n1_nodal * stride;
		for(int i=0; i<n1_coeff; i++){
            const T * nodal_pos = n1_nodal_data + i * stride;
            T * coeff_pos = n1_coeff_data + i * stride;
            // TODO: optimize average computation
            T * nodal_coeff_pos = coeff_pos;	// coeffcients in nodal rows
            T * coeff_coeff_pos = coeff_pos + n2_nodal;	// coefficients in coeffcients rows
            for(int j=0; j<n2_coeff; j++){
                // coefficients in nodal columns
                *(nodal_coeff_pos++) -= (nodal_pos[j] + nodal_pos[stride + j]) / 2;
                // coefficients in centers
                *(coeff_coeff_pos++) -= (nodal_pos[j] + nodal_pos[j + 1] + nodal_pos[stride + j] + nodal_pos[stride + j + 1]) / 4;
            }
            // compute the last (or second last if n2 is even) nodal column
            *(nodal_coeff_pos ++) -= (nodal_pos[n2_coeff] + nodal_pos[stride + n2_coeff]) / 2;
            if(even_n2){
                // compute the last nodal column
                *(nodal_coeff_pos ++) -= (nodal_pos[n2_coeff + 1] + nodal_pos[stride + n2_coeff + 1]) / 2;
            }
		}
	}
	void compute_interpolant_difference_2D(T * data_pos, size_t n1_nodal, size_t n1_coeff, size_t n2_nodal, size_t n2_coeff, size_t stride){
		size_t n1 = n1_nodal + n1_coeff;
		size_t n2 = n2_nodal + n2_coeff;
		// compute horizontal difference
		const T * nodal_pos = data_pos;
		T * coeff_pos = data_pos + n2_nodal;
		for(int i=0; i<n1_nodal; i++){
			compute_interpolant_difference_1D(n2_coeff, nodal_pos, coeff_pos);
			nodal_pos += stride, coeff_pos += stride;
		}
		// compute vertical difference
		compute_interpolant_difference_2D_vertical(data_pos, n1, n2, stride);
	}	
	// decompose n1 x n2 data into coarse level (n1/2 x n2/2)
	void decompose_level_2D(T * data_pos, size_t n1, size_t n2, T h, size_t stride){
		cerr << "decompose, h = " << h << endl; 
		size_t n1_nodal = (n1 >> 1) + 1;
		size_t n1_coeff = n1 - n1_nodal;
		size_t n2_nodal = (n2 >> 1) + 1;
		size_t n2_coeff = n2 - n2_nodal;
		// if(h > 1) print(data, n1, n2, "data_entry");
		data_reorder_2D(data_pos, n1_nodal, n1_coeff, n2_nodal, n2_coeff, stride);
		compute_interpolant_difference_2D(data_pos, n1_nodal, n1_coeff, n2_nodal, n2_coeff, stride);
		// if(h > 1) print(data, n1, n2, "data_before_correction");
		compute_and_add_correction_2D(data_pos, n1_nodal, n1_coeff, n2_nodal, n2_coeff, h, stride);
		// if(h > 1) print(data, n1, n2, "data_after_correction");
	}
};


}

#endif