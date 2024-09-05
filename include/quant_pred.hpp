#ifndef _MGARD_QUANT_PRED
#define _MGARD_QUANT_PRED

#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/quantizer/IntegerQuantizer.hpp"
#include "utils.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <iostream>


namespace SZ = SZ3;

namespace MGARD {

template <class T> class QuantPred {
public:
  QuantPred(int* quant_inds, int N, const size_t* dims, int start_level, int unpred, int quant_radius, int decompse = 1) {
    N = N; 
    num_elements = 1;
    for(int i=0; i<N; i++){
      num_elements *= dims[i];
      global_dims.push_back(dims[i]);
      // std::cout << "global_dims[" << i << "] = " << global_dims[i] << std::endl;  
    }
    if(decompse == 1)
    {
      aux_quant_inds.resize(num_elements);
      std::copy(quant_inds, quant_inds + num_elements, aux_quant_inds.begin());
    }
    orig_quant_inds = quant_inds;
    unpred = unpred; 
    quant_radius = quant_radius; 
    init();
  }

  ~QuantPred() {} 

  inline int  get_quant_compensate (size_t idx, size_t stride1, size_t stride2)
  {
    //   if(orig_quant_inds[idx] == unpred){
    //   return 0;
    //  }
     int A = aux_quant_inds[idx - stride1 - stride2];
     int B = aux_quant_inds[idx - stride1];
     int C = aux_quant_inds[idx - stride2];
     if(A == unpred || B == unpred || C == unpred){
      return 0;
     }
     if( B > quant_radius && C > quant_radius){
       return B+C-A-quant_radius;
     }
     else if (B< quant_radius && C < quant_radius)
     {
      return B+C-A-quant_radius; 
     }
     else{
      return 0;
    }
  }

  inline int  get_quant_compensate_recover (size_t idx, size_t stride1, size_t stride2)
  {
    //  if(orig_quant_inds[idx] == unpred){
    //   return 0;
    //  }
     int A = orig_quant_inds[idx - stride1 - stride2];
     int B = orig_quant_inds[idx - stride1];
     int C = orig_quant_inds[idx - stride2];
     if(A == unpred || B == unpred || C == unpred){
      return 0;
     }
     if( B > quant_radius && C > quant_radius){
       return B+C-A-quant_radius;
     }
     else if (B< quant_radius && C < quant_radius)
     {
      return B+C-A-quant_radius; 
     }
     else{
      return 0;
    }
  }



  void quant_pred_level_3D(int target_level)
  {
    // if(target_level > start_level){
    //   return;
    // }
    int level_stride = 1<<(target_level -1);
    // std::cout << "level_stride = " << level_stride << std::endl;
    int dim_starts[3] = {0, 0, 0};
    size_t dim0_stride = global_dims[1] * global_dims[2];
    size_t dim1_stride = global_dims[2];
    size_t dim2_stride = 1;

    // std::cout << "dim0_stride = " << dim0_stride << std::endl;
    // std::cout << "dim1_stride = " << dim1_stride << std::endl;
    // std::cout << "dim2_stride = " << dim2_stride << std::endl;

    int quant_compensate = 0 ;

    size_t cur_pos = 0;
    std::vector<int> push_data; 

    auto timer = MGARD::Timer();
    timer.start();

    // 1. first reorganize the loop order
    for (int i = 2*level_stride; i < global_dims[0]-level_stride; i+=level_stride*2){
      for(int j = 2*level_stride; j < global_dims[1]-level_stride; j+=level_stride*2 ){
        for(int k = level_stride; k < global_dims[2]; k+=2*level_stride){
          cur_pos = i * dim0_stride + j * dim1_stride + k * dim2_stride; 
          quant_compensate = get_quant_compensate(cur_pos, dim0_stride*2*level_stride, dim1_stride*2*level_stride);
          orig_quant_inds[cur_pos] = orig_quant_inds[cur_pos] -  quant_compensate;
        }
      }
    }
    // timer.stop("first dim");

    // 2. second dim 
    for (int i = level_stride*2; i < global_dims[0]-level_stride; i+=2*level_stride){
      for (int j = level_stride; j < global_dims[1]; j+=2*level_stride){
        for(int k = level_stride; k < global_dims[2]-level_stride; k+=level_stride){
          cur_pos = i * dim0_stride + j * dim1_stride + k * dim2_stride; 
          quant_compensate = get_quant_compensate(cur_pos, dim0_stride*2*level_stride, dim2_stride*level_stride);
          orig_quant_inds[cur_pos] = orig_quant_inds[cur_pos] -  quant_compensate;
        }
      }

    }

    // 3. third dim
    for(int i = level_stride; i < global_dims[0]; i+=2*level_stride){
      for(int j = level_stride; j < global_dims[1]-level_stride; j+=level_stride ){
        for(int k = level_stride; k < global_dims[2]-level_stride; k+=level_stride){
          cur_pos = i * dim0_stride + j * dim1_stride + k * dim2_stride; 
          quant_compensate = get_quant_compensate(cur_pos, dim1_stride*level_stride, dim2_stride*level_stride);
          orig_quant_inds[cur_pos] = orig_quant_inds[cur_pos] -  quant_compensate;
        }
      }
    }


    // writefile("push_data.dat", push_data.data(), push_data.size());
  }

  void quant_pred_level_3D_recover(int target_level)
  {
    // if(target_level > start_level){
    //   return;
    // }
    int level_stride = 1<<(target_level -1);
    // std::cout << "[recover]level_stride = " << level_stride << std::endl;
    int dim_starts[3] = {0, 0, 0};
    size_t dim0_stride = global_dims[1] * global_dims[2];
    size_t dim1_stride = global_dims[2];
    size_t dim2_stride = 1;

    // std::cout << "[recover]dim0_stride = " << dim0_stride << std::endl;
    // std::cout << "[recover]dim1_stride = " << dim1_stride << std::endl;
    // std::cout << "[recover]dim2_stride = " << dim2_stride << std::endl;

    int quant_compensate = 0 ;

    size_t cur_pos = 0;
    
    auto timer = MGARD::Timer();
    timer.start();

    // 1. first reorganize the loop order
    for (int i = 2*level_stride; i < global_dims[0]-level_stride; i+=level_stride*2){
      for(int j = 2*level_stride; j < global_dims[1]-level_stride; j+=level_stride*2 ){
        for(int k = level_stride; k < global_dims[2]; k+=2*level_stride){
          cur_pos = i * dim0_stride + j * dim1_stride + k * dim2_stride; 
          quant_compensate = get_quant_compensate_recover(cur_pos, dim0_stride*2*level_stride, dim1_stride*2*level_stride);
          orig_quant_inds[cur_pos] = orig_quant_inds[cur_pos] +  quant_compensate;
        }
      }
    }
    // timer.stop("first dim");

    // 2. second dim 
    for (int i = level_stride*2; i < global_dims[0]-level_stride; i+=2*level_stride){
      for (int j = level_stride; j < global_dims[1]; j+=2*level_stride){
        for(int k = level_stride; k < global_dims[2]-level_stride; k+=level_stride){
          cur_pos = i * dim0_stride + j * dim1_stride + k * dim2_stride; 
          quant_compensate = get_quant_compensate_recover(cur_pos, dim0_stride*2*level_stride, dim2_stride*level_stride);
          orig_quant_inds[cur_pos] = orig_quant_inds[cur_pos] +  quant_compensate;
        }
      }

    }

    // 3. third dim
    for(int i = level_stride; i < global_dims[0]; i+=2*level_stride){
      for(int j = level_stride; j < global_dims[1]-level_stride; j+=level_stride ){
        for(int k = level_stride; k < global_dims[2]-level_stride; k+=level_stride){
          cur_pos = i * dim0_stride + j * dim1_stride + k * dim2_stride; 
          quant_compensate = get_quant_compensate_recover(cur_pos, dim1_stride*level_stride, dim2_stride*level_stride);
          orig_quant_inds[cur_pos] = orig_quant_inds[cur_pos] +  quant_compensate;
        }
      }
    }

  }
 
private:



  void init() 
  {

  }

  

  size_t num_elements; 
  int N; 
  std::vector<size_t>  global_dims; 
  int start_level = -1;
  std::vector<int> aux_quant_inds; 
  int* orig_quant_inds; 
  int unpred = 0;
  int quant_radius = 2<<15;
};

} // namespace MGARD
#endif
