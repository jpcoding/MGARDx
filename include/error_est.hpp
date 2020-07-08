#ifndef _REFACTOR_ERROR_EST_HPP
#define _REFACTOR_ERROR_EST_HPP

#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <bitset>
#include <queue>

namespace REFACTOR{

using namespace std;

#define MAX_2(a, b) (a > b) ? (a) : (b)
// compute maximum value in level
/*
@params data: level data
@params n: number of level data points
*/
template <class T>
T record_level_max_value(const T * data, size_t n){
    T max_val = 0;
    for(int i=0; i<n; i++){
        T val = fabs(data[i]);
        if(val > max_val) max_val = val;
    }
    return max_val;
}

union FloatingInt32{
    float f;
    uint32_t i;
};
union FloatingInt64{
    double f;
    uint64_t i;
};
// compute max_e indicator in level
/*
@params data: level data
@params n: number of level data points
@params num_bitplanes: number of encoded bitplanes (include sign)
@params level_exp: aligned level exponent
*/
template <class T>
vector<double> record_level_max_e(const T * data, size_t n, int num_bitplanes, int level_exp);
template <>
vector<double> record_level_max_e(const float * data, size_t n, int num_bitplanes, int level_exp){
    const int prec = 23;
    int encode_prec = num_bitplanes - 1;
    // if(encode_prec > prec) encode_prec = prec;
    vector<double> max_e = vector<double>(encode_prec + 1, 0);
    FloatingInt32 fi;
    for(int i=0; i<n; i++){
        int data_exp = 0;
        frexp(data[i], &data_exp);
        auto val = data[i];
        fi.f = val;
        int exp_diff = level_exp - data_exp + prec - encode_prec;
        int index = encode_prec;
        if(exp_diff > 0){
            // zeroing out unrecorded bitplanes
            for(int b=0; b<exp_diff; b++){
                fi.i &= ~(1u << b);            
            }
        }
        else{
            // skip padding 0s (no errors)
            index += exp_diff;
            exp_diff = 0;
        }
        for(int b=exp_diff; b<prec; b++){
            // change b-th bit to 0
            fi.i &= ~(1u << b);
            max_e[index] = MAX_2(max_e[index], fabs(data[i] - fi.f));
            index --;
        }
        while(index >= 0){
            max_e[index] = MAX_2(max_e[index], fabs(data[i]));
            index --;
        }
    }
    return max_e;
}
template <>
vector<double> record_level_max_e(const double * data, size_t n, int num_bitplanes, int level_exp){
    cout << "Not implemented yet...\nExit -1.\n";
    exit(-1);
    if(num_bitplanes > 52) num_bitplanes = 52;
    vector<double> mse = vector<double>(num_bitplanes, 0);
    // FloatingInt64 fi;
    return mse;
}
// compute mse indicator in level
/*
@params data: level data
@params n: number of level data points
@params num_bitplanes: number of encoded bitplanes (include sign)
@params level_exp: aligned level exponent
*/
template <class T>
vector<double> record_level_mse(const T * data, size_t n, int num_bitplanes, int level_exp);
template <>
vector<double> record_level_mse(const float * data, size_t n, int num_bitplanes, int level_exp){
    const int prec = 23;
    const int encode_prec = num_bitplanes - 1;
    vector<double> mse = vector<double>(num_bitplanes, 0);
    FloatingInt32 fi;
    for(int i=0; i<n; i++){
        int data_exp = 0;
        frexp(data[i], &data_exp);
        auto val = data[i];
        fi.f = val;
        int exp_diff = level_exp - data_exp + prec - encode_prec;
        int index = encode_prec;
        if(exp_diff > 0){
            // zeroing out unrecorded bitplanes
            for(int b=0; b<exp_diff; b++){
                fi.i &= ~(1u << b);            
            }
        }
        else{
            // skip padding 0s (no errors)
            index += exp_diff;
            exp_diff = 0;
        }
        for(int b=exp_diff; b<prec; b++){
            // change b-th bit to 0
            fi.i &= ~(1u << b);
            mse[index] += (data[i] - fi.f)*(data[i] - fi.f);
            index --;
        }
        while(index >= 0){
            mse[index] += data[i] * data[i];
            index --;
        }
    }
    return mse;
}
template <>
vector<double> record_level_mse(const double * data, size_t n, int num_bitplanes, int level_exp){
    cout << "Not implemented yet...\nExit -1.\n";
    exit(-1);
    if(num_bitplanes > 52) num_bitplanes = 52;
    vector<double> mse = vector<double>(num_bitplanes, 0);
    // FloatingInt64 fi;
    return mse;
}

struct Efficiency{
    double efficiency;
    int level;
    Efficiency(double e, int l) : efficiency(e), level(l) {}
};
struct CompareEfficiency { 
    bool operator()(Efficiency const& e1, Efficiency const& e2) { 
        return e1.efficiency < e2.efficiency; 
    } 
}; 
// reorganize refactored data with respect to error estimators
// use a greedy algorithm to pick up the most efficient bitplane
// return reorganized data
/*
@params level_components: level bitplanes
@params level_sizes: level encoded sizes
@params level_errors: level error estimators
@params reordered_size: size of aggregated data
*/
unsigned char * refactored_data_reorganization(const vector<vector<unsigned char*>>& level_components, const vector<vector<size_t>>& level_sizes, const vector<vector<double>>& level_errors, size_t& reordered_size){
    const int num_levels = level_components.size();
    size_t total_size = 0;
    // init error_gain: reduced error by including current bitplane
    // init sizes: the corresponding sizes of bitplanes
    // NOTE: the sizes of the two vector is 1 less than component sizes
    //         because sign is grouped with the first bitplane
    vector<vector<double>> error_gain;
    vector<vector<double>> sizes;
    for(int i=0; i<num_levels; i++){
        error_gain.push_back(vector<double>(level_components[i].size() - 1));
        sizes.push_back(vector<double>(level_components[i].size() - 1));
    }
    // compute erorr gain and sizes for bitplanes
    for(int i=0; i<num_levels; i++){
        for(int j=0; j<error_gain[i].size(); j++){
            error_gain[i][j] = level_errors[i][j] - level_errors[i][j+1];
            sizes[i][j] = (j == 0) ? (level_sizes[i][0] + level_sizes[i][1]) : level_sizes[i][j+1];
        }
    }
    for(int i=0; i<num_levels; i++){
        for(int j=0; j<error_gain[i].size(); j++){
            cout << error_gain[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    for(int i=0; i<num_levels; i++){
        for(int j=0; j<error_gain[i].size(); j++){
            cout << sizes[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    vector<size_t> index(num_levels, 0);
    // metric for greedy algorithm: how much error gain per byte
    priority_queue<Efficiency, vector<Efficiency>, CompareEfficiency> efficiency_heap;
    for(int i=0; i<num_levels; i++){
        efficiency_heap.push(Efficiency(error_gain[i][0] / sizes[i][0], i));
    }
    while(!efficiency_heap.empty()){
        auto eff = efficiency_heap.top();
        efficiency_heap.pop();
        auto level = eff.level;
        cout << "Encode level " << level << " component " << index[level] << ", efficiency = " << eff.efficiency << endl;
        index[level] ++;
        if(index[level] != level_components[level].size() - 1){
            efficiency_heap.push(Efficiency(error_gain[level][index[level]] / sizes[level][index[level]], level));
        }
    }
    exit(0);
    return NULL;
}

}
#endif