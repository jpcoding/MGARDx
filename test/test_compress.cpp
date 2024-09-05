#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include "utils.hpp"
#include "decompose.hpp"
#include "recompose.hpp"

using namespace std;

template <class T>
unsigned char * test_compress(vector<T>& data, const vector<size_t>& dims, int target_level, double eb, size_t& compressed_size, bool use_sz){
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    MGARD::Decomposer<T> decomposer(use_sz);
    auto compressed_data = decomposer.compress(data.data(), dims, target_level, eb, compressed_size);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Compression time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    return compressed_data;
}

template <class T>
T * test_decompress(const unsigned char * compressed, size_t compressed_size, const vector<size_t>& dims){
    struct timespec start, end;
    int err = 0;
    err = clock_gettime(CLOCK_REALTIME, &start);
    MGARD::Recomposer<T> recomposer;
    auto data_dec = recomposer.decompress(compressed, compressed_size, dims);
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Decompression time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    return data_dec;
}

template <class T>
void test(string ori_filename, string compressed_filename, const vector<size_t>& dims){
    size_t num_elements = 0;
    size_t compressed_size = 0;
    auto compressed = MGARD::readfile<unsigned char>(compressed_filename.c_str(), compressed_size);
    auto data_ori = MGARD::readfile<T>(ori_filename.c_str(), num_elements);
    auto data_dec = test_decompress<T>(compressed.data(), compressed_size, dims);
    MGARD::print_statistics(data_ori.data(), data_dec, num_elements, compressed_size);
    MGARD::writefile((compressed_filename + ".out").c_str(), data_dec, num_elements);
    free(data_dec);    
}



int main(int argc, char ** argv){
    
    string filename = string(argv[1]);
    int type = atoi(argv[2]); // 0 for float, 1 for double

    double rel_eb = atof(argv[3]);
   
    double tolerance = atof(argv[3]);


    const int target_level = atoi(argv[4]);
    bool use_sz = atoi(argv[5]);    // whether to use sz
    const int num_dims = atoi(argv[6]);
    vector<size_t> dims(num_dims);
    for(int i=0; i<dims.size(); i++){
       dims[i] = atoi(argv[7 + i]);
       cout << dims[i] << " ";
    }
    cout << endl;
    double eb = tolerance;
    cerr << "Data = " << filename << endl;
    // cerr << "Required eb = " << tolerance << endl;
    std::cout << "rel_eb = " << rel_eb << std::endl; 
    size_t num_elements = 0;
    size_t compressed_size = 0;
    unsigned char * compressed = NULL;
    int data_size = 0;
    switch(type){
        case 0:
            {
                auto data = MGARD::readfile<float>(filename.c_str(), num_elements);
                
                eb = rel_eb * (*std::max_element(data.begin(), data.end()) - *std::min_element(data.begin(), data.end()));
                std::cout<<"abs eb = "<<eb<<std::endl;
                compressed = test_compress(data, dims, target_level, eb, compressed_size, use_sz);
                data_size = sizeof(float);  
                break;
            }
        case 1:
            {
                auto data = MGARD::readfile<double>(filename.c_str(), num_elements);
                eb = rel_eb * (*std::max_element(data.begin(), data.end()) - *std::min_element(data.begin(), data.end()));
                std::cout<<"abs eb = "<<eb<<std::endl;
                compressed = test_compress(data, dims, target_level, eb, compressed_size, use_sz);
                data_size = sizeof(double); 
                break;
            }
        default:
            cerr << "Only 0 (float) and 1 (double) are implemented in this test\n";
            exit(0);
    }
    std::cout<<"Compressed size = "<<compressed_size<<std::endl;
    std::cout<<"Compression ratio = "<<(double)num_elements*data_size/compressed_size<<std::endl;

    MGARD::writefile((filename + ".mgardpred").c_str(), compressed, compressed_size);
    free(compressed);


    switch(type){
    case 0:
        {
            test<float>(filename, (filename + ".mgardpred"), dims);
            break;
        }
    case 1:
        {
            test<double>(filename, (filename + ".mgardpred"), dims);
            break;
        }
    default:
            cerr << "Only 0 (float) and 1 (double) are implemented in this test\n";
            exit(0);
    }








    return 0;
}