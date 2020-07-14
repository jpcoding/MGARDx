#ifndef _REFACTOR_DATA_ENC_HPP
#define _REFACTOR_DATA_ENC_HPP

#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <bitset>
#include <iomanip>
#include "data_enc_opt.hpp"

namespace REFACTOR{

using namespace std;

#define ENCODING_DEFAULT 0
#define ENCODING_RLE 1
#define ENCODING_HYBRID 2

class EncoderInterface{
public:
    virtual void encode(bool) = 0;
    virtual void flush() = 0;
    virtual size_t size() = 0;
    virtual unsigned char * save() = 0;
};

class DecoderInterface{
public:
    virtual ~DecoderInterface() = default;
    virtual bool decode() = 0;
    virtual void load(const unsigned char *) = 0;
};
// A class to read data bit by bit
// Currently does not consider the case when decoding size is larger than capacity
class BitDecoder : public DecoderInterface{
public:
    BitDecoder() = default;
    ~BitDecoder() = default;
    bool decode(){
        if (!(buffer >> 1u)) buffer = 0x100u + *current++;
        bool bit = buffer & 1u;
        buffer >>= 1u;
        return bit;
    }
    size_t size(){ return (buffer == 1u) ? current - start : current - start + 1; }
    void load(const unsigned char * encoded_data){
        start = encoded_data;
        current = start;
        buffer = 1u;
    }
private:
    unsigned char const * start = NULL;
    unsigned char const * current = NULL;
    unsigned int buffer = 1u;
};
/*******************************************/

// Runlength encoder
class RunlengthEncoder : public EncoderInterface{
public:
    RunlengthEncoder(){}
    void encode(bool bit){
        if(lastbit == bit){
            count ++;
            if(count == 256){
                length.push_back(count);
                count = 0;
                lastbit = !bit;
            }
        }
        else {
            length.push_back(count);
            count = 1;
            lastbit = bit;
        }
    }
    void flush(){
        if(count != 0){
            length.push_back(count);
            count = 0;
            lastbit = false;
        }
    }
    size_t size(){
        return encoded_size;
    }
    unsigned char * save(){
        // Huffman
        unsigned char * encoded = (unsigned char *) malloc(length.size() * sizeof(int));
        auto encoded_pos = encoded;
        *reinterpret_cast<size_t*>(encoded_pos) = length.size();
        encoded_pos += sizeof(size_t);
        auto encoder = SZ::HuffmanEncoder<int>();
        encoder.preprocess_encode(length, 2*256);
        encoder.save(encoded_pos);
        encoder.encode(length, encoded_pos);
        encoder.postprocess_encode();
        encoded_size = encoded_pos - encoded;
        // return encoded;
        unsigned char * lossless_encoded = lossless_encode(encoded, encoded_size);
        free(encoded);
        return lossless_encoded;
    }
private:
    unsigned char * lossless_encode(unsigned char * encoded, size_t& encoded_size){
        unsigned char * lossless_compressed = NULL;
        size_t lossless_length = MGARD::sz_lossless_compress(ZSTD_COMPRESSOR, 3, encoded, encoded_size, &lossless_compressed);
        encoded_size = lossless_length + sizeof(size_t);
        unsigned char * lossless_compressed_with_size = (unsigned char *) malloc(encoded_size);
        *reinterpret_cast<size_t*>(lossless_compressed_with_size) = lossless_length;
        memcpy(lossless_compressed_with_size + sizeof(size_t), lossless_compressed, lossless_length);
        free(lossless_compressed);
        return lossless_compressed_with_size;
    }
    vector<int> length;
    bool lastbit = false;
    int count = 0;
    int index = 0;
    size_t encoded_size = 0;
};

class RunlengthDecoder : public DecoderInterface{
public:
    RunlengthDecoder() = default;
    ~RunlengthDecoder() = default;
    bool decode(){
        if(count){
            count --;
            return lastbit;
        }
        else{
            count = length[index ++];
            lastbit = !lastbit;
            return decode();
        }
    }
    // void load(const unsigned char * encoded){
    void load(const unsigned char * lossless_encoded){
        auto encoded = lossless_decode(lossless_encoded);
        const unsigned char * encoded_pos = encoded;
        size_t n = *reinterpret_cast<const size_t*>(encoded_pos);
        encoded_pos += sizeof(size_t);
        auto encoder = SZ::HuffmanEncoder<int>();
        size_t remaining_length = INT_MAX;
        encoder.load(encoded_pos, remaining_length);
        length = encoder.decode(encoded_pos, n);
        encoder.postprocess_decode();
        free(encoded);
        count = 0;
        index = 0;
        // toggle lastbit to true such that the first decode would be false since count=0
        lastbit = true;
    }
private:
    unsigned char * lossless_decode(const unsigned char * lossless_encoded){
        unsigned char * encoded = NULL;
        size_t lossless_length = *reinterpret_cast<const size_t*>(lossless_encoded);
        size_t compressed_length = MGARD::sz_lossless_decompress(ZSTD_COMPRESSOR, lossless_encoded + sizeof(size_t), lossless_length, &encoded);
        return encoded;
    }
    vector<int> length;
    bool lastbit = false;
    int count = 0;
    int index = 0;
};

// encode the intra level components progressively
/*
@params data: coefficient data
@params n: number of coefficients in current level
@params level_exp: exponent of max level element
@params num_level_component: number of encoded bitplanes
@params encoded_sizes: size of encoded data
*/
template <class T>
vector<unsigned char*> progressive_encoding(T const * data, size_t n, int level_exp, int num_level_component, vector<size_t>& encoded_sizes){
    vector<unsigned char*> intra_level_components;
    size_t level_component_size = (n * sizeof(T) - 1) / num_level_component + 1 + 8;
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    cout << "level_component_size = " << level_component_size << endl;
    vector<unsigned char *> byte_encoders;
    for(int i=0; i<num_level_component; i++){
        unsigned char * buffer = (unsigned char *) malloc(level_component_size);
        intra_level_components.push_back(buffer);
        byte_encoders.push_back(buffer);
    }    
    size_t buffer_index = byte_wise_direct_encoding_unrolled(data, n, level_exp, num_level_component, byte_encoders);
    for(int k=0; k<num_level_component; k++){
        encoded_sizes.push_back(buffer_index);
    }
    return intra_level_components;
}


template <class T>
T * progressive_decoding(const vector<const unsigned char*>& level_components, size_t n, int level_exp, int num_level_component){
    T * level_data = (T *) malloc(n * sizeof(T));
    size_t level_component_size = (n * sizeof(T) - 1) / num_level_component + 1 + 8;
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    cout << "level_component_size = " << level_component_size << endl;
    vector<BitDecoder> decoders;
    // vector<bitstream*> decoders;
    for(int i=0; i<num_level_component; i++){
        decoders.push_back(BitDecoder());
        decoders[i].load(level_components[i]);
        // decoders.push_back(stream_open(level_components[i], level_component_size));
    }
    T * data_pos = level_data;
    for(int i=0; i<n; i++){
        // decode each bit of the data for each level component
        bool sign = decoders[0].decode();
        // bool sign = stream_read_bit(decoders[0]);
        unsigned int fp = 0;
        for(int j=1; j<num_level_component; j++){
            unsigned int current_bit = decoders[j].decode();
            // unsigned int current_bit = stream_read_bit(decoders[j]);
            fp = (fp << 1) + current_bit;
        }
        long int fix_point = fp;
        if(sign) fix_point = -fix_point;
        *data_pos = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        data_pos ++;
    }
    // for(int i=0; i<num_level_component; i++){
    //     stream_close(decoders[i]);
    // }
    return level_data;
}

// encode the intra level components progressively, with runlength encoding on each bit-plane
/*
@params data: coefficient data
@params n: number of coefficients in current level
@params level_exp: exponent of max level element
@params num_level_component: number of encoded bitplanes
@params encoded_sizes: size of encoded data
*/
template <class T>
vector<unsigned char*> progressive_encoding_with_rle_compression(T const * data, size_t n, int level_exp, int num_level_component, vector<size_t>& encoded_sizes){
    vector<unsigned char*> intra_level_components;
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    struct timespec start, end;
    int err = clock_gettime(CLOCK_REALTIME, &start);        
    vector<RunlengthEncoder> encoders;
    for(int i=0; i<num_level_component; i++){
        encoders.push_back(RunlengthEncoder());
    }
    for(int i=0; i<n; i++){
        T cur_data = ldexp(data[i], num_level_component - 1 - level_exp);
        long int fix_point = (long int) cur_data;
        // encode each bit of the data for each level component
        bool sign = data[i] < 0;
        unsigned int fp = sign ? -fix_point : +fix_point;
        encoders[0].encode(sign);
        for(int j=num_level_component - 1; j>0; j--){
            encoders[j].encode(fp & 1);
            fp >>= 1;
        }
    }
    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "RLE encoding time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
    size_t count = 0;
    for(int i=0; i<num_level_component; i++){
        encoders[i].flush();
        err = clock_gettime(CLOCK_REALTIME, &start);        
        intra_level_components.push_back(encoders[i].save());
        err = clock_gettime(CLOCK_REALTIME, &end);
        encoded_sizes.push_back(encoders[i].size());
        count += encoders[i].size();
        // cout << "The " << i << "-th bitplanes encoding time = " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << " s, ratio = " << (n / 8) * 1.0 /encoders[i].size()<< " , progressive ratio = " << ((i+1) * (n / 8)) * 1.0 / count << endl;
    }
    return intra_level_components;
}

template <class T>
T * progressive_decoding_with_rle_compression(const vector<const unsigned char*>& level_components, size_t n, int level_exp, int num_level_component){
    T * level_data = (T *) malloc(n * sizeof(T));
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    vector<RunlengthDecoder> decoders;
    for(int i=0; i<num_level_component; i++){
        decoders.push_back(RunlengthDecoder());
        decoders[i].load(level_components[i]);
    }
    T * data_pos = level_data;
    for(int i=0; i<n; i++){
        // decode each bit of the data for each level component
        bool sign = decoders[0].decode();
        unsigned int fp = 0;
        for(int j=1; j<num_level_component; j++){
            unsigned int current_bit = decoders[j].decode();
            fp = (fp << 1) + current_bit;
        }
        long int fix_point = fp;
        if(sign) fix_point = -fix_point;
        *data_pos = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        data_pos ++;
    }
    return level_data;
}

// encode the intra level components progressively using alternatives of direct encoding and rle encoding
/*
@params data: coefficient data
@params n: number of coefficients in current level
@params level_exp: exponent of max level element
@params num_level_component: number of encoded bitplanes
@params encoded_sizes: size of encoded data
@params bitplane_indicator: indicator for encoder selection
*/
template <class T>
vector<unsigned char*> progressive_hybrid_encoding(T const * data, size_t n, int level_exp, int num_level_component, vector<size_t>& encoded_sizes, vector<unsigned char>& bitplane_indicator){
    vector<unsigned char*> intra_level_components;
    size_t level_component_size = (n * sizeof(T) - 1) / num_level_component + 1 + 8;
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    cout << "level_component_size = " << level_component_size << endl;
    vector<unsigned char *> byte_encoders;
    for(int i=0; i<num_level_component; i++){
        unsigned char * buffer = (unsigned char *) malloc(level_component_size);
        intra_level_components.push_back(buffer);
        byte_encoders.push_back(buffer);
        // set indicator
        bitplane_indicator.push_back(0);
    }    
    size_t buffer_index = byte_wise_direct_encoding_unrolled(data, n, level_exp, num_level_component, byte_encoders);
    for(int k=0; k<num_level_component; k++){
        encoded_sizes.push_back(buffer_index);
    }
    struct timespec start, end;
    bool use_rle = true;
    for(int k=1; k<num_level_component; k++){
        if((k > 25) || use_rle){
            // int err = clock_gettime(CLOCK_REALTIME, &start);        
            RunlengthEncoder rle;
            for(int i=0; i<buffer_index; i++){
                unsigned char datum = byte_encoders[k][i];
                for(int j=0; j<8; j++){
                    rle.encode(datum & 1);
                    datum >>= 1;            
                }
            }
            rle.flush();
            free(intra_level_components[k]);
            // change content of level components, encoded size and indicator
            intra_level_components[k] = rle.save();
            encoded_sizes[k] = rle.size();
            bitplane_indicator[k] = 1;
            if(rle.size() * 1.5 > level_component_size) use_rle = false;
            // cout << "RLE index = " << k << endl;  
            // err = clock_gettime(CLOCK_REALTIME, &end);
            // cout << "bitplane " << k << " runlength encoding time = " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000;
            // cout << "s, encoded size = " << rle.size() << endl;
        }
    }
    return intra_level_components;
}

template <class T>
T * progressive_hybrid_decoding(const vector<const unsigned char*>& level_components, size_t n, int level_exp, int num_level_component, const vector<unsigned char>& bitplane_indicator){
    T * level_data = (T *) malloc(n * sizeof(T));
    cout << "level element = " << n << endl;
    cout << "num_level_component = " << num_level_component << endl;
    vector<DecoderInterface*> decoders;
    for(int i=0; i<num_level_component; i++){
        switch(bitplane_indicator[i]){
            case 0:{
                decoders.push_back(new BitDecoder());
                break;
            }
            case 1:{
                decoders.push_back(new RunlengthDecoder());
                break;
            }
            default:{
                cerr << "Only direct encoding (indicator = 0) and RLE (indicator = 1) are supported\n";
                exit(0);
            }
        }
        decoders[i]->load(level_components[i]);
    }
    T * data_pos = level_data;
    for(int i=0; i<n; i++){
        // decode each bit of the data for each level component
        bool sign = decoders[0]->decode();
        unsigned int fp = 0;
        for(int j=1; j<num_level_component; j++){
            unsigned int current_bit = decoders[j]->decode();
            fp = (fp << 1) + current_bit;
        }
        long int fix_point = fp;
        if(sign) fix_point = -fix_point;
        *data_pos = ldexp((float)fix_point, - num_level_component + 1 + level_exp);
        data_pos ++;
    }
    for(int i=0; i<num_level_component; i++){
        delete decoders[i];
    }
    return level_data;
}

}
#endif