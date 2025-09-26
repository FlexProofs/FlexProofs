#include <vector>
#include <string>
#include <cstring>
#include <stdexcept>
#include <algorithm>
#include <map>
#include <openssl/sha.h>
#include <unordered_map>
// #include <thread>
#include <mcl/ec.hpp>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <cstring>
#include <basic.h>
#include <cmath>
#include <iomanip>
using namespace std;

std::vector<uint8_t> Fr_serialize(const Fr& element) {
    size_t size = element.getByteSize();
    std::vector<uint8_t> buf(size);
    element.serialize(buf.data(), size,   IoSerialize);
    return buf;
}

std::vector<uint8_t> G1_serialize(const G1& element) {
    size_t size = element.getSerializedByteSize();
    std::vector<uint8_t> buf(size);
    element.serialize(buf.data(), size, IoSerialize);
    return buf;
}

std::vector<uint8_t> GT_serialize(const Fp12& element) {
    const size_t GT_SERIALIZE_SIZE = 384;
    std::vector<uint8_t> buf(GT_SERIALIZE_SIZE);
    element.serialize(buf.data(), GT_SERIALIZE_SIZE, IoSerialize);
    return buf;
}

Fr hash_to_x(const Fr& x_loop, const Fp12& left1, const G1& left2, const Fp12 right1, const G1 right2) {
    std::vector<uint8_t> ser_x_loop = Fr_serialize(x_loop);
    std::vector<uint8_t> left1s = GT_serialize(left1);
    std::vector<uint8_t> left2s = G1_serialize(left2);
    std::vector<uint8_t> right1s = GT_serialize(right1);
    std::vector<uint8_t> right2s = G1_serialize(right2);

    std::vector<uint8_t> input;
    input.insert(input.end(), ser_x_loop.begin(), ser_x_loop.end());
    input.insert(input.end(), left1s.begin(), left1s.end());
    input.insert(input.end(), left2s.begin(), left2s.end());
    input.insert(input.end(), right1s.begin(), right1s.end());
    input.insert(input.end(), right2s.begin(), right2s.end());

    // for (size_t i = 0; i < input.size(); ++i) {
    //     printf("%02x ", input[i]);
    // }
    // printf("\n");
    Fr res;
    res.setHashOf(input.data(), input.size());
    return res;
}

std::vector<Fr> hash_to_r(const Fp12 com, const std::vector<std::vector<Fr>> b_matrix, const std::vector<G1> y_vec) {
    std::vector<uint8_t> tmp;
    std::vector<uint8_t> ser_com;
    ser_com  = GT_serialize(com);
    tmp.insert(tmp.end(), ser_com.begin(), ser_com.end());

    for (size_t i = 0; i < b_matrix.size(); ++i) {
        const auto& row = b_matrix[i];
        for (size_t j = 0; j < row.size(); ++j) {
            const Fr& elem = row[j];
            if (elem == Fr(0)) continue;

            uint8_t i_bytes[sizeof(size_t)];
            memcpy(i_bytes, &i, sizeof(size_t));
            tmp.insert(tmp.end(), i_bytes, i_bytes + sizeof(size_t));

            uint8_t j_bytes[sizeof(size_t)];
            memcpy(j_bytes, &j, sizeof(size_t));
            tmp.insert(tmp.end(), j_bytes, j_bytes + sizeof(size_t));

            std::vector<uint8_t> elem_bytes = Fr_serialize(elem);
            tmp.insert(tmp.end(), elem_bytes.begin(), elem_bytes.end());
        }
    }

    for (const auto& g1 : y_vec) {
        std::vector<uint8_t> ser_g1;
        ser_g1 = G1_serialize(g1);
        tmp.insert(tmp.end(), ser_g1.begin(), ser_g1.end());
    }

    std::vector<uint8_t> digest;
    // （uint8_t* is unsigned char*）
    unsigned char* tmp_data = tmp.data();
    unsigned char* hash1 = SHA512(tmp_data, sizeof(tmp), nullptr);
    digest = std::vector<uint8_t>(hash1, hash1 + SHA512_DIGEST_LENGTH);

    std::vector<Fr> r_vec;
    r_vec.reserve(b_matrix.size());
    for (size_t i = 0; i < b_matrix.size(); ++i) {
        const auto& row = b_matrix[i];
        std::vector<uint8_t> b_vec_ser;

        for (size_t j = 0; j < row.size(); ++j) {
            const Fr& elem = row[j];
            if (elem == Fr(0)) continue;

            uint8_t i_bytes[sizeof(size_t)];
            memcpy(i_bytes, &i, sizeof(size_t));
            b_vec_ser.insert(b_vec_ser.end(), i_bytes, i_bytes + sizeof(size_t));

            uint8_t j_bytes[sizeof(size_t)];
            memcpy(j_bytes, &j, sizeof(size_t));
            b_vec_ser.insert(b_vec_ser.end(), j_bytes, j_bytes + sizeof(size_t));

            std::vector<uint8_t> elem_bytes = Fr_serialize(elem);
            b_vec_ser.insert(b_vec_ser.end(), elem_bytes.begin(), elem_bytes.end());
        }

        std::vector<uint8_t> input(b_vec_ser);
        input.insert(input.end(), digest.begin(), digest.end());
        Fr r;
        r.setHashOf(input.data(), input.size());
        r_vec.push_back(r);
    }

    return r_vec;
}

std::vector<Fr> hash_to_r_svc(Fp12 com, std::vector<size_t> I, std::vector<G1> y_vec) {
    std::vector<uint8_t> tmp(1024);
    std::vector<uint8_t> ser_com;

    ser_com = GT_serialize(com);
    tmp.insert(tmp.end(), ser_com.begin(), ser_com.end());
    Fr one(1);

    for (size_t i = 0; i < I.size(); ++i) {
        size_t j = I[i];

        uint8_t i_bytes[sizeof(size_t)];
        memcpy(i_bytes, &i, sizeof(size_t));
        tmp.insert(tmp.end(), i_bytes, i_bytes + sizeof(size_t));

        uint8_t j_bytes[sizeof(size_t)];
        memcpy(j_bytes, &j, sizeof(size_t));
        tmp.insert(tmp.end(), j_bytes, j_bytes + sizeof(size_t));

        std::vector<uint8_t> elem_bytes = Fr_serialize(one);
        tmp.insert(tmp.end(), elem_bytes.begin(), elem_bytes.end());

    }

    for (G1 g1 : y_vec) {
        std::vector<uint8_t> ser_g1;
        ser_g1 = G1_serialize(g1);
        tmp.insert(tmp.end(), ser_g1.begin(), ser_g1.end());
    }

    std::vector<uint8_t> digest;
    // （uint8_t* is unsigned char*）
    unsigned char* tmp_data = tmp.data();
    unsigned char* hash1 = SHA512(tmp_data, sizeof(tmp), nullptr);
    digest = std::vector<uint8_t>(hash1, hash1 + SHA512_DIGEST_LENGTH);

    std::vector<Fr> r_vec;
    r_vec.reserve(I.size());
    for (size_t i = 0; i < I.size(); ++i) {
        size_t j = I[i];
        std::vector<uint8_t> b_vec_ser;

        uint8_t i_bytes[sizeof(size_t)];
        memcpy(i_bytes, &i, sizeof(size_t));
        b_vec_ser.insert(b_vec_ser.end(), i_bytes, i_bytes + sizeof(size_t));

        uint8_t j_bytes[sizeof(size_t)];
        memcpy(j_bytes, &j, sizeof(size_t));
        b_vec_ser.insert(b_vec_ser.end(), j_bytes, j_bytes + sizeof(size_t));

        std::vector<uint8_t> elem_bytes = Fr_serialize(one);
        b_vec_ser.insert(b_vec_ser.end(), elem_bytes.begin(), elem_bytes.end());

        std::vector<uint8_t> input(b_vec_ser);
        input.insert(input.end(), digest.begin(), digest.end());
        Fr r;
        r.setHashOf(input.data(), input.size());
        r_vec.push_back(r);
    }
    return r_vec;
}

std::vector<uint8_t> to_fixed_length_binary_vec(size_t value, size_t length) {
    std::string binary_str;
    if (value == 0) {
        binary_str = "0";
    } else {
        while (value > 0) {
            binary_str = (value % 2 ? "1" : "0") + binary_str;
            value /= 2;
        }
    }
    if (binary_str.size() < length) {
        binary_str = std::string(length - binary_str.size(), '0') + binary_str;
    } else if (binary_str.size() > length) {
        binary_str = binary_str.substr(binary_str.size() - length);
    }

    std::vector<uint8_t> result;
    for (char c : binary_str) {
        result.push_back(c - '0');
    }
    return result;
}

size_t binary_vec_to_decimal(std::vector<uint8_t> bits) {
    size_t result = 0;
    for (uint8_t b : bits) {
        result = (result << 1) | b;
    }
    return result;
}

Fr find_root_of_unity(size_t log_n) {
    Fr root;
    root.setStr("5", 10);
    return root;
}

std::vector<Fr> polynomial_coefficients_from_x(const std::vector<Fr>& x_vec) {
    // for (size_t i = 0; i < x_vec.size(); ++i) {
    //     std::cout << "i, x_vec[i]===" << i << "," << x_vec[i] << std::endl;
    // }

    std::vector<Fr> coefficients = {Fr::one()};
    for (size_t i = 0; i < x_vec.size(); ++i) {
        size_t k = 1 << i;
        Fr x = x_vec[i];
        for (size_t j = 0; j < k; ++j) {
            Fr tmp;
            Fr::mul(tmp, x, coefficients[j]);
            coefficients.push_back(tmp);
        }
    }

    // for (size_t i = 0; i < coefficients.size(); ++i) {
    //     std::cout << "i, coefficients[i]===" << i << "," << coefficients[i] << std::endl;
    // }

    std::vector<Fr> result;
    result.reserve(coefficients.size() * 2 - 1);

    for (size_t i = 0; i < coefficients.size(); ++i) {
        result.push_back(coefficients[i]);
        if (i < coefficients.size() - 1) {
            result.push_back(Fr(0));
        }
    }
    // for (size_t i = 0; i < result.size(); ++i) {
    //     std::cout << "i, result[i]===" << i << "," << result[i] << std::endl;
    // }

    return result;
}

#define ROUNDS 80
int partitions;
std::vector<F> Common;
F K_MIMC;

void init_hash(){
    partitions = 4;

    Common.resize(ROUNDS);
    for(int i = 0; i < ROUNDS; i++){
        Common[i] = random();
    }
    K_MIMC = random();
}

int merkle_tree::size_after_padding;

void merkle_tree::merkle_tree_prover::MT_commit(vector<F> &leafs,int size,vector<vector<hDigest>> &hashes){
    vector<hDigest> leaf_hashes(size/2);
    hDigest data[2], ret;
    //data[0] = prev_hash;
    F element[2];
    char buff[32];
    //#pragma omp parallel for
    for(int i = 0; i < size/2; i++){
        int n = leafs[2*i].getStr(buff,32);
        leafs[2*i].serialize(buff,32);

        memcpy(&data[0], buff, 32);

        leafs[2*i+1].serialize(buff,32);
        memcpy(&data[1], buff, 32);

        my_hhash(data, &leaf_hashes[i]);
    }


    if(hashes.size() == 0){
        hashes.resize((int) ::log2(leaf_hashes.size())+1);
    }
    hashes[0] = leaf_hashes;
    merkle_tree::merkle_tree_prover::create_tree( leaf_hashes.size() , leaf_hashes, hashes, sizeof(hDigest), true);
}

hDigest merkle_tree::merkle_tree_prover::create_tree(int ele_num, vector<hDigest> &dst_ptr, vector<vector<hDigest>> &hashes, const int element_size = 256 / 8, bool alloc_required = false)
{
    //assert(element_size == sizeof(prime_field::field_element) * 2);

    int current_lvl_size = ele_num;

    int level = 1;
    current_lvl_size /= 2;
    hDigest data[2];
    int lvl = 1;

    while(current_lvl_size >= 1)
    {
        if(lvl >= 0){
            hashes[lvl].resize(current_lvl_size);
        }
        #pragma omp parallel for
        for(int i = 0; i < current_lvl_size; ++i)
        {
            //printf("%d,%d\n",start_ptr + current_lvl_size + i * 2,start_ptr + current_lvl_size + i * 2+1);
            data[0] = dst_ptr[ i * 2];
            data[1] = dst_ptr[ i * 2 + 1];
            my_hhash(data, &dst_ptr[i]);
            if(lvl >= 0){
                hashes[lvl][i] = dst_ptr[i];
            }
        }
        lvl++;
        current_lvl_size /= 2;
    }
    //dst = dst_ptr.data();
    return dst_ptr[0];
}


vector<hDigest> merkle_tree::merkle_tree_prover::open_tree(vector<vector<hDigest>> &MT_hashes, vector<size_t> c,  int collumns){
    int pos = (c[1]/4)*collumns + c[0];
    vector<hDigest> path;
    for(int i = 0; i < MT_hashes.size(); i++){
        path.push_back(MT_hashes[i][2*(pos/2) + (1-(pos%2))]);
        pos = pos/2;
    }
    return path;
}

bool merkle_tree::merkle_tree_verifier::verify_claim(hDigest root_hhash, const hDigest* tree, hDigest leaf_hash, int pos_element_arr, int N)
{
    //check N is power of 2
    // assert((N & (-N)) == N);
    //
    // int pos_element = pos_element_arr + N;
    // __hhash_digest data[2];
    // while(pos_element != 1)
    // {
    //     data[pos_element & 1] = leaf_hash;
    //     data[(pos_element & 1) ^ 1] = tree[pos_element ^ 1];
    //     my_hhash(data, &leaf_hash);
    //     pos_element /= 2;
    // }
    // return equals(root_hhash, leaf_hash);

    // Check if N is a power of 2
    assert((N & (-N)) == N);
    // Check if index is within valid range
    assert(pos_element_arr >= 0 && pos_element_arr < N && "Index out of valid range");

    // Start calculation from leaf hash layer position (each hash corresponds to 2 leaves)
    int pos_element = pos_element_arr / 2;
    // Calculate number of iterations needed (tree height)
    int tree_height = static_cast<int>(log2(N)) - 1;

    hDigest data[2];
    hDigest current_hash = leaf_hash;  // Use temporary variable to store current hash

    // Iterate according to tree height instead of fixing to pos_element=1
    for (int i = 0; i < tree_height; ++i)
    {
        // Correctly handle left and right node order
        if (pos_element % 2 == 0) {
            // Current node is left child, proof is right child
            data[0] = current_hash;
            data[1] = tree[i];
        } else {
            // Current node is right child, proof is left child
            data[0] = tree[i];
            data[1] = current_hash;
        }

        my_hhash(data, &current_hash);
        pos_element /= 2;  // Move up to parent node
    }

    return equals(root_hhash, current_hash);
}