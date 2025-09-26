#include <mcl/bls12_381.hpp>
#include <vector>
#include <string>
#include <cstring>
#include <mcl/ec.hpp>
#include <cstdint>
using namespace mcl;

#define F Fr

struct FCpp {
    std::vector<G2> v_vec;
};

struct PI {
    std::vector<std::pair<Fp12, G1>> L_vec;
    std::vector<std::pair<Fp12, G1>> R_vec;
    G1 finalA;
    G2 finalv;
    G2 finalv_proof;
};

std::string fr_repr_to_string(const Fr& fr_repr);

std::vector<uint8_t> Fr_serialize(const Fr& element);

std::vector<uint8_t> G1_serialize(const G1& element);

std::vector<uint8_t> GT_serialize(const Fp12& element);

Fr hash_to_x(const Fr& x_loop, const Fp12& left1, const G1& left2, const Fp12 right1, const G1 right2);

std::vector<Fr> hash_to_r(const Fp12 com, const std::vector<std::vector<Fr>> b_matrix, const std::vector<G1> y_vec);

std::vector<Fr> hash_to_r_svc(Fp12 com, std::vector<size_t> I, std::vector<G1> y_vec);

std::vector<uint8_t> to_fixed_length_binary_vec(size_t value, size_t length);

size_t binary_vec_to_decimal(std::vector<uint8_t> bits);

Fr find_root_of_unity(size_t log_n);

std::vector<Fr> polynomial_coefficients_from_x(const std::vector<Fr>& x_vec);

#include <immintrin.h>
#include <cassert>

#define USESHA3
#ifdef USESHA3
extern "C"{
#include "../lib/libXKCP.a.headers/SimpleFIPS202.h"
}
#endif
#include <cstring>
inline void my_hhash(const void* src, void* dst)
{
#ifdef USESHA3
    SHA3_256((unsigned char*)dst, (const unsigned char*)src, 64);
#else
    //memset(dst, 0, sizeof(__hhash_digest));
    sha256_update_shani((const unsigned char*)src, 64, (unsigned char*)dst);
#endif
}
class hDigest
{
public:
    __m128i h0, h1;
};

inline bool equals(const hDigest &a, const hDigest &b)
{
    __m128i v0 = _mm_xor_si128(a.h0, b.h0);
    __m128i v1 = _mm_xor_si128(a.h1, b.h1);
    return _mm_test_all_zeros(v0, v0) && _mm_test_all_zeros(v1, v1);
}

namespace merkle_tree
{
    extern int size_after_padding;

    namespace merkle_tree_prover
    {
        //Merkle tree functions used by the prover
        //void create_tree(void* data_source_array, int lengh, void* &target_tree_array, const int single_element_size = 256/8)
        hDigest create_tree(int ele_num, std::vector<hDigest> &dst_ptr, std::vector<std::vector<hDigest>> &hashes, const int element_size, bool alloc_required);
        void MT_commit(std::vector<F> &leafs,int size,std::vector<std::vector<hDigest>> &hashes);
        std::vector<hDigest> open_tree(std::vector<std::vector<hDigest>> &MT_hashes, std::vector<size_t> c,  int collumns);

    }

    namespace merkle_tree_verifier
    {
        //Merkle tree functions used by the verifier
        bool verify_claim(hDigest root_hhash, const hDigest* tree, hDigest hhash_element, int pos_element, int N);
    }
}

void init_hash();