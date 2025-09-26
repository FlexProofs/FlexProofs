#ifndef KZG_HPP
#define KZG_HPP

#include <mcl/bls12_381.hpp>
#include <vector>
#include <string>
#include <cassert>

using namespace mcl::bls12;

struct KZGpp {
    G1 g1_gen;
    G2 g2_gen;
    G1 g1_alpha;
    std::vector<G2> g2_vec;
};

// set KZG's public parameter
KZGpp kzg_setup(size_t n, const Fr alpha);

G2 kzg_commit(KZGpp pp, std::vector<Fr> poly_values);

G2 kzg_open(KZGpp pp, std::vector<Fr> poly, Fr z);

bool kzg_verify(
    KZGpp pp,
    G2 commitment,
    Fr point,
    Fr value,
    G2 pi
);

// evaluation of polynomial
Fr kzg_evaluate(std::vector<Fr> coefficients, Fr point);

// compute quotient
std::vector<Fr> compute_quotient_coeffs(std::vector<Fr> coeffs,  Fr z, Fr y);




#endif // KZG_HPP
