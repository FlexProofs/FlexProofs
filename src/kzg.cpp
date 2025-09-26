#include <math.h>

// #include <thread>
#include <mcl/ec.hpp>
#include <openssl/sha.h>
#include <kzg.h>

#include "basic.h"


// set KZG's public parameter
KZGpp kzg_setup(size_t n, const Fr alpha) {
    KZGpp kzgpp;

    kzgpp.g2_vec.reserve(n);
    Fr alpha_power = Fr(1);
    mcl::bn::hashAndMapToG2(kzgpp.g2_gen, "g2_gen");// get a generate in G2

    for (size_t i = 0; i < n; ++i) {
        G2 tmp;
        G2::mul(tmp, kzgpp.g2_gen, alpha_power); // tmp = alpha_power · g2_gen
        kzgpp.g2_vec.push_back(tmp);
        // cout << i << " " << tmp << endl;
        Fr::mul(alpha_power, alpha_power, alpha);
        // alpha_power *= alpha;
    }

    // g1_alpha =  g1_gen * alpha
    mcl::bn::hashAndMapToG1(kzgpp.g1_gen, "g1_gen");
    G1::mul(kzgpp.g1_alpha, kzgpp.g1_gen, alpha);

    return kzgpp;
}

G2 kzg_commit(KZGpp pp, std::vector<Fr> poly_values) {
    G2 commitment;
    G2::mulVec(commitment, pp.g2_vec.data(), poly_values.data(), poly_values.size());

    return commitment;
}

G2 kzg_open(KZGpp pp, std::vector<Fr> poly, Fr z) {
    Fr y = kzg_evaluate(poly, z);

    // divide (x-z)
    std::vector<Fr> quotient = compute_quotient_coeffs(poly, z, y);

    G2 pi;
    std::vector<G2> sub_pp(pp.g2_vec.begin(), pp.g2_vec.begin() + quotient.size());
    G2::mulVec(pi, sub_pp.data(), quotient.data(), quotient.size());

    return pi;
}

bool kzg_verify(
    KZGpp pp,
    G2 commitment,
    Fr point,
    Fr value,
    G2 pi
) {
    // compute left：e(G1, commitment - G2^value)
    G2 tmp;
    G2::mul(tmp, pp.g2_gen, value);
    G2 com_sub_rx = commitment - tmp;

    Fp12 lhs;
    mcl::bn::pairing(lhs, pp.g1_gen, com_sub_rx);

    // compute right：e(alpha*G1 - point*G1, pi) = e((alpha - point)*G1, pi)
    G1 tmp2;
    G1::mul(tmp2, pp.g1_gen, point);
    G1 gx_sub_point = pp.g1_alpha - tmp2;

    Fp12 rhs;
    pairing(rhs, gx_sub_point, pi);

    return lhs == rhs;
}

// evaluation of polynomial
Fr kzg_evaluate(std::vector<Fr> coefficients, Fr point) {
    Fr result = Fr(0);
    Fr power = Fr(1);

    for (const Fr coeff : coefficients) {
        Fr term = coeff * power;
        result += term;
        power *= point;
    }

    return result;
}

// compute quotient
std::vector<Fr> compute_quotient_coeffs(std::vector<Fr> coeffs,  Fr z, Fr y) {
    // f(x) - f(z)
    std::vector<Fr> a = coeffs;
    a[0] -= y;

    size_t degree = a.size() - 1;
    std::vector<Fr> q_coeffs(degree, Fr(0));  // degree-

    // (f(x) - f(z)) / (x - z)
    q_coeffs[degree - 1] = a[degree];
    for (int i = degree - 1; i >= 1; --i) {
        Fr tmp = z * q_coeffs[i];
        q_coeffs[i - 1] = a[i] + tmp;
    }

    // check whether remainder is zero
    Fr remainder = a[0] + z * q_coeffs[0];
    assert(remainder.isZero());

    return q_coeffs;
}


