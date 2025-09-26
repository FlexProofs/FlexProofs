#include <vector>
#include <cmath>
#include <svc.h>

using namespace std;

PI svc_batch_open(const FCpp pp, const KZGpp kzgpp, const Fp12 com, const vector<G1> a_vec,
    const vector<size_t> I, const ::vector<G1> y_vec) {
    size_t n = pp.v_vec.size();
    vector<Fr> r_vec = hash_to_r_svc(com, I, y_vec);

    size_t l = static_cast<size_t>(log2(n));
    size_t I_size = I.size();
    vector<vector<uint8_t>> I_binary_vec;
    I_binary_vec.reserve(I_size);
    vector<Fr> final_ur_i(I_size, Fr::one());
    vector<Fr> b_loop(n, Fr(0));

    for (size_t index = 0; index < I_size; ++index) {
        I_binary_vec.push_back(to_fixed_length_binary_vec(I[index], l));
        final_ur_i[index]=r_vec[I[index]];
        b_loop[I[index]] = r_vec[I[index]];
    }

    vector<Fr> x_inverse_vec;
    x_inverse_vec.reserve(l);
    vector<pair<Fp12, G1>> L_vec;
    L_vec.reserve(l);
    vector<pair<Fp12, G1>> R_vec;
    R_vec.reserve(l);
    vector<G1> a_loop(a_vec);
    vector<G2> v_loop(pp.v_vec);

    Fr x_loop(0);
    Fr final_b_acc = Fr::one();
    size_t mid = n >> 1;

    for (size_t j = 1; j <= l; ++j) {
        vector<G1> a_left(a_loop.begin(), a_loop.begin() + mid);
        vector<G1> a_right(a_loop.begin() + mid, a_loop.end());
        vector<G2> v_left(v_loop.begin(), v_loop.begin() + mid);
        vector<G2> v_right(v_loop.begin() + mid, v_loop.end());
        vector<Fr> b_left(b_loop.begin(), b_loop.begin() + mid);
        vector<Fr> b_right(b_loop.begin() + mid, b_loop.end());

        Fp12 left1Temp, left1;
        mcl::bn::millerLoopVec(left1Temp, a_right.data(), v_left.data(), mid);
        mcl::bn::finalExp(left1, left1Temp);
        G1 left2;
        G1::mulVec(left2, a_right.data(), b_left.data(), a_right.size());
        L_vec.emplace_back(left1, left2);

        Fp12 right1Temp, right1;
        mcl::bn::millerLoopVec(right1Temp, a_left.data(), v_right.data(),mid);
        mcl::bn::finalExp(right1, right1Temp);
        G1 right2;
        G1::mulVec(right2, a_left.data(), b_right.data(), a_left.size());
        R_vec.emplace_back(right1, right2);

        x_loop = hash_to_x(x_loop, left1, left2, right1, right2);
        Fr x_loop_inverse;
        x_loop.inv(x_loop_inverse, x_loop);
        x_inverse_vec.push_back(x_loop_inverse);

        a_loop.clear();
        a_loop.resize(mid);
        v_loop.clear();
        v_loop.resize(mid);
        b_loop.clear();
        b_loop.resize(mid, Fr(0));

        for (size_t k = 0; k < mid; ++k) {
            G1 tmp;
            G1::mul(tmp, a_right[k], x_loop);
            G1::add(a_loop[k], a_left[k], tmp);

            G2 tmp2;
            G2::mul(tmp2, v_right[k], x_loop_inverse);
            G2::add(v_loop[k], v_left[k], tmp2);
        }

        for (size_t index = 0; index < I_size; ++index) {
            if (I_binary_vec[index][j-1] == 1) {
                final_ur_i[index] *= x_loop_inverse;
            }

            vector<uint8_t> tenForm(I_binary_vec[index].begin() + j, I_binary_vec[index].end());
            size_t idx = binary_vec_to_decimal(tenForm);
            Fr tmp_fr;
            tmp_fr=b_loop[idx];
            tmp_fr += final_ur_i[index];
            b_loop[idx] = tmp_fr;
        }

        mid >>= 1;
    }

    G1 finalA = a_loop[0];
    G2 finalv = v_loop[0];

    reverse(x_inverse_vec.begin(), x_inverse_vec.end());
    vector<Fr> coeffs = polynomial_coefficients_from_x(x_inverse_vec);

    string text = "KZG point";
    Fr point;
    point.setHashOf(text.data(), text.size());
    Fr value = kzg_evaluate(coeffs, point);
    G2 finalv_proof = kzg_open(kzgpp, coeffs, point);

    return {L_vec, R_vec, finalA, finalv, finalv_proof};
}

bool svc_batch_verify(const FCpp pp, const KZGpp kzgpp, const Fp12 com, const vector<size_t> I,
                 const vector<G1> y_vec, const PI pi) {
    size_t n = pp.v_vec.size();
    vector<Fr> r_vec = hash_to_r_svc(com, I, y_vec);

    G1 y;
    G1::mulVec(y, const_cast<G1*>(y_vec.data()), r_vec.data(), y_vec.size());

    const size_t l = static_cast<size_t>(log2(n));
    const size_t I_size = I.size();
    vector<vector<uint8_t>> I_binary_vec;
    I_binary_vec.reserve(I_size);
    vector<Fr> final_ur_i(I_size, Fr::one());

    for (size_t index = 0; index < I_size; ++index) {
        I_binary_vec.push_back(to_fixed_length_binary_vec(I[index], l));
        final_ur_i[index]=r_vec[I[index]];
    }

    Fp12 c_loop1 = com;
    G1 c_loop2 = y;
    Fr x_loop(0);
    vector<Fr> x_inverse_vec;
    x_inverse_vec.reserve(l);

    for (size_t j = 1; j <= l; ++j) {
        auto& [left1, left2] = pi.L_vec[j-1];
        auto& [right1, right2] = pi.R_vec[j-1];
        x_loop = hash_to_x(x_loop, left1, left2, right1, right2);

        Fr x_loop_inverse;
        x_loop.inv(x_loop_inverse, x_loop);
        x_inverse_vec.push_back(x_loop_inverse);

        Fp12 left1_exp;
        Fp12::pow(left1_exp, left1, x_loop);
        Fp12 right1_exp;
        Fp12::pow(right1_exp, right1, x_loop_inverse);

        // c_loop1 *= left1_exp;
        // c_loop1 *= right1_exp;
        Fp12::mul(c_loop1, c_loop1, left1_exp);
        Fp12::mul(c_loop1, c_loop1, right1_exp);

        // c_loop2 += left2 * x_loop;
        G1 left2_exp;
        G1::mul(left2_exp, left2, x_loop);
        G1::add(c_loop2, c_loop2, left2_exp);
        // c_loop2 += right2 * x_loop_inverse;
        G1 right2_exp;
        G1::mul(right2_exp, right2, x_loop_inverse);
        G1::add(c_loop2, c_loop2, right2_exp);

        for (size_t index = 0; index < I_size; ++index) {
            if (I_binary_vec[index][j-1] == 1) {
                final_ur_i[index] *= x_loop_inverse;
            }
        }
    }

    Fr final_b(0);
    for (size_t index = 0; index < I_size; ++index) {
        final_b += final_ur_i[index];
    }

    string text = "KZG point";
    Fr point;
    point.setHashOf(text.data(), text.size());
    Fr value = Fr::one();
    Fr term = point;

    for (size_t i = 0; i < l; ++i) {
        Fr::sqr(term, term);
        Fr coffi = x_inverse_vec[l-i-1];
        coffi *= term;
        coffi += Fr::one();
        value *= coffi;
    }

    bool kzg_result = kzg_verify(kzgpp, pi.finalv, point, value, pi.finalv_proof);

    Fp12 pairing;
    mcl::bls12::pairing(pairing, pi.finalA, pi.finalv);
    bool r1 = (c_loop1 == pairing);

    G1 finalA_mul = pi.finalA * final_b;
    bool r2 = (c_loop2 == finalA_mul);

    return kzg_result && r1 && r2;
}