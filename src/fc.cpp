#include <vector>
#include <cmath>
#include <fc.h>

using namespace std;

FCpp fc_gen(const size_t n, const Fr beta) {
    vector<G2> v_vec;
    v_vec.reserve(n);

    Fr beta_square;
    Fr::sqr(beta_square, beta);

    Fr beta_power = Fr::one();
    G2 g2_gen2;
    mcl::bn::hashAndMapToG2(g2_gen2, "g2_gen");

    for (size_t i = 0; i < n; ++i) {
        G2 tmp2;
        G2::mul(tmp2, g2_gen2, beta_power);
        // cout << i << " " << tmp2 << endl;
        v_vec.push_back(tmp2);

        Fr::mul(beta_power, beta_power, beta_square);
        // beta_power *= beta_square;
    }

    return {v_vec};
}

Fp12 fc_com(const FCpp pp, const vector<G1> a_vec) {
    Fp12 comTmp, com;
    mcl::bn::millerLoopVec(comTmp, a_vec.data(), pp.v_vec.data(),a_vec.size());
    mcl::bn::finalExp(com, comTmp);
    return com;
}

PI fc_batch_open(const FCpp pp, const KZGpp kzgpp, const Fp12 com, vector<G1> a_vec,
    const vector<vector<Fr>> b_matrix, vector<G1> y_vec) {
    size_t n = pp.v_vec.size();
    vector<Fr> r_vec = hash_to_r(com, b_matrix, y_vec);

    vector<Fr> b_vec;
    b_vec.reserve(n);
    for (size_t i = 0; i < n; ++i) { //column index of b_matrix
        Fr tmpi(0);
        for (size_t j = 0; j < y_vec.size(); ++j) { // row index of b_matrix
            tmpi += r_vec[j] * b_matrix[j][i];
        }
        b_vec.push_back(tmpi);
    }

    size_t l = static_cast<size_t>(log2(n));

    vector<Fr> x_inverse_vec;
    x_inverse_vec.reserve(l);
    vector<pair<Fp12, G1>> L_vec;
    L_vec.reserve(l);
    vector<pair<Fp12, G1>> R_vec;
    R_vec.reserve(l);
    vector<G1> a_loop(a_vec);
    vector<G2> v_loop(pp.v_vec);
    vector<Fr> b_loop(b_vec);

    Fr x_loop(0);
    size_t mid = n >> 1;

    for (size_t j = 1; j <= l; ++j) {
        vector<G1> a_left(a_loop.begin(), a_loop.begin() + mid);
        vector<G1> a_right(a_loop.begin() + mid, a_loop.end());
        vector<G2> v_left(v_loop.begin(), v_loop.begin() + mid);
        vector<G2> v_right(v_loop.begin() + mid, v_loop.end());
        vector<Fr> b_left(b_loop.begin(), b_loop.begin() + mid);
        vector<Fr> b_right(b_loop.begin() + mid, b_loop.end());

        Fp12 tmpleft1, left1;
        mcl::bn::millerLoopVec(tmpleft1, a_right.data(), v_left.data(), mid);
        mcl::bn::finalExp(left1, tmpleft1);
        G1 left2;
        G1::mulVec(left2, a_right.data(), b_left.data(), mid);
        L_vec.emplace_back(left1, left2);

        Fp12 tmpright1, right1;
        mcl::bn::millerLoopVec(tmpright1, a_left.data(), v_right.data(),mid);
        mcl::bn::finalExp(right1, tmpright1);
        G1 right2;
        G1::mulVec(right2, a_left.data(), b_right.data(), mid);
        R_vec.emplace_back(right1, right2);

        x_loop = hash_to_x(x_loop, left1, left2, right1, right2);
        Fr x_loop_inverse;
        Fr::inv(x_loop_inverse, x_loop);
        x_inverse_vec.push_back(x_loop_inverse);

        a_loop.clear();
        a_loop.resize(mid);
        v_loop.clear();
        v_loop.resize(mid);
        b_loop.clear();
        b_loop.resize(mid);

        for (size_t k = 0; k < mid; ++k) {
            G1 tmp;
            G1::mul(tmp, a_right[k], x_loop);
            G1::add(a_loop[k], a_left[k], tmp);

            G2 tmp2;
            G2::mul(tmp2, v_right[k], x_loop_inverse);
            G2::add(v_loop[k], v_left[k], tmp2);

            Fr tmp3 ;
            Fr::mul(tmp3,  b_right[k], x_loop_inverse);
            Fr::add(b_loop[k], b_left[k], tmp3);
        }
        mid >>= 1;
    }

    G1 finalA = a_loop[0];
    G2 finalv = v_loop[0];

    reverse(x_inverse_vec.begin(), x_inverse_vec.end());
    vector<Fr> coeffs = polynomial_coefficients_from_x(x_inverse_vec);
    // G2  comTemp= kzg_commit(kzgpp, coeffs);
    // cout << "comTemp === " << comTemp << endl;
    // cout << "finalv === " << finalv << endl;

    string text = "KZG point";
    Fr point;
    point.setHashOf(text.data(), text.size());
    Fr value = kzg_evaluate(coeffs, point);
    G2 finalv_proof = kzg_open(kzgpp, coeffs, point);
    // bool kzg_result = kzg_verify(kzgpp, comTemp, point, value, finalv_proof);
    // cout << "kzg_result in prove=== " << kzg_result << endl;

    return {L_vec, R_vec, finalA, finalv, finalv_proof};
}

bool fc_batch_verify(const FCpp pp, const KZGpp kzgpp, const Fp12 com, const vector<vector<Fr>> b_matrix,
                  vector<G1> y_vec, const PI pi) {
    size_t n = pp.v_vec.size();
    vector<Fr> r_vec = hash_to_r(com, b_matrix, y_vec);

    vector<Fr> b_vec;
    b_vec.reserve(n);
    for (size_t i = 0; i < n; ++i) { //column index of b_matrix
        Fr tmpi(0);
        for (size_t j = 0; j < y_vec.size(); ++j) { // row index of b_matrix
            tmpi += r_vec[j] * b_matrix[j][i];
        }
        b_vec.push_back(tmpi);
    }

    G1 y;
    G1::mulVec(y, y_vec.data(), r_vec.data(), y_vec.size());

    const size_t l = static_cast<size_t>(log2(n));

    Fp12 c_loop1 = com;
    G1 c_loop2 = y;
    Fr x_loop(0);
    vector<Fr> x_inverse_vec;
    x_inverse_vec.reserve(l);
    vector<Fr> b_loop(b_vec);
    size_t mid = n >> 1;

    for (size_t j = 1; j <= l; ++j) {
        vector<Fr> b_left(b_loop.begin(), b_loop.begin() + mid);
        vector<Fr> b_right(b_loop.begin() + mid, b_loop.end());
        auto& [left1, left2] = pi.L_vec[j-1];
        auto& [right1, right2] = pi.R_vec[j-1];

        x_loop = hash_to_x(x_loop, left1, left2, right1, right2);

        Fr x_loop_inverse;
        Fr::inv(x_loop_inverse, x_loop);
        x_inverse_vec.push_back(x_loop_inverse);

        Fp12 left1_exp;
        Fp12::pow(left1_exp, left1, x_loop);
        Fp12 right1_exp;
        Fp12::pow(right1_exp, right1, x_loop_inverse);

        Fp12::mul(c_loop1, c_loop1, left1_exp);
        Fp12::mul(c_loop1, c_loop1, right1_exp);
        // c_loop1 *= left1_exp;
        // c_loop1 *= right1_exp;

        // c_loop2 += left2 * x_loop;
        G1 left2_exp;
        G1::mul(left2_exp, left2, x_loop);
        G1::add(c_loop2, c_loop2, left2_exp);
        // c_loop2 += right2 * x_loop_inverse;
        G1 right2_exp;
        G1::mul(right2_exp, right2, x_loop_inverse);
        G1::add(c_loop2, c_loop2, right2_exp);

        b_loop.clear();
        b_loop.insert(b_loop.end(), b_left.begin(), b_left.end());

        for (size_t k = 0; k < mid; ++k) {
            Fr tmp3 = b_right[k] * x_loop_inverse;
            b_loop[k] += tmp3;
        }
        mid = mid >> 1;
    }

    Fr final_b = b_loop[0];

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
    mcl::bn::pairing(pairing, pi.finalA, pi.finalv);
    bool r1 = (c_loop1 == pairing);

    G1 finalA_mul = pi.finalA * final_b;
    bool r2 = (c_loop2 == finalA_mul);

    // cout << "kzg_result === " << kzg_result << endl;
    // cout << "r1 === " << r1 << endl;
    // cout << "r2 === " << r2 << endl;

    return kzg_result && r1 && r2;
}