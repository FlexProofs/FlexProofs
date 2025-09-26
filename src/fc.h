#include <vector>
#include <cmath>
#include <stdexcept>
#include <kzg.h>
#include <basic.h>

FCpp fc_gen(const size_t n);
Fp12 fc_com(const FCpp pp, const std::vector<G1> a_vec);

FCpp fc_gen(const size_t n, const Fr beta);

Fp12 fc_com(const FCpp pp, const std::vector<G1> a_vec);

PI fc_batch_open(const FCpp pp, const KZGpp kzgpp, const Fp12 com, std::vector<G1> a_vec,
    const std::vector<std::vector<Fr>> b_matrix, std::vector<G1> y_vec);
bool fc_batch_verify(const FCpp pp, const KZGpp kzgpp, const Fp12 com, const std::vector<std::vector<Fr>> b_matrix,
        std::vector<G1> y_vec, const PI pi);