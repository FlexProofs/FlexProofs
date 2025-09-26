#include <vector>
#include <kzg.h>
#include <fc.h>

PI svc_batch_open(const FCpp pp, const KZGpp kzgpp, const Fp12 com, const std::vector<G1> a_vec,
    const std::vector<size_t> I, const std::vector<G1> y_vec);
bool svc_batch_verify(const FCpp pp, const KZGpp kzgpp, const Fp12 com, const std::vector<size_t> I,
                 const std::vector<G1> y_vec, const PI pi);