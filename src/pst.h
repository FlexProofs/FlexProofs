#ifndef PROOFOFLEARNING_PST_H
#define PROOFOFLEARNING_PST_H

#include <mcl/bls12_381.hpp>
#include <vector>
#include <mcl/ec.hpp>
using namespace mcl;
using namespace std;

#define F Fr

struct Comm{
    vector<G1> pp1;
    vector<G2> pp2;
    vector<G1> Proof;
    vector<F> secret;
    vector<vector<G1>> base;
    G1 G,C;
    G2 H;
    int bit_size;
};

void pst_commit(vector<Fr> poly, struct Comm &commitment);
void pst_open(vector<F> poly, vector<F> x, vector<vector<G1>> &pp, vector<G1> &Proof);
void pst_verify(vector<F> x, vector<G2> &pp2,G2 H, vector<G1> &Proof);
void pst_generatePP(struct Comm &commitment,int bit_size);
void pst_openall(vector<F> poly,vector<vector<G1>> &pp);

#endif //PROOFOFLEARNING_PST_H