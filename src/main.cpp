#include <stdio.h>
#include <mcl/bls12_381.hpp>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <chrono>
#include <svc.h>
#include <cmath>
#include <pst.h>

using namespace std::chrono;
using namespace mcl::bn;

extern int partitions;
using namespace std;
double proving_time = 0.0;
double verification_time;
int comm_size = 0;
int ps = 0;

double vt = 0.0;

// void test_Hydraproof(int logN) {
// 	cout << "========================Hydraproofs=====================" << endl;
// 	int N = 1ULL<<logN;
// 	cout << "N = " << N << endl;
// 	vector<F> rand(10);
// 	for(int i = 0; i < 10; i++){
// 		rand[i].setByCSPRNG();
// 	}
// 	vector<F> data(N);
// 	data[0] = rand[9];
// 	for(int i = 1; i < N; i++){
// 		data[i] = data[i-1]*rand[0];
// 	}
//
// 	VC_pp pp;
// 	MIPP_Comm pp_pc,pp_input;
//
// 	generate_pp_VC((int)log2(N),pp);
// 	pp_pc.pp1 = pp.pp1;
// 	pp_pc.pp2 = pp.pp2;
// 	pp_pc.precompute = false;
// 	pp_pc.bits = false;
// 	pp_pc.base = pp.base;
// 	pp_pc.G = pp.G;
// 	pp_pc.H = pp.H;
// 	pp_pc.u = pp.u;
// 	pp_pc.w = pp.w;
//
// 	pp_pc.bit_size = pp.bit_size;
// 	pp_pc.commit_bit_size = pp.commit_bit_size;
// 	vector<G1> partial_commitments;
// 	GT C;
// 	vector<int>v_int;
//
// 	auto start = std::chrono::high_resolution_clock::now();
// 	VC_Commit(data,partial_commitments,C,pp);
// 	auto stop = std::chrono::high_resolution_clock::now();
// 	auto duration = duration_cast<milliseconds>(stop - start);
// 	cout << "com generated: " << duration.count() << "ms" << endl;
//
// 	int ps = 0;
// 	vector<F> poly = data;
// 	start = std::chrono::high_resolution_clock::now();
// 	int n = (int)log2(poly.size());
//     int segments,segment_size;
//
//     if(n%2==0){
//         segments = 1ULL<<(n/2);
//         segment_size = 1ULL<<(n/2);
//     }else{
//         segments = 1ULL<<((n-1)/2);
//         segment_size = 1ULL<<((n+1)/2);
//         if(segments*segment_size!= poly.size()){
//             printf("Error\n");
//             exit(-1);
//         }
//     }
//     vector<F> r;
//     vector<F> r1((int)log2(segment_size)),r2((int)log2(segments));
//     // merkle commit the partital commitments
//     //printf("Phase 1\n");
//     // Phase 1:
//     // Get randomness r1;
// 	vector<Fr> MTr1_leafs;
// 	vector<vector<uint8_t>> partial_commitments_ser;
// 	for (int i=0; i<segments; i++) {
// 		partial_commitments_ser.push_back(G1_serialize(partial_commitments[i]));
// 	}
// 	for(int i = 0; i < segments*segment_size; i++) {
// 		Fr res;
// 		vector<uint8_t> leafData = Fr_serialize(poly[i]);
// 		int comIdx = i/segment_size;
// 		leafData.insert(leafData.end(), partial_commitments_ser[comIdx].begin(), partial_commitments_ser[comIdx].end());
// 		res.setHashOf(leafData.data(), leafData.size());
// 		MTr1_leafs.push_back(res);
// 	}
//    	vector<vector<hDigest>> MTr1;
// 	merkle_tree::merkle_tree_prover::MT_commit(MTr1_leafs,segments*segment_size,MTr1);
// 	ps += ((int)log2(segments*segment_size))*32;
//     for(int i = 0; i < r1.size(); i++){
//         r1[i].setByCSPRNG();
//         r.push_back(r1[i]);
//     }
//
//     vector<F> poly_compressed(segments);// polynomials to be compressed
//     vector<F> buff(segment_size);
//     for(int i = 0; i < segments; i++){
//         for(int j = 0; j < segment_size; j++){
//             buff[j] = poly[i*segment_size+j];
//         }
//         poly_compressed[i] = evaluate_vector(buff,r1);
//         //MY CODE
//         vector<G1> buff_proof;
//         open(buff,r1,pp.base_compressed,buff_proof);
//         //MY CODE
//     }
// 	ps += ((int)log2(buff.size()))*32+32;
//     // Commit+ openall
//     G1 C_compressed;
//     G1::mulVec(C_compressed, pp.pp1_compressed.data(), poly_compressed.data(), poly_compressed.size());
//
//     hyperproofs_openall(poly_compressed,pp.base_compressed);
//     ps += (int)log2((poly_compressed.size()))*32+2*32;
//
//     // Sample a second point
// 	vector<Fr> MTr2_leafs;
// 	std::vector<uint8_t> ser2 = G1_serialize(C_compressed);
// 	Fr res;
// 	res.setHashOf(ser2.data(), ser2.size());
//     for(int i = 0; i < r2.size(); i++){
//         r2[i].setByCSPRNG();
//         r.push_back(r2[i]);
//     }
//
//     vector<G1> proof_compressed;
//     open(poly_compressed,r2,pp.base_compressed,proof_compressed);
//
//     ps += (proof_compressed.size()+1)*32;
//     // verify_proof(r2,pp.pp2_compressed,pp.H,proof_compressed);
//     // verify_proof(r2,pp.pp2_compressed,pp.H,proof_compressed);
//     // if(application == 3){
//
//         pp_pc.Commitments = partial_commitments;
//         MIPP_open(poly,r,pp_pc);
// 	ps += (pp_pc.Proof.size()+1)*32;
//     // }
//     printf("Phase 2\n");
//
//     vector<vector<vector<hDigest>>> MT(r2.size());
//     // vector<vector<vector<__hhash_digest>>> MT2(r2.size());
//     // vector<vector<vector<__hhash_digest>>> MT3(r2.size());
//     // Phase 2:
//     for(int k = 0; k < r2.size(); k++){
//         // Merkle commit the following:
//         // 1) The current polynomial
//         // 2) the partial evaluations
//         // 3) the partial commitmnets
//         // From 1+2+3 compute the new randomness
//     	vector<Fr> MT_leafs;
//     	vector<vector<uint8_t>> partial_commitments_ser;
//     	for (int i=0; i<segments; i++) {
//     		partial_commitments_ser.push_back(Fr_serialize(poly_compressed[i]));
//     	}
//     	for(int i = 0; i < segments*segment_size; i++) {
//     		Fr res;
//     		std::vector<uint8_t> leafData = Fr_serialize(poly[i]);
//     		int comIdx = i/segment_size;
//     		leafData.insert(leafData.end(), partial_commitments_ser[comIdx].begin(), partial_commitments_ser[comIdx].end());
//     		leafData.insert(leafData.end(), partial_commitments_ser[comIdx].begin(), partial_commitments_ser[comIdx].end());
//     		res.setHashOf(leafData.data(), leafData.size());
//     		MT_leafs.push_back(res);
//     	}
//         merkle_tree::merkle_tree_prover::MT_commit(MT_leafs,segments*segment_size,MT[k]);
//         // merkle_tree::merkle_tree_prover::MT_commit(poly_compressed,segments,MT2[k]);
//         //
//         // merkle_tree::merkle_tree_prover::MT_commit(MTr1_leafs,segments,MT3[k]);
//
//         ps += ((int)log2(segments*segment_size))*32;// + 2*(int)log2(segments))*32;
//         F a,b;
//         a.setByCSPRNG();b.setByCSPRNG();
//         for(int i = 0; i < segments/2; i++){
//             for(int j = 0; j < segment_size; j++){
//                 poly[i*segment_size +j] = poly[(2*i)*segment_size + j]*a + poly[(2*i+1)*segment_size + j]*b;
//             }
//             poly_compressed[i] = poly_compressed[2*i]*a+poly_compressed[2*i+1]*b;
//             partial_commitments[i] = partial_commitments[2*i]*a + partial_commitments[2*i+1]*b;
//         }
//     	ps += 2*32;
//         segments = segments/2;
//     }
//     vector<F> final_poly(segment_size);
//     for(int i = 0; i < segment_size; i++){
//         final_poly[i] = poly[i];
//     }
//
//     vector<G1> proof_phase2_1;
//     open(final_poly,r1,pp.base,proof_phase2_1);
//     // verify_proof(r1,pp.pp2,pp.H,proof_phase2_1);
//     // verify_proof(r1,pp.pp2,pp.H,proof_phase2_1);
//     hyperproofs_openall(final_poly,pp.hyper_proofs_base);
// 	ps += proof_phase2_1.size()*32;
// 	stop = std::chrono::high_resolution_clock::now();
// 	duration = duration_cast<milliseconds>(stop - start);
// 	cout << "FlexProofs in Hydraproofs" << " in " << duration.count() << "ms" << endl;
// 	std::cout << "Proof size in Hydraproofs: " << ps/1024.0 << " KB" << std::endl;
//
//
// 	hDigest root_hash_r1 = MTr1.back()[0];
// 	int target_index_r1 = 0;
// 	// get the column number
// 	int collumns_r1 = MTr1[0].size();
// 	vector<size_t> c_r1 = {
// 		(size_t)(target_index_r1 % collumns_r1),
// 		(size_t)(target_index_r1 / collumns_r1 * 4)
// 	};
// 	if (c_r1[0] >= MTr1[0].size()) {
// 		throw out_of_range("out of range");
// 	}
// 	vector<hDigest> proof_r1 = merkle_tree::merkle_tree_prover::open_tree(MTr1, c_r1, collumns_r1);
// 	// verify
// 	int leaf_hash_idx_r1 = target_index_r1 / 2;
// 	hDigest target_leaf_hash_r1 = MTr1[0][leaf_hash_idx_r1];
//
// 	start = high_resolution_clock::now();
// 	//verify all Merkle proofs
// 	int num=(logN-1)/2+2;
// 	for (int i=0; i<num; i++) {
// 		bool result = merkle_tree::merkle_tree_verifier::verify_claim(
// 			root_hash_r1,
// 			proof_r1.data(),
// 			target_leaf_hash_r1,
// 			target_index_r1,
// 			MTr1_leafs.size()
// 		);
// 	}
// 	verify_proof(r2,pp.pp2_compressed,pp.H,proof_compressed);
// 	verify_proof(r2,pp.pp2_compressed,pp.H,proof_compressed);
// 	verify_proof(r1,pp.pp2,pp.H,proof_phase2_1);
// 	verify_proof(r1,pp.pp2,pp.H,proof_phase2_1);
//
// 	F a,b;
// 	a.setByCSPRNG();b.setByCSPRNG();
// 	Fr verifierData = data[0];
// 	G1 verifierdCom = partial_commitments[0] * a;
// 	for(int k = 0; k < logN/2; k++){
// 		verifierData = verifierData*a+verifierData*b;
// 		verifierdCom = verifierdCom*a+verifierdCom*b;
// 	}
// 	verify_proof(r1,pp.pp2,pp.H,proof_phase2_1);
// 	stop = high_resolution_clock::now();
// 	duration = duration_cast<milliseconds>(stop - start);
// 	cout << "Verify in Hydraproofs in: " << duration.count() << "ms" << endl;
//
// }

void test_ours(int logN){
	cout << "========================Our VC=====================" << endl;
	//generate random numbers
	size_t N = 1ULL<<logN;
	int logmu = logN/2;
	size_t mu = 1ULL << logmu;
	F r;
	r.setByCSPRNG();

	vector<F> data(N);
	data[0].setByCSPRNG();
	for(int i = 1; i < N; i++){
		data[i] = data[i-1]*r;
	}
	// generate public parameters
	Fr rand;
	rand.setByCSPRNG();
	FCpp fc_pp = fc_gen(mu, rand);
	size_t kzp_size = 2 * mu - 1;
	KZGpp kzg_pp = kzg_setup(kzp_size, rand);
	Comm pp_pc;
	pst_generatePP(pp_pc,logmu);
	cout << "N = " << N << ", parameters generated" << endl;

	// generate commitment
	vector<vector<F>> partitions(mu);
	for(int i = 0; i < mu; i++){
		partitions[i].resize(mu);
		for(int j = 0; j < mu; j++){
			partitions[i][j] = data[i*mu+j];
		}
	}
	auto start = high_resolution_clock::now();
	vector<G1> c_vec;
	c_vec.reserve(mu);
	for (size_t i = 0; i < mu; ++i) {
		pst_commit(partitions[i],pp_pc);
		c_vec.push_back(pp_pc.C);
	}
	Fp12 com = fc_com(fc_pp, c_vec);
	auto end = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(end - start);
	cout << "com generated: " << duration.count() << "ms" << endl;

	int ps = 0;
	// cout << "ps:" << ps << endl;
	//OpanAll-Step 1: generate proof for (\log \mu)-sized batch
	size_t batchSize = logN * logN;//static_cast<int>(sqrt(logN));
	cout << "batchSize:" << "logN * logN" << endl;
	std::vector<std::vector<size_t>> index_batches;
	std::vector<std::vector<G1>> value_batches;

	for (size_t i = 0; i < mu; i += batchSize) {
		size_t end = std::min(i + batchSize, mu);
		// cout << "end:" << end << endl;

		std::vector<size_t> idx_batch;
		std::vector<G1> val_batch;

		for (size_t j = i; j < end; ++j) {
			idx_batch.push_back(j);
			val_batch.push_back(c_vec[j]);
		}
		index_batches.push_back(idx_batch);
		value_batches.push_back(val_batch);
	}
	std::vector<PI> batch_proofs;
	size_t batchNum = index_batches.size();
	start = high_resolution_clock::now();
	for (size_t i = 0; i < batchNum; ++i) {
		PI batch_pi = svc_batch_open(fc_pp, kzg_pp, com, c_vec, index_batches[i], value_batches[i]);
		batch_proofs.push_back(batch_pi);
	}
	ps += batchSize*32;
	int batch_pi_size = 2*(logmu*(32+384))+32+2*64;//L_vec\|R_vec\|finalA\|finalv\|finalv_proof
	ps += batch_pi_size;

	//OpanAll-Step 2
	//random numbers in the non-interactive version
	vector<vector<uint8_t>> c_vec_ser;
	for (int i=0; i<mu; i++) {
		c_vec_ser.push_back(G1_serialize(c_vec[i]));
	}
	vector<Fr> MTrand_leaves;
	for(int i = 0; i < data.size(); i++) {
		Fr res;
		int comIdx = i/mu;
		std::vector<uint8_t> leafData = Fr_serialize(data[i]);
		leafData.insert(leafData.end(), c_vec_ser[comIdx].begin(), c_vec_ser[comIdx].end());
		res.setHashOf(leafData.data(), leafData.size());
		MTrand_leaves.push_back(res);
	}

	vector<vector<hDigest>> MTrand;
	merkle_tree::merkle_tree_prover::MT_commit(MTrand_leaves,N,MTrand);

	ps += logN*32;
	vector<vector<F>> random_polys(mu);
	vector<G1> random_c_vec(mu);
	for(int i = 0; i < mu; i++){
		F a;
		a.setByCSPRNG();
		random_polys[i].resize(mu);
		for(int j = 0; j < mu; j++){
			random_polys[i][j] = a * partitions[i][j];
		}
		random_c_vec[i] = c_vec[i] * a;
	}
	auto segments = mu;
	for(int k = 0; k < logmu; k++){
		for(int i = 0; i < segments/2; i++){
			for(int j = 0; j < mu; j++){
				random_polys[i][j] = random_polys[2*i][j] + random_polys[2*i+1][j];
			}
			random_c_vec[i] = random_c_vec[2*i] + random_c_vec[2*i+1];
		}
		segments = segments/2;
		ps += 2*32;
	}
	vector<F> final_poly(mu);
	for(int i = 0; i < mu; i++){
		final_poly[i] = random_polys[0][i];
	}
	pst_openall(final_poly,pp_pc.base);
	ps += ((int)log2(final_poly.size()))*32;
	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(end - start);
	cout << "FlexProofs in Our VC in: " << duration.count() << "ms" << endl;
	std::cout << "Proof size in our VC: " << ps/1024.0 << " KB" << std::endl;

	hDigest root_hash = MTrand.back()[0];
	int target_index = 0;
	// get the column number
	int collumns = MTrand[0].size();
	vector<size_t> c = {
		(size_t)(target_index % collumns),
		(size_t)(target_index / collumns * 4)
	};
	if (c[0] >= MTrand[0].size()) {
		throw out_of_range("out of range");
	}
	vector<hDigest> proof = merkle_tree::merkle_tree_prover::open_tree(MTrand, c, collumns);

	vector<F> r1(logmu);
	for(int i = 0; i < r1.size(); i++){
		r1[i].setByCSPRNG();
	}
	vector<G1> proof_phase2_1;
	pst_open(final_poly,r1,pp_pc.base,proof_phase2_1);
	// verify
	int leaf_hash_idx = target_index / 2;
	hDigest target_leaf_hash = MTrand[0][leaf_hash_idx];

	start = high_resolution_clock::now();
	//Verification in Step 1
	svc_batch_verify(fc_pp, kzg_pp, com, index_batches[0], value_batches[0], batch_proofs[0]);
	//Verification in Step 2
	merkle_tree::merkle_tree_verifier::verify_claim(
		root_hash,
		proof.data(),
		target_leaf_hash,
		target_index,
		N
	);
	F a;
	a.setByCSPRNG();
	Fr verifierData = data[0] * a;
	G1 verifierdCom = value_batches[0][0] * a;
	for(int k = 0; k < logmu; k++){
		verifierData += verifierData;
		verifierdCom += verifierdCom;
	}
	pst_verify(r1,pp_pc.pp2,pp_pc.H,proof_phase2_1);
	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(end - start);
	cout << "Verify in Our VC in: " << duration.count() << "ms" << endl;
}

void test_svc(uint32_t exp, size_t batchsize) {
    cout << "========================SVC=====================" << endl;
    size_t n = 1 << exp;  // 2^exp
    size_t t = batchsize; // batch size > 1

	Fr rand;
	rand.setByCSPRNG();
    FCpp pp = fc_gen(n, rand);
    size_t kzp_size = 2 * n - 1;
    KZGpp kzg_pp = kzg_setup(kzp_size, rand);
    cout << "n = " << n << ", t = " << t << ", parameters generated" << endl;

    vector<G1> a_vec;
    a_vec.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        Fr randomNum;
    	randomNum.setByCSPRNG();

        G1 g1;
    	mcl::bn::hashAndMapToG1(g1, "g_1");// get a generate in G2
        g1 *= randomNum;  //  scalar multiplication
        a_vec.push_back(g1);
    }
    cout << "values generated" << endl;

    // generate a commitment for a_vec
    auto start = high_resolution_clock::now();
    Fp12 com = fc_com(pp, a_vec);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);
    cout << "com generated: " << duration.count() << "ms" << endl;

    // open in batch
    vector<size_t> vecI;
    vector<G1> subvec_value;

    for (size_t j = 0; j < batchsize; ++j) {
        vecI.push_back(j);
        subvec_value.push_back(a_vec[j]);
    }

    start = high_resolution_clock::now();
    PI batch_pi = svc_batch_open(pp, kzg_pp, com, a_vec, vecI, subvec_value);
    end = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(end - start);
    cout << "Open subvector of size " << batchsize << " in " << duration.count() << "ms" << endl;

    // verify in batch
    start = high_resolution_clock::now();
    bool res = svc_batch_verify(pp, kzg_pp, com, vecI, subvec_value, batch_pi);
    end = high_resolution_clock::now();
    duration = duration_cast<milliseconds>(end - start);
    cout << "Batch verify result is " << boolalpha << res
         << ", the required time is " << duration.count() << "ms" << endl;
}


void test_fc(uint32_t exp, size_t batchsize) {
	cout << "========================FC=====================" << endl;
	size_t n = 1 << exp;  // 2^exp
	size_t s = batchsize; // batch size > 1

	Fr rand;
	rand.setByCSPRNG();
	FCpp pp = fc_gen(n, rand);
	size_t kzp_size = 2 * n - 1;
	KZGpp kzg_pp = kzg_setup(kzp_size, rand);
	cout << "n = " << n << ", s = " << s << ", parameters generated" << endl;

	vector<G1> a_vec;
	a_vec.reserve(n);
	for (size_t i = 0; i < n; ++i) {
		Fr randomNum;
		randomNum.setByCSPRNG();

		G1 g1;
		mcl::bn::hashAndMapToG1(g1, "g_1");// get a generate in G2
		g1 *= randomNum;  //  scalar multiplication
		a_vec.push_back(g1);
	}

	std::vector<std::vector<Fr>> b_matrix;
	b_matrix.reserve(s);
	for (size_t i = 0; i < s; ++i) {
		std::vector<Fr> row;
		row.reserve(n);

		for (size_t j = 0; j < n; ++j) {
			Fr randomNum;
			randomNum.setByCSPRNG();
			row.push_back(randomNum);
		}

		b_matrix.push_back(std::move(row));
	}
	cout << "values generated" << endl;

	// generate a commitment for a_vec
	auto start = high_resolution_clock::now();
	Fp12 com = fc_com(pp, a_vec);
	auto end = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(end - start);
	cout << "com generated: " << duration.count() << "ms" << endl;

	//Direct computation
	vector<G1> y_vec;
	y_vec.reserve(s);
	start = high_resolution_clock::now();
	for (size_t i = 0; i < s; ++i) {
		G1 y;
		G1::mulVec(y, a_vec.data(), b_matrix[i].data(), n);
		y_vec.push_back(y);
	}
	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(end - start);
	// cout << "The time cost by " << s << " direct computation is " << duration.count() << "ms" << endl;

	// FC.BOpen
	start = high_resolution_clock::now();
	PI batch_pi = fc_batch_open(pp, kzg_pp, com, a_vec, b_matrix, y_vec);
	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(end - start);
	cout << "Generate batch proof for " << batchsize << " in " << duration.count() << "ms" << endl;

	// FC.BVerify
	start = high_resolution_clock::now();
	bool res = fc_batch_verify(pp, kzg_pp, com, b_matrix, y_vec, batch_pi);
	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(end - start);
	cout << "Batch verify result is " << boolalpha << res
		 << ", the required time is " << duration.count() << "ms" << endl;
}

int main() {
	initPairing(mcl::BN_SNARK1);
	init_hash();

	// const size_t fr_size = Fr::getByteSize();
	// cout<< "Fr size = " << fr_size << endl;
	// G1 P;
	// G2 Q;
	// GT R;
	// mapToG1(P, 1);
	// mapToG2(Q, 1);
	// pairing(R, P, Q);
	// size_t G1_size = P.getSerializedByteSize();
	// cout<< "G1 size = " << G1_size << endl;
	// size_t G2_size = Q.getSerializedByteSize();
	// cout<< "G2 size = " << G2_size << endl;

	// // test_Hyperproofs(18);
	for(int i = 16; i <= 24; i=i+2){
		cout << i << endl;
		test_ours(i);
		// test_Hydraproof(i);
	}

	// test_svc(10, 1);
	// test_fc(12, 16);

	return 0;
}

int test_serialize() {
	initPairing(BN_SNARK1);

	Fr x;
	x.setStr("11440409931983508671476017804487278721750219707868305868837093894100494299711", 10);

	std::cout << "x = " << x.getStr() << std::endl;
	std::cout << "x.isValid() = " << x.isValid() << std::endl;

	const size_t size = Fr::getByteSize();
	std::vector<uint8_t> buf(size);
	size_t written = x.serialize(buf.data(), size, mcl::IoSerialize);

	if (written == 0) {
		std::cerr << "Serialization failed!" << std::endl;
		return 1;
	}

	std::cout << "Serialization success: written = " << written << std::endl;
}