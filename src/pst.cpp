#include "pst.h"

#include <cmath>
#include <vector>
using namespace std;

// int threads = 1;
void precomputeBeta(F* destMatrix, F* tempMatrix, F randomVal, vector<F>& betaVec, int index, int length) {
    int elemCount = length;
    int startPos = index * elemCount;

    if (betaVec.size() / 2 == elemCount) {
        for (int j = startPos; j < startPos + elemCount; j++) {
            int vecIndex = j << 1;
            F tempVal = randomVal * tempMatrix[j];
            betaVec[vecIndex] = tempMatrix[j] - tempVal;
            betaVec[vecIndex + 1] = tempVal;
        }
    } else {
        for (int j = startPos; j < startPos + elemCount; j++) {
            int matIndex = j << 1;
            F tempVal = randomVal * tempMatrix[j];
            destMatrix[matIndex] = tempMatrix[j] - tempVal;
            destMatrix[matIndex + 1] = tempVal;
        }
    }
}

void computeBinary(const vector<F>& inputVec, vector<F>& resultVec) {
    // Calculate the size of result vector as 2^N where N is the size of input vector
    const size_t resultSize = 1 << inputVec.size();
    resultVec.resize(resultSize);

    // Iterate over all possible binary combinations (from 0 to 2^N - 1)
    for (size_t i = 0; i < resultSize; ++i) {
        size_t bitMask = i;
        resultVec[i] = F(1);

        // Check each bit in the current bit mask to compute the product
        for (size_t j = 0; j < inputVec.size(); ++j) {
            if ((bitMask & 1) == 1) {
                resultVec[i] *= inputVec[j];
            }

            bitMask >>= 1;
        }
    }
}


void generatePP1(vector<G1>& commitments, const vector<F>& betaVals, G1 baseG1, int index) {
    const int segmentSize = betaVals.size() / 16;
    const int startPos = index * segmentSize;

    for (int i = startPos; i < startPos + segmentSize; ++i) {
        commitments[i] = baseG1 * betaVals[i];
    }
}

// Precomputes beta values into outputB using inputR, with iterative expansion
void precomputeBeta(const vector<F>& inputR, vector<F>& outputB) {
    outputB.resize(1 << inputR.size());  // Resize to 2^N (N = inputR size)
    F* tempArray = (F*)malloc(1 * sizeof(F));
    outputB[0] = F(1);
    tempArray[0] = F(1);  // Initialize base case

    for (int i = 0; i < inputR.size(); ++i) {
        const size_t nextTempSize = 1 << (i + 1);  // Next temp size: 2^(i+1)
        F* nextTempArray = (F*)malloc(nextTempSize * sizeof(F));

        // Threshold: handle small i (i <=10) directly, larger via helper
        if ((1 << 10) >= (1 << i)) {
            // Check if outputB is target (half its size == 2^i)
            if (outputB.size() / 2 == (1 << i)) {
                // Update outputB: split each tempArray[j] into two elements
                for (int j = 0; j < (1 << i); ++j) {
                    const int idx = j << 1;  // j*2
                    const F rVal = inputR[inputR.size() - 1 - i];  // Reverse index from inputR
                    const F tempVal = rVal * tempArray[j];
                    outputB[idx] = tempArray[j] - tempVal;
                    outputB[idx + 1] = tempVal;
                }
            } else {
                // Update nextTempArray instead of outputB
                for (int j = 0; j < (1 << i); ++j) {
                    const int idx = j << 1;
                    const F rVal = inputR[inputR.size() - 1 - i];
                    const F tempVal = rVal * tempArray[j];
                    nextTempArray[idx] = tempArray[j] - tempVal;
                    nextTempArray[idx + 1] = tempVal;
                }
            }
        } else {
            // Delegate to helper for large i
            precomputeBeta(nextTempArray, tempArray, inputR[inputR.size() - 1 - i], outputB, 0, 1 << i);
        }

        tempArray = nextTempArray;  // Move to next temp array
    }
}

// Generates public parameters and commitments for PST (with specified bit size)
void pst_generatePP(Comm& commitment, int bitSize) {
    vector<F> randomVals;
    char randBytes[256];

    // Initialize random byte array (simple pattern for demonstration)
    for (int i = 0; i < 256; ++i) {
        randBytes[i] = i + '4';  // '4' is ASCII base for initial values
    }

    // Hash bytes to G1/G2 base points for commitment
    mcl::bn::hashAndMapToG1(commitment.G, randBytes, 256);
    mcl::bn::hashAndMapToG2(commitment.H, randBytes, 256);
    commitment.bit_size = bitSize;  // Store bit size

    // Generate random secret values (CSPRNG = cryptographically secure PRNG)
    for (int i = 0; i < commitment.bit_size; ++i) {
        F temp;
        temp.setByCSPRNG();
        randomVals.push_back(temp);
    }
    commitment.secret = randomVals;  // Store secrets

    // Precompute beta values using helper functions
    vector<F> betas, betas2;
    precomputeBeta(randomVals, betas2);
    computeBinary(randomVals, betas);

    // Resize and generate pp1 commitments (split into 16 segments)
    commitment.pp1.resize(betas.size());
    for (int i = 0; i < 16; ++i) {
        generatePP1(commitment.pp1, betas2, commitment.G, i);
    }

    // Build base layers by combining consecutive elements
    vector<G1> tempG1 = commitment.pp1;
    commitment.base.resize(randomVals.size());
    for (int i = 0; i < randomVals.size(); ++i) {
        for (int j = 0; j < tempG1.size() / 2; ++j) {
            commitment.base[i].push_back(tempG1[2 * j] + tempG1[2 * j + 1]);
        }
        tempG1 = commitment.base[i];  // Update temp for next layer
    }

    // Generate pp2 commitments (H multiplied by each secret value)
    for (int i = 0; i < randomVals.size(); ++i) {
        commitment.pp2.push_back(commitment.H * randomVals[i]);
    }
}

void pst_commit(vector<Fr> poly, struct Comm &commitment){
    G1::mulVec(commitment.C, commitment.pp1.data(), poly.data(), poly.size());
}

vector<F> arrayDifference(F x1, F x2, const vector<F>& arr) {
    vector<F> diff(arr.size() / 2);  // Result size is half the input array

    if (x1 == F(1) && x2 == F(1)) {
        // Special case: element = (second element of pair) - (first element of pair)
        for (int i = 0; i < arr.size() / 2; ++i) {
            diff[i] = arr[2 * i + 1] - arr[2 * i];
        }
    } else {
        // General case: element = x1*(pair difference) + first element of pair
        for (int i = 0; i < arr.size() / 2; ++i) {
            diff[i] = x1 * (arr[2 * i + 1] - arr[2 * i]) + arr[2 * i];
        }
    }

    return diff;
}

vector<F> computeCoefficients(const vector<F>& arr) {
    // Base case: if array has one element, return it as the sole coefficient
    if (arr.size() == 1) {
        return { arr[0] };
    }

    // Split array into left and right halves
    const int halfSize = arr.size() / 2;
    vector<F> tempArr(halfSize);

    // Process left half
    for (int i = 0; i < halfSize; ++i) {
        tempArr[i] = arr[i];
    }
    vector<F> leftCoeffs = computeCoefficients(tempArr);

    // Process right half
    for (int i = 0; i < halfSize; ++i) {
        tempArr[i] = arr[i + halfSize];
    }
    vector<F> rightCoeffs = computeCoefficients(tempArr);

    // Combine results: left coefficients followed by (right - left) coefficients
    vector<F> result(2 * leftCoeffs.size());
    for (int i = 0; i < rightCoeffs.size(); ++i) {
        result[i] = leftCoeffs[i];
    }
    for (int i = 0; i < leftCoeffs.size(); ++i) {
        result[i + leftCoeffs.size()] = rightCoeffs[i] - leftCoeffs[i];
    }

    return result;
}

vector<vector<F>> findQuotients(const vector<F>& arr, const vector<F>& x) {
    vector<vector<F>> quotients(x.size());
    vector<F> currentArr = arr;

    // Process first (x.size() - 1) elements to compute intermediate quotients
    for (int i = 0; i < x.size() - 1; ++i) {
        vector<F> diff = arrayDifference(F(1), F(1), currentArr);
        quotients[i] = computeCoefficients(diff);

        currentArr = arrayDifference(x[i], x[i] - F(1), currentArr);
    }

    // Process last element: compute coefficients of final array and extract specific value
    vector<F> finalCoeffs = computeCoefficients(currentArr);
    quotients[x.size() - 1].push_back(finalCoeffs[1]);

    return quotients;
}

void pst_open(vector<F> poly, vector<F> x, vector<vector<G1>> &pp, vector<G1> &Proof){
    vector<vector<F>> quotients = findQuotients(poly,x);
    for(int i = 0; i < x.size(); i++){
        G1 proofElem;
        G1::mulVec(proofElem, pp[i].data(),quotients[i].data(), quotients[i].size());
        Proof.push_back(proofElem);
    }
}

void pst_verify(vector<F> x, vector<G2> &pp2,G2 H, vector<G1> &Proof){
    vector<G2> diffs;
    GT Pairing_prod;

    for(int i = 0; i < x.size(); i++){
        diffs.push_back(pp2[i] + H*(F(0)-x[i]));
    }
    millerLoopVec(Pairing_prod,Proof.data(),diffs.data(),x.size()+1);
    finalExp(Pairing_prod,Pairing_prod);
}

void computeProofs(vector<F> &poly, vector<vector<G1>> &pp,vector<vector<G1>> &proof_tree, int level){
    if(poly.size() == 1){
        return;
    }
    const int halfSize = poly.size() / 2;
    vector<F> poly1(halfSize),poly2(halfSize);
    vector<F> quotient(halfSize);
    for(int i = 0; i < halfSize; i++){
        poly1[i] = poly[2*i];
        poly2[i] = poly[2*i + 1];
        quotient[i] = poly2[i] - poly1[i];
    }
    G1 proofElem;
    G1::mulVec(proofElem, pp[level].data(),quotient.data(), quotient.size());

    proof_tree[level].push_back(proofElem);

    computeProofs(poly1,pp,proof_tree,level+1);
    computeProofs(poly2,pp,proof_tree,level+1);
}

void pst_openall(vector<F> poly,vector<vector<G1>> &pp){
    const int numLevels = static_cast<int>(log2(poly.size()));
    vector<vector<G1>> proofTree(numLevels);
    computeProofs(poly,pp,proofTree,0);
}
