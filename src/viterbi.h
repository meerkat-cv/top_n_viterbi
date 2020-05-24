#pragma once

#include <stdio.h>
#include <queue>
#include <vector>
#include <unordered_map>


// This defines the maximum viterbi search
#define ALPHABET_LEN 100
#define TIME_LEN 1000
#define TOP_N 20

struct ProbState {
    ProbState() {};
    ProbState(float prob, int state) {
        this->prob = prob;
        this->state = state;
    };
    float prob;
    int state;
    // UGLY HACK TO DO A REVERSE ORDER SORT!!!
    bool operator<(const ProbState& rhs) const
    {
        return prob > rhs.prob;
    }
};

struct ProbStateRank {
    ProbStateRank() {};
    ProbStateRank(float prob, int state1, int state2) {
        this->prob = prob;
        this->state1 = state1;
        this->state2 = state2;
    };
    float prob;
    int state1, state2;
    // UGLY HACK TO DO A REVERSE ORDER SORT!!!
    bool operator<(const ProbStateRank& rhs) const
    {
        return prob > rhs.prob;
    }
};

extern "C" {
    // outPath memory allocation/deallocation is responsability of the caller.
    // It must have topK*nStates length.
    void ocr_viterbi(const float *pi, const float *a, const float *b, int nStates, int T, int *outPath, float *outProbs);
    void ocr_viterbi_topk(const float *pi, const float *a, const float *b, int nStates, int T, int topK, int *outPaths, float *outProbs);
    void free_variables();
}

// Allocating all those variables only once to avoid several
// allocations/deallocations on what should be a very fast
// function. This makes it NOT THREAD SAFE!!!!
extern float *g_delta, *g_phi, *g_delta_top, *g_phi_top, *g_rank;
extern std::vector<ProbState> g_state_vec;
extern std::vector<ProbStateRank> g_state_rank_vec;

// Main functions

// Helper functions
void alloc();
int get_coord(int H, int h, int W, int w, int C, int c);
