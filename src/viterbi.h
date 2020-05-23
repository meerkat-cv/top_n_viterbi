#pragma once

#include <stdio.h>
#include <queue>
#include <vector>
#include <unordered_map>


// This defines the maximum viterbi search
#define ALPHABET_LEN 200
#define TIME_LEN 5000
#define TOP_N 50

struct ProbState {
    ProbState(float prob, int state1, int state2) {
        this->prob = prob;
        this->state1 = state1;
        this->state2 = state2;
    };
    float prob;
    int state1, state2;
    bool operator<(const ProbState& rhs) const
    {
        return prob < rhs.prob;
    }
};

extern "C" {
    void print_stuff();
    void print_numpy_float32_array(float *array, int size);
}

// Allocating all those variables only once to avoid several
// allocations/deallocations on what should be a very fast
// function. This makes it NOT THREAD SAFE!!!!
extern float *g_delta, *g_phi, *g_delta_top, *g_phi_top, *g_rank;

// Main functions
std::vector<int> ocr_viterbi(const float *pi, const float *a, const float *b, int nStates, int T);
std::vector<std::vector<int>> ocr_viterbi_topk(const float *pi, const float *a, const float *b, int nStates, int T, int topK);
void free();

// Helper functions
void alloc();
int get_coord(int H, int h, int W, int w, int C, int c);
