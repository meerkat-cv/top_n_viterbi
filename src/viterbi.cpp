#include "viterbi.h"
#include <vector>
#include <algorithm>

using namespace std;


float *g_delta=NULL, *g_phi=NULL, *g_delta_top=NULL, *g_phi_top=NULL, *g_rank=NULL;
vector<ProbState> g_state_vec = vector<ProbState>(ALPHABET_LEN*ALPHABET_LEN,ProbState());
vector<ProbStateRank> g_state_rank_vec = vector<ProbStateRank>(ALPHABET_LEN*TOP_N,ProbStateRank());

int get_coord(int H, int h, int W, int w, int C, int c)
{
    return h*W*C + w*C + c;
}

void ocr_viterbi(const float *pi, const float *a, const float *b, int nStates, int T, int *outPath, float *outProbs)
{
    if (g_delta == NULL) {
        alloc();
    }

    for (int i=0; i<nStates; i++) {
        for (int t=1; t<T; t++) {
            g_delta[i*T+t] = 0;
            g_phi[i*T+t] = 0;
        }
        g_delta[i*T] = pi[i] * b[i*T+0];
        g_phi[i*T] = 0;
    }

    for (int t=1; t<T; t++) {
        for (int s=0; s<nStates; s++) {
            float max = 0;
            int max_idx = 0;
            for (int k=0; k<nStates; k++) {
                float val = g_delta[k*T + t-1] * a[k*nStates + s];
                if (val > max) {
                    max = val;
                    max_idx = k;
                }
            }
            g_delta[s*T + t] = max * b[s*T + t];
            g_phi[s*T + t] = max_idx;
        }
    }

    int max_delta_idx = -1;
    float max = 0;
    for (int i=0; i<nStates; i++) {
        if (g_delta[i*T + T-1] > max) {
            max = g_delta[i*T + T-1];
            max_delta_idx = i;
        }
    }
    outProbs[0] = max;
    outPath[T-1] = max_delta_idx;
    for (int t=T-2; t>-1; t--) {
        outPath[t] = g_phi[outPath[t+1]*T + t + 1];
    }
}


void ocr_viterbi_topk(const float *pi, const float *a, const float *b, int nStates, int T, int topK, int *outPaths, float *outProbs)
{
    if (g_delta == NULL) {
        alloc();
    }
    if (topK == 1) {
        ocr_viterbi(pi, a, b, nStates, T, outPaths, outProbs);
    }

    for (int t=0; t<T; t++) {
        for (int s=0; s<nStates; s++) {
            for (int k=0; k<topK; k++) {
                int idx = get_coord(T, t, nStates, s, topK, k);
                g_delta_top[idx] = 0;
                g_phi_top[idx] = 0;
                g_rank[idx] = 0;
            }
        }
    }

    for (int i=0; i<nStates; i++) {
        int idx = get_coord(T, 0, nStates, i, topK, 0);
        g_delta_top[idx] = pi[i] * b[i*T];
        g_phi_top[idx] = i;

        for (int k=1; k<topK; k++) {
            int idx = get_coord(T, 0, nStates, i, topK, k);
            g_phi_top[idx] = i;
        }
    }


    for (int t=1; t<T; t++) {
        for (int s1=0; s1<nStates; s1++) {
            int vec_idx=0;
            for (int s2=0; s2<nStates; s2++) {
                for (int k=0; k<topK; k++) {
                    int idx = get_coord(T, t-1, nStates, s2, topK, k);
                    float prob = g_delta_top[idx] * a[s2*nStates+s1] * b[s1*T + t];
                    g_state_vec[vec_idx].prob = prob;
                    g_state_vec[vec_idx++].state = s2;
                }
            }

            partial_sort(g_state_vec.begin(), g_state_vec.begin()+topK, g_state_vec.begin()+vec_idx);

            unordered_map<int,int> rankDict;

            for (int k=0; k<topK; k++) {
                int idx = get_coord(T, t, nStates, s1, topK, k);
                ProbState ps = g_state_vec[k];
                g_delta_top[idx] = ps.prob;
                g_phi_top[idx] = ps.state;
                int state = ps.state;

                if (rankDict.find(state) != rankDict.end()) {
                    rankDict[state] += 1;
                } else {
                    rankDict[state] = 0;
                }

                g_rank[idx] = rankDict[state];
            }
        }
    }

    // print_matrix(string("Rank"), g_rank, T, nStates, topK);
    // print_matrix(string("Delta"), g_delta_top, T, nStates, topK);
    // print_matrix(string("Phi"), g_phi_top, T, nStates, topK);

    int vec_idx = 0;
    for (int s1=0; s1<nStates; s1++) {
        for (int k=0; k<topK; k++) {
            int idx = get_coord(T, T-1, nStates, s1, topK, k);
            float prob = g_delta_top[idx];
            g_state_rank_vec[vec_idx].prob = prob;
            g_state_rank_vec[vec_idx].state1 = s1;
            g_state_rank_vec[vec_idx++].state2 = k;
        }
    }

    partial_sort(g_state_rank_vec.begin(), g_state_rank_vec.begin()+topK, g_state_rank_vec.begin()+vec_idx);

    // Now backtrace for k and each time stamp
    for (int k=0; k<topK; k++) {
        // The maximum probability and the state it came from
        float max_prob = g_state_rank_vec[k].prob;
        int state = g_state_rank_vec[k].state1;
        int rankK = g_state_rank_vec[k].state2;

        // Assign to output arrays
        outPaths[k*T+T-1] = state;
        outProbs[k] = max_prob;

        // Then from t down to 0 store the correct sequence for t+1
        for (int t=T-2; t>-1; t--) {
            int nextState = outPaths[k*T+t+1];

            int idx = get_coord(T, t+1, nStates, nextState, topK, rankK);
            float p = g_phi_top[idx];

            outPaths[k*T+t] = p;

            rankK = g_rank[idx];
        }
    }
}

void alloc() {
    g_delta = new float[ALPHABET_LEN*TIME_LEN];
    g_phi = new float[ALPHABET_LEN*TIME_LEN];
    g_delta_top = new float[ALPHABET_LEN*TIME_LEN*TOP_N];
    g_phi_top = new float[ALPHABET_LEN*TIME_LEN*TOP_N];
    g_rank = new float[ALPHABET_LEN*TIME_LEN*TOP_N];
}

void free_variables() {
    if (g_delta != NULL) {
        delete []g_delta;
        delete []g_phi;
        delete []g_delta_top;
        delete []g_phi_top;
        delete []g_rank;
        g_delta = g_phi = g_delta_top = NULL;
        g_phi_top = g_rank = NULL;
    }
}
