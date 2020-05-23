#include <stdio.h>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>
#include "viterbi.h"

using namespace std;

int print_matrix(string name, float* matrix, int H, int W, int C) {
    printf("Matrix %s: %d %d %d\n", name.c_str(), H, W, C);
    for (int h=0; h<H; h++) {
        printf("\t| ");
        for (int w=0; w<W; w++) {
            printf("\t| ");
            for (int c=0; c<C; c++) {
                int idx = h*W*C + w*C + c;
                // printf("(%d - %d %d %d)%f ", idx, h, w, c, matrix[idx]);
                printf("%f ", matrix[idx]);
            }
            printf("| ");
        }
        printf("|\n");
    }
}

vector<string> split_string(string str, string delimiter) {
    size_t pos = 0;
    vector<string> splits;
    while ((pos = str.find(delimiter)) != string::npos) {
        splits.push_back(str.substr(0, pos));
        str.erase(0, pos + delimiter.length());
    }
    splits.push_back(str);
    return splits;
}

int read_matrix(string filename, float* outMatrix) {
    ifstream fin (filename);
    string line;
    getline(fin, line);

    vector<string> splits = split_string(line, string(" "));
    for (int i=0; i<splits.size(); i++) {
        outMatrix[i] = atof(splits.at(i).c_str());
    }
    return splits.size();
}



int main(int argc, char* agv[]) {
    float* pi = new float[ALPHABET_LEN];
    float* a = new float[ALPHABET_LEN*ALPHABET_LEN];
    float* b = new float[ALPHABET_LEN*TIME_LEN];
    alloc();

    int alphabet_size = read_matrix(string("../data/pi.txt"), pi);
    read_matrix(string("../data/a.txt"), a);
    int b_size = read_matrix(string("../data/b.txt"), b);
    int time_size = b_size/alphabet_size;

    print_matrix(string("pi"), pi, alphabet_size, 1, 1);
    print_matrix(string("a"), a, alphabet_size, alphabet_size, 1);
    print_matrix(string("b"), b, alphabet_size, time_size, 1);

    printf("Runing viterbi with:\n\talphabet_size: %d\n\ttime_size: %d\n--------",
            alphabet_size, time_size);

    auto start = chrono::steady_clock::now();
    for (int i=0; i<1; i++) {
        vector<vector<int>> paths = ocr_viterbi_topk(pi, a, b, alphabet_size, time_size, 10);
    }
    auto end = chrono::steady_clock::now();
    printf("Elapsed time in milliseconds : %ld ms\n", 
		chrono::duration_cast<chrono::milliseconds>(end - start).count());

    delete []pi;
    delete []a;
    delete []b;
    free();

    return 0;
}
