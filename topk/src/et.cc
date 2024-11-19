/**
    USI: useful string indexing
    Copyright (C) 2024 Solon P. Pissis.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <omp.h>
#include <sys/time.h>

#include <chrono>
#include <cmath>
#include <ctime>
#include <iostream>
#include <sdsl/suffix_trees.hpp>
#include <string>
#include <tuple>
#include <vector>

#include "defs.h"
#include "krfp.h"
#include "unordered_dense.h"
#include "useful_methods.h"
#include "utils.h"

using namespace sdsl;
using namespace std;

int main(int argc, char *argv[]) {
    if (argc != 4) {
        cout << "Usage: " << argv[0] << " text-file " << " weights-file " << " K " << endl;
        return -1;
    }

    string second_arg(argv[2]);
    INT K = atoll(argv[3]);

    ifstream seq(argv[1], ios::in | ios::binary);
    unsigned char *input_seq_char = NULL;
    unsigned char *sequence;
    INT n = 0;
    char c;
    while (seq.get(c)) {
        if (n == 0 || n % ALLOC_SIZE) input_seq_char = (unsigned char *)realloc(input_seq_char, (n + ALLOC_SIZE) * sizeof(unsigned char));
        input_seq_char[n] = c;
        n++;
    }
    seq.close();
    if (n == 0 || n % ALLOC_SIZE) input_seq_char = (unsigned char *)realloc(input_seq_char, (n + ALLOC_SIZE) * sizeof(unsigned char));
    if (c == '\n')
        input_seq_char[--n] = 255;
    else
        input_seq_char[n] = 255;
    sequence = input_seq_char;
    sequence[++n] = '\0';

    cout << "Text is of length n = " << n - 1 << "." << endl;

    std::chrono::steady_clock::time_point construct_begin = std::chrono::steady_clock::now();

    std::chrono::steady_clock::time_point tup_construct_begin = std::chrono::steady_clock::now();

    INT *SA = (INT *)malloc((n) * sizeof(INT));
    if ((SA == NULL)) {
        fprintf(stderr, " Error: Cannot allocate memory for SA.\n");
        return (0);
    }

    INT *LCP = (INT *)calloc(n, sizeof(INT));
    if ((LCP == NULL)) {
        fprintf(stderr, " Error: Cannot allocate memory for LCP.\n");
        return (0);
    }

    B *b = (B *)realloc(b, (n) * sizeof(B));
    if ((b == NULL)) {
        fprintf(stderr, " Error: Cannot allocate memory for b.\n");
        return (0);
    }

    if (divsufsort64(sequence, SA, n) != 0) {
        fprintf(stderr, " Error: SA computation failed.\n");
        exit(EXIT_FAILURE);
    }

    INT *invSA = (INT *)calloc(n, sizeof(INT));
    if ((invSA == NULL)) {
        fprintf(stderr, " Error: Cannot allocate memory for invSA.\n");
        return (0);
    }
    for (INT i = 0; i < n; i++) invSA[SA[i]] = i;

    if (LCParray(sequence, n, SA, invSA, LCP) != 1) {
        fprintf(stderr, " Error: LCP computation failed.\n");
        exit(EXIT_FAILURE);
    }

    free(invSA);
    cout << "SA and LCP array are constructed." << endl;

    construct_tuples(sequence, n, SA, LCP, b);

    std::chrono::steady_clock::time_point tup_construct_end = std::chrono::steady_clock::now();
    cout << "Tuple construction time: " << std::chrono::duration_cast<std::chrono::milliseconds>(tup_construct_end - tup_construct_begin).count() << "[ms]." << std::endl;

    INT bsize = 0;
    INT tau = n;
    INT L = 0;
    INT w = 0;
    find_longest(b, n, K, tau, L, bsize, w);
    b = (B *)realloc(b, (bsize) * sizeof(B));

    vector<double> PS;
    if (second_arg.compare("random") == 0) {
        cout << "Random weights generation.\n";
        random_weight_generation(PS, n - 1);
    } else {
        cout << "The weights are read from file.\n";
        weight_from_file(argv[2], PS, n - 1);
    }

    std::chrono::steady_clock::time_point util_comp_begin = std::chrono::steady_clock::now();

    auto H = ankerl::unordered_dense::map<INT, double>();
    substrings_utility(SA, sequence, PS, b, bsize, H, string(argv[1]), K);
    cout << "The computed tau is " << tau << "." << endl;
    cout << "The longest frequent pattern is of length L = " << L << "." << endl;
    free(sequence);
    free(LCP);
    free(SA);
    return 0;
}
