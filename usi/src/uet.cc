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

#include <malloc.h>
#include <omp.h>
#include <stdint.h>
#include <stdlib.h>
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

INT global_space = 0;

int main(int argc, char *argv[]) {
    if (argc < 5) {
        cout << "Usage: " << argv[0] << " text-file " << " weights-file " << " K " << " patterns-file" << endl;
        cout << "or" << endl;
        cout << "Usage: " << argv[0] << " text-file " << " random " << " K " << " no_of_patterns min_pattern_length max_pattern_length" << endl;
        return -1;
    }
    string second_arg(argv[2]);

    INT K = atoll(argv[3]);
    INT no_of_patterns;
    INT pat_minlen;
    INT pat_maxlen;

    if (second_arg.compare("random") != 0) {
        if (argc < 5) {
            cout << "Usage: " << argv[0] << " text-file " << " weights-file " << " K " << " patterns-file" << endl;
            return -1;
        }
    } else {
        if (argc < 7) {
            cout << "Usage: " << argv[0] << " text-file " << " random " << " K " << " no_of_patterns min_pattern_length max_pattern_length" << endl;
            return -1;
        }

        no_of_patterns = (INT)atoi(argv[4]);
        pat_minlen = (INT)atoi(argv[5]);
        pat_maxlen = (INT)atoi(argv[6]);
    }

    ifstream seq(argv[1], ios::in | ios::binary);
    unsigned char *input_seq_char = NULL;
    unsigned char *sequence;
    INT n = 0;
    char c;

    INT mall_str = display_mallinfo2();
    while (seq.get(c)) {
        if (n == 0 || n % ALLOC_SIZE) input_seq_char = (unsigned char *)realloc(input_seq_char, (n + ALLOC_SIZE) * sizeof(unsigned char));
        input_seq_char[n] = c;
        n++;
    }
    seq.close();

    INT mallinfo_str = display_mallinfo2() - mall_str;
    global_space += mallinfo_str;

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

    INT mall_sa = display_mallinfo2();
    INT *SA = (INT *)malloc((n) * sizeof(INT));
    if ((SA == NULL)) {
        fprintf(stderr, " Error: Cannot allocate memory for SA.\n");
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

    INT *LCP = (INT *)calloc(n, sizeof(INT));
    if ((LCP == NULL)) {
        fprintf(stderr, " Error: Cannot allocate memory for LCP.\n");
        return (0);
    }

    if (LCParray(sequence, n, SA, invSA, LCP) != 1) {
        fprintf(stderr, " Error: LCP computation failed.\n");
        exit(EXIT_FAILURE);
    }

    free(invSA);
    cout << "SA and LCP array are constructed." << endl;

    INT mallinfo_sa_lcp = display_mallinfo2() - mall_sa;
    global_space += mallinfo_sa_lcp;

    B *b = (B *)malloc((n) * sizeof(B));
    if ((b == NULL)) {
        fprintf(stderr, " Error: Cannot allocate memory for b.\n");
        return (0);
    }

    construct_tuples(sequence, n, SA, LCP, b);

    std::chrono::steady_clock::time_point tup_construct_end = std::chrono::steady_clock::now();
    cout << "Tuple construction time: " << std::chrono::duration_cast<std::chrono::milliseconds>(tup_construct_end - tup_construct_begin).count() << "[ms]." << std::endl;

    INT bsize = 0;
    INT tau = n;
    INT L = 0;
    INT w = 0;
    find_longest(b, n, K, tau, L, bsize, w);
    b = (B *)realloc(b, (bsize) * sizeof(B));

    INT mall_ps = display_mallinfo2();
    vector<double> PS;
    if (second_arg.compare("random") == 0) {
        cout << "Random weights generation.\n";
        random_weight_generation(PS, n - 1);
    } else {
        cout << "The weights are read from file.\n";
        weight_from_file(argv[2], PS, n - 1);
    }

    INT mallinfo_ps = display_mallinfo2() - mall_ps;
    global_space += mallinfo_ps;

    std::chrono::steady_clock::time_point util_comp_begin = std::chrono::steady_clock::now();

    cout << "We are computing the utility of the top K = " << K << " frequent substrings." << endl;

    auto H = ankerl::unordered_dense::map<INT, double>();
    substrings_utility(SA, sequence, PS, b, bsize, H);
    free(b);

    cout << "The computed tau is " << tau << "." << endl;
    cout << "The longest frequent pattern is of length L = " << L << "." << endl;
    char longest_substr[L + 1];
    strncpy(longest_substr, (char *)&sequence[SA[w]], L);
    longest_substr[L] = 0;
    cout << "The utilities of the top K frequent substrings are computed." << endl;

    std::chrono::steady_clock::time_point grouping_end = std::chrono::steady_clock::now();

    std::chrono::steady_clock::time_point util_comp_end = std::chrono::steady_clock::now();
    cout << "Utility computation time: " << std::chrono::duration_cast<std::chrono::milliseconds>(util_comp_end - util_comp_begin).count() << "[ms]." << std::endl;

    std::chrono::steady_clock::time_point ind_construct_begin = std::chrono::steady_clock::now();

    INT mall_rmq = display_mallinfo2();
    int_vector<> v(n, 0);
    for (INT i = 0; i < n; i++) v[i] = LCP[i];
    rmq_succinct_sct<> rmq(&v);
    util::clear(v);

    INT mallinfo_rmq = display_mallinfo2() - mall_rmq;
    global_space += mallinfo_rmq;

    std::chrono::steady_clock::time_point ind_construct_end = std::chrono::steady_clock::now();
    cout << "RMQ construction time: " << std::chrono::duration_cast<std::chrono::milliseconds>(ind_construct_end - ind_construct_begin).count() << "[ms]." << std::endl;

    std::chrono::steady_clock::time_point construct_end = std::chrono::steady_clock::now();
    cout << "Total construction time: " << std::chrono::duration_cast<std::chrono::milliseconds>(construct_end - construct_begin).count() << "[ms]." << std::endl;

    vector<vector<unsigned char> > all_patterns;
    if (second_arg.compare("random") == 0) {
        cout << "Random pattern generation" << endl;
        random_pattern_generation(pat_minlen, pat_maxlen, no_of_patterns, all_patterns, n, sequence);
    } else {
        cout << "The patterns are read from file." << endl;
        pattern_from_file(all_patterns, n, argv[4]);
    }

    std::chrono::steady_clock::time_point pt_begin = std::chrono::steady_clock::now();
    for (auto &pattern : all_patterns) {
        INT m = pattern.size();
        INT fp = 0;
        for (INT i = 0; i < m; i++) fp = karp_rabin_hashing::concat(fp, pattern[i], 1);
        double U = 0.0;

        if (H.find(fp) != H.end()) {
            std::chrono::steady_clock::time_point pt_hash_end = std::chrono::steady_clock::now();

        } else {
            pair<INT, INT> interval = pattern_matching(pattern, sequence, SA, LCP, rmq, n, m);
            INT occs = 0;
            if (interval.second >= interval.first) {
                occs = interval.second - interval.first + 1;
                for (int i = 0; i < occs; i++) {
                    if (SA[interval.first + i] == 0)
                        U += PS[m - 1];
                    else
                        U += PS[SA[interval.first + i] + m - 1] - PS[SA[interval.first + i] - 1];
                }
            }
        }
    }

    std::chrono::steady_clock::time_point pt_end = std::chrono::steady_clock::now();
    cout << "====Index size in MBytes: " << (double)global_space / 1000000 << endl;
    cout << "Total pattern matching time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(pt_end - pt_begin).count() << "[ns]." << std::endl;
    cout << "Avg pattern matching time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(pt_end - pt_begin).count() / (double)all_patterns.size() << "[ns]." << std::endl;
    free(sequence);
    free(LCP);
    free(SA);
    return 0;
}
