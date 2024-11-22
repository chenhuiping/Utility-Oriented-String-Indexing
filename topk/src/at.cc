/**
    USI: useful string indexing
    Copyright (C) 2024 Huiping Chen.

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

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <execution>
#include <iostream>
#include <sdsl/suffix_trees.hpp>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "bitv.hpp"
#include "defs.h"
#include "includes.hpp"
#include "krfp.h"
#include "rk_lce.hpp"
#include "rk_lce_bin.hpp"
#include "unordered_dense.h"
#include "useful_methods.h"
#include "utils.h"

using namespace sdsl;
using namespace std;
using namespace rklce;

int main(int argc, char *argv[]) {
    if (argc != 5) {
        cout << "Usage: " << argv[0] << " text-file " << " K " << " tau " << endl;
        return -1;
    }

    INT K = atoll(argv[2]);
    INT interval = atoll(argv[3]);

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

    if (interval >= n - 1) {
        cout << "Please find an interval that smaller then the length of the string " << n - 1 << endl;
        return 0;
    }

    cout << "Text is of length n = " << n - 1 << "." << endl;
    cout << "Parameters: K = " << K << ", tau = " << interval << endl;
    INT threads = omp_get_num_procs();
    cout << "Threads: " << threads << endl;

    std::chrono::steady_clock::time_point construct_begin = std::chrono::steady_clock::now();

    std::chrono::steady_clock::time_point rl_lce_start = std::chrono::steady_clock::now();
    auto lce = rk_lce(string(argv[1]));
    std::chrono::steady_clock::time_point rl_lce_end = std::chrono::steady_clock::now();
    cout << "LCE data structure construction time: " << std::chrono::duration_cast<std::chrono::milliseconds>(rl_lce_end - rl_lce_start).count() << "[ms]." << std::endl;

    karp_rabin_hashing::init();
    unordered_map<INT, Util_a *> util_index;
    vector<Util_a *> util_list;

    INT num_iter = 0;
    double tuple_construt_time_total = 0.0;

    while (num_iter < interval) {
        std::chrono::steady_clock::time_point tup_construct_begin = std::chrono::steady_clock::now();
        vector<INT> SSA;
        SSA = vector<INT>();

        INT curt_pos = num_iter;
        INT num = 0;

        while (curt_pos < n) {
            SSA.push_back(curt_pos);
            num++;
            curt_pos += interval;
        }

        SSA.resize(num);

        std::sort(std::execution::par, SSA.begin(), SSA.end(), lce.lex_less_than());

        vector<INT> SLCP(SSA.size() + 1, 0);

        for (INT i = 1; i < SSA.size(); i++) {
            uint64_t lcp = lce.LCE(SSA[i - 1], SSA[i]);
            SLCP[i] = lcp;
        }

        auto b = new B[SSA.size() + 1];

        if (b == nullptr) {
            fprintf(stderr, " Error: Cannot allocate memory for b.\n");
            return (0);
        }

        construct_tuples(sequence, SSA.size() + 1, SSA, SLCP, b, K);

        std::chrono::steady_clock::time_point tup_construct_end = std::chrono::steady_clock::now();

        tuple_construt_time_total += std::chrono::duration_cast<std::chrono::milliseconds>(tup_construct_end - tup_construct_begin).count();

        INT bsize = 0;
        INT tau = n;
        INT L = 0;
        INT w = 0;
        find_longest(b, num + 1, K, tau, L, bsize, w);

        auto b_beta = new B[bsize];
        if (b_beta == nullptr) {
            fprintf(stderr, " Error: Cannot allocate memory for b.\n");
            return (0);
        }
        copy(b, b + bsize, b_beta);
        delete[] b;
        b = b_beta;

        gen_util_list(util_list, util_index, sequence, SSA, b, bsize, n, K);

        num_iter += 1;
    }
    util_index.clear();

    if (util_list.size() < K) {
        cout << "final output size is " << util_list.size() << ", smaller than " << K << endl;
    }
    cout << "Tuple construction time in total: " << tuple_construt_time_total << "[ms]." << std::endl;

    INT *SA;
    INT *LCP;
    INT *invSA;

    SA = (INT *)malloc((n) * sizeof(INT));
    if ((SA == NULL)) {
        fprintf(stderr, " Error: Cannot allocate memory for SA.\n");
        return (0);
    }

    if (divsufsort64(sequence, SA, n) != 0) {
        fprintf(stderr, " Error: SA computation failed.\n");
        exit(EXIT_FAILURE);
    }

    invSA = (INT *)malloc((n) * sizeof(INT));

    if ((invSA == NULL)) {
        fprintf(stderr, " Error: Cannot allocate memory for invSA.\n");
        return (0);
    }

    for (INT i = 0; i < n; i++) {
        invSA[SA[i]] = i;
    }

    LCP = (INT *)malloc((n) * sizeof(INT));
    if ((LCP == NULL)) {
        fprintf(stderr, " Error: Cannot allocate memory for LCP.\n");
        return (0);
    }

    LCParray(sequence, n, SA, invSA, LCP);

    free(invSA);

    int_vector<> v(n, 0);
    for (INT i = 0; i < n; i++) v[i] = LCP[i];
    rmq_succinct_sct<> rmq(&v);
    util::clear(v);

    unordered_map<string, INT> pattern_act_freq;
    vector<pair<string, INT>> pattern_freq_to_add(util_list.size());


    for (INT i = 0; i < util_list.size(); i++) {
        INT start_pos = util_list[i]->start;
        INT ell = util_list[i]->ell;
        vector<unsigned char> pattern;
        string temp_str = "";

        for (INT curt_pos = 0; curt_pos < ell; curt_pos++) {
            pattern.emplace_back(sequence[start_pos + curt_pos]);
            temp_str += sequence[start_pos + curt_pos];
        }

        pair<INT, INT> interval = pattern_matching(pattern, sequence, SA, LCP, rmq, n, pattern.size());
        INT occs = interval.second - interval.first + 1;

        pattern_freq_to_add[i] = {temp_str, occs};
    }

    for (const auto &p : pattern_freq_to_add) {
        pattern_act_freq.insert(p);
    }

    free(sequence);

    return 0;
}
