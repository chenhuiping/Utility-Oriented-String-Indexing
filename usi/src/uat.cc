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

#include <chrono>
#include <cmath>
#include <ctime>
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

INT global_space = 0;

unsigned int calculate_utility(unsigned char *sequence, INT n, INT *SA, INT *LCP, rmq_succinct_sct<> &rmq, vector<Util_a *> &util_list, vector<double> PS, ankerl::unordered_dense::map<INT, double> &H,
                               INT K) {
    vector<tuple<INT, double>> H_insert(K);

    for (INT i = 0; i < K; i++) {
        INT start_pos = util_list[i]->start;
        INT ell = util_list[i]->ell;
        vector<unsigned char> pattern;

        for (INT curt_pos = 0; curt_pos < ell; curt_pos++) {
            pattern.emplace_back(sequence[start_pos + curt_pos]);
        }

        pair<INT, INT> interval = pattern_matching(pattern, sequence, SA, LCP, rmq, n, pattern.size());
        INT occs = interval.second - interval.first + 1;
        double U = 0;
        for (INT j = 0; j < occs; j++) {
            INT start = SA[interval.first + j];
            if (start == 0)
                U += PS[util_list[i]->ell - 1];
            else
                U += PS[start + util_list[i]->ell - 1] - PS[start - 1];
        }

        H_insert[i] = make_tuple(util_list[i]->fp, U);
    }

    for (auto &t : H_insert) {
        H[get<0>(t)] = get<1>(t);
    }

    return (1);
}

int main(int argc, char *argv[]) {
    if (argc < 6) {
        cout << "Usage: " << argv[0] << " text-file " << " weights-file " << " K " << " tau " << " patterns-file" << endl;
        cout << "or" << endl;
        cout << "Usage: " << argv[0] << " text-file " << " random " << " K " << " tau " << " no_of_patterns min_pattern_length max_pattern_length" << endl;
        return -1;
    }

    string second_arg(argv[2]);
    INT K = atoll(argv[3]);
    INT no_of_patterns;
    INT pat_minlen;
    INT pat_maxlen;
    INT interval = atoll(argv[4]);

    if (second_arg.compare("random") != 0) {
        if (argc < 6) {
            cout << "Usage: " << argv[0] << " text-file " << " weights-file " << " K " << " tau " << " patterns-file" << endl;
            return -1;
        }
    } else {
        if (argc < 8) {
            cout << "Usage: " << argv[0] << " text-file " << " random " << " K " << " tau " << " no_of_patterns min_pattern_length max_pattern_length" << endl;
            return -1;
        }

        no_of_patterns = (INT)atoi(argv[5]);
        pat_minlen = (INT)atoi(argv[6]);
        pat_maxlen = (INT)atoi(argv[7]);
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

    global_space += display_mallinfo2() - mall_str;

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

        std::sort(SSA.begin(), SSA.end(), lce.lex_less_than());

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

    cout << "Tuple construction time in total: " << tuple_construt_time_total << "[ms]." << std::endl;

    INT mall_ps = display_mallinfo2();

    vector<double> PS;

    if (second_arg.compare("random") == 0) {
        cout << "Random weights generation.\n";
        random_weight_generation(PS, n - 1);
    } else {
        cout << "The weights are read from file.\n";
        weight_from_file(argv[2], PS, n - 1);
    }

    global_space += display_mallinfo2() - mall_ps;

    std::chrono::steady_clock::time_point ind_construct_begin = std::chrono::steady_clock::now();

    INT *SA;
    INT *LCP;
    INT *invSA;
    INT mall_sa = display_mallinfo2();

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
    global_space += display_mallinfo2() - mall_sa;

    std::chrono::steady_clock::time_point ind_construct_end = std::chrono::steady_clock::now();
    cout << "Text index construction time: " << std::chrono::duration_cast<std::chrono::milliseconds>(ind_construct_end - ind_construct_begin).count() << "[ms]." << std::endl;

    std::chrono::steady_clock::time_point util_comp_begin = std::chrono::steady_clock::now();

    cout << "We are computing the utility of the top " << K << " frequent substrings." << endl;
    auto H = ankerl::unordered_dense::map<INT, double>();
    if (util_list.size() < K) {
        K = util_list.size();
    }
    calculate_utility(sequence, n, SA, LCP, rmq, util_list, PS, H, K);

    for (auto &pu : util_list) {
        delete pu;
    }
    vector<Util_a *>().swap(util_list);

    std::chrono::steady_clock::time_point util_comp_end = std::chrono::steady_clock::now();
    cout << "Utility computation time: " << std::chrono::duration_cast<std::chrono::milliseconds>(util_comp_end - util_comp_begin).count() << "[ms]." << std::endl;

    std::chrono::steady_clock::time_point construct_end = std::chrono::steady_clock::now();
    cout << "Total construction time: " << std::chrono::duration_cast<std::chrono::milliseconds>(construct_end - construct_begin).count() << "[ms]." << std::endl;

    cout << endl;

    vector<vector<unsigned char>> all_patterns;
    if (second_arg.compare("random") == 0) {
        cout << "Random pattern generation" << endl;
        random_pattern_generation(pat_minlen, pat_maxlen, no_of_patterns, all_patterns, n, sequence);
    } else {
        cout << "The patterns are read from file." << endl;
        pattern_from_file(all_patterns, n, argv[5]);
    }

    std::chrono::steady_clock::time_point pt_begin = std::chrono::steady_clock::now();


    for (auto &pattern : all_patterns) {
        INT m = pattern.size();
        INT fp = 0;
        for (INT i = 0; i < m; i++) fp = karp_rabin_hashing::concat(fp, pattern[i], 1);

        if (H.find(fp) != H.end())
            ;
        else {
            double U = 0;
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

    free(SA);
    free(LCP);
    free(sequence);

    return 0;
}
