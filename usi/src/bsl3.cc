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

#include <malloc.h>
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
#include <unordered_map>
#include <vector>

#include "defs.h"
#include "krfp.h"
#include "spacesavingplus.hpp"
#include "unordered_dense.h"
#include "useful_methods.h"
#include "utils.h"

using namespace sdsl;
using namespace std;

int64_t g_global_space = 0;

double get_utility(SpaceSavingPlus &SSP, unsigned char *sequence, INT *SA, INT *LCP, rmq_succinct_sct<> &rmq, INT n, vector<unsigned char> &pattern, vector<double> &PS) {
    INT fp = 0;
    for (INT i = 0; i < pattern.size(); i++) {
        fp = karp_rabin_hashing::concat(fp, pattern[i], 1);
    }

    double util;
    if (SSP.existsKey(fp)) {
        util = SSP.getUtility(fp);
    } else {
        util = cal_utility(pattern, PS, sequence, SA, LCP, rmq, n);
    }

    SSP.processKey(fp, util);

    return util;
}

int main(int argc, char *argv[]) {
    if (argc < 5) {
        cout << "Usage: " << argv[0] << " text-file " << " weights-file " << " K " << " patterns-file" << endl;
        return 1;
    }

    ifstream seq(argv[1], ios::in | ios::binary);
    unsigned char *input_seq_char = NULL;
    INT n = 0;
    char c;

    int64_t mallinfo_seq_begin = display_mallinfo2();

    while (seq.get(c)) {
        if (n == 0 || n % ALLOC_SIZE) input_seq_char = (unsigned char *)realloc(input_seq_char, (n + ALLOC_SIZE) * sizeof(unsigned char));
        input_seq_char[n] = c;
        n++;
    }
    seq.close();

    int64_t mallinfo_seq_end = display_mallinfo2();
    g_global_space = g_global_space + mallinfo_seq_end - mallinfo_seq_begin;

    if (c == '\n')
        input_seq_char[--n] = '\0';
    else
        input_seq_char[n] = '\0';
    unsigned char *sequence = input_seq_char;
    n = strlen((const char *)sequence);
    cout << "Text is of length n = " << n << "." << endl;

    INT K = atoll(argv[3]);

    int64_t mallinfo_PS_begin = display_mallinfo2();
    vector<double> PS;
    ifstream is;
    is.open(argv[2], ios::binary);
    std::istream_iterator<double> start(is), end;
    std::chrono::steady_clock::time_point ps_construct_begin = std::chrono::steady_clock::now();
    for (auto it = start; it != end; ++it) PS.push_back(*it);
    is.close();
    if (PS.size() != n) {
        cout << " Error: Weights size do not match the text size: " << PS.size() << endl;
        return (1);
    }
    for (INT i = 1; i < n; ++i) PS[i] += PS[i - 1];

    int64_t mallinfo_PS_end = display_mallinfo2();
    g_global_space = g_global_space + mallinfo_PS_end - mallinfo_PS_begin;

    std::chrono::steady_clock::time_point ps_construct_end = std::chrono::steady_clock::now();
    cout << "Prefix-sum construction time: " << std::chrono::duration_cast<std::chrono::milliseconds>(ps_construct_end - ps_construct_begin).count() << "[ms]" << std::endl;

    std::chrono::steady_clock::time_point ind_construct_begin = std::chrono::steady_clock::now();
    INT *SA;
    INT *LCP;
    INT *invSA;

    int64_t mallinfo_SA_LCP_rmq_begin = display_mallinfo2();

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

    int64_t mallinfo_SA_LCP_rmq_end = display_mallinfo2();
    g_global_space = g_global_space + mallinfo_SA_LCP_rmq_end - mallinfo_SA_LCP_rmq_begin;

    std::chrono::steady_clock::time_point ind_construct_end = std::chrono::steady_clock::now();
    cout << "Total construction time: " << std::chrono::duration_cast<std::chrono::milliseconds>(ind_construct_end - ind_construct_begin).count() << "[ms]" << std::endl;

    cout << "SA is constructed as the text index." << endl;
    cout << "Total construction time: " << std::chrono::duration_cast<std::chrono::milliseconds>(ind_construct_end - ps_construct_begin).count() << "[ms]" << std::endl;

    SpaceSavingPlus SSP(K);

    INT hash_variable = karp_rabin_hashing::init();
    std::chrono::steady_clock::time_point pt_begin = std::chrono::steady_clock::now();

    ifstream is2;
    is2.open(argv[4], ios::in | ios::binary);
    size_t pattern_num = 0;
    vector<unsigned char> pattern;

    int64_t mallinfo_cache_begin = display_mallinfo2();

    while (is2.read(reinterpret_cast<char *>(&c), 1)) {
        if (c == '\n') {
            if (pattern.empty()) break;

            double U = get_utility(SSP, sequence, SA, LCP, rmq, n, pattern, PS);
            pattern_num++;
            pattern.clear();
        } else
            pattern.push_back((unsigned char)c);
    }
    is2.close();
    pattern.clear();

    int64_t mallinfo_cache_end = display_mallinfo2();
    g_global_space = g_global_space + mallinfo_cache_end - mallinfo_cache_begin;

    cout << endl;
    std::chrono::steady_clock::time_point pt_end = std::chrono::steady_clock::now();

    cout << "Index size in MBytes: " << (double)g_global_space / 1000000 << endl;
    cout << "Total pattern matching time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(pt_end - pt_begin).count() << "[ns]" << std::endl;
    cout << "Avg pattern matching time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(pt_end - pt_begin).count() / (double)pattern_num << "[ns]" << std::endl;

    return 0;
}
