#ifndef __USEFUL_METHODS_INCLUDED
#define __USEFUL_METHODS_INCLUDED

#include <sdsl/suffix_trees.hpp>
#include <tuple>

#include "defs.h"
#include "unordered_dense.h"

int64_t display_mallinfo2(void);

INT lcp(unsigned char *x, INT M, vector<unsigned char> y, INT l, INT a_size, INT w_size);

unsigned int LCParray(unsigned char *text, INT n, INT *SA, INT *ISA, INT *LCP);

pair<INT, INT> pattern_matching(vector<unsigned char> p, unsigned char *T, INT *SA, INT *LCP, sdsl::rmq_succinct_sct<> &rmq, INT n, INT p_size);

double cal_utility(vector<unsigned char> &pattern, vector<double> &PS, unsigned char *sequence, INT *SA, INT *LCP, rmq_succinct_sct<> &rmq, INT n);

struct B {
    INT l = 0;
    INT r = 0;
    INT lcp = 0;
    vector<INT> ch;
};

bool tuples_sorter(B const &lhs, B const &rhs);

unsigned int construct_tuples(unsigned char *seq, INT n, INT *SA, INT *LCP, B *b);
unsigned int construct_tuples(unsigned char *seq, INT n, vector<INT> &SSA, vector<INT> &SLCP, B *b, INT each_sample_K);

unsigned int find_longest(B *b, INT n, INT K, INT &tau, INT &L, INT &bsize, INT &w);

struct Util_a {
    INT fp = 0;
    INT freq = 0;
    uint32_t ell = 0;
    INT start = -1;
};

void gen_util_list(vector<Util_a *> &util_list, unordered_map<INT, Util_a *> &util_index, unsigned char *sequence, vector<INT> &SSA, B *b, INT bsize, INT n, INT K);

struct Util_e {
    INT fp = 0;
    INT l = 0;
    INT r = 0;
    uint32_t ell = 0;
    double util = 0;
};

unsigned int substrings_utility(INT *SA, unsigned char *sequence, vector<double> PS, B *b, INT bsize, ankerl::unordered_dense::map<INT, double> &H, string inputfile, INT K);
unsigned int substrings_utility(INT *SA, unsigned char *sequence, vector<double> PS, B *b, INT bsize, ankerl::unordered_dense::map<INT, double> &H);

void random_weight_generation(vector<double> &PS, INT n);
void weight_from_file(char *weights_fname, vector<double> &PS, INT n);

void random_pattern_generation(INT pat_minlen, INT pat_maxlen, INT no_of_patterns, vector<vector<unsigned char>> &all_patterns, INT n, unsigned char *sequence);
void pattern_from_file(vector<vector<unsigned char>> &all_patterns, INT n, char *patterns_fname);

#endif