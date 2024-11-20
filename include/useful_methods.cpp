#include "useful_methods.h"

#include <malloc.h>
#include <omp.h>

#include "krfp.h"

int64_t display_mallinfo2(void) {
    struct mallinfo2 mi;
    mi = mallinfo2();
    return (int64_t)(mi.uordblks + mi.hblkhd);
}

INT lcp(unsigned char *x, INT M, vector<unsigned char> y, INT l, INT a_size, INT w_size) {
    INT xx = a_size;
    if (M >= xx) return 0;
    INT yy = w_size;
    if (l >= yy) return 0;

    INT i = 0;
    while ((M + i < xx) && (l + i < yy)) {
        if (x[M + i] != y[l + i]) break;
        i++;
    }
    return i;
}

unsigned int LCParray(unsigned char *text, INT n, INT *SA, INT *ISA, INT *LCP) {
    INT i = 0, j = 0;

    LCP[0] = 0;
    for (i = 0; i < n; i++)
        if (ISA[i] != 0) {
            if (i == 0)
                j = 0;
            else
                j = (LCP[ISA[i - 1]] >= 2) ? LCP[ISA[i - 1]] - 1 : 0;
            while (text[i + j] == text[SA[ISA[i] - 1] + j]) j++;
            LCP[ISA[i]] = j;
        }

    return (1);
}

pair<INT, INT> pattern_matching(vector<unsigned char> p, unsigned char *T, INT *SA, INT *LCP, rmq_succinct_sct<> &rmq, INT n, INT p_size) {
    INT m = p_size;
    INT N = n;
    INT d = -1;
    INT ld = 0;
    INT f = n;
    INT lf = 0;

    pair<INT, INT> interval;

    while (d + 1 < f) {
        INT i = (d + f) / 2;

        INT lcpif;

        if (f == n)
            lcpif = 0;
        else
            lcpif = LCP[rmq(i + 1, f)];

        INT lcpdi;

        if (i == n)
            lcpdi = 0;
        else
            lcpdi = LCP[rmq(d + 1, i)];

        if ((ld <= lcpif) && (lcpif < lf)) {
            d = i;
            ld = lcpif;
        } else if ((ld <= lf) && (lf < lcpif))
            f = i;
        else if ((lf <= lcpdi) && (lcpdi < ld)) {
            f = i;
            lf = lcpdi;
        } else if ((lf <= ld) && (ld < lcpdi))
            d = i;
        else {
            INT l = std::max(ld, lf);
            l = l + lcp(T, SA[i] + l, p, l, n, p_size);
            if (l == m) {
                INT e = i;
                while (d + 1 < e) {
                    INT j = (d + e) / 2;

                    INT lcpje;

                    if (e == n)
                        lcpje = 0;
                    else
                        lcpje = LCP[rmq(j + 1, e)];

                    if (lcpje < m)
                        d = j;
                    else
                        e = j;
                }

                INT lcpde;

                if (e == n)
                    lcpde = 0;
                else
                    lcpde = LCP[rmq(d + 1, e)];

                if (lcpde >= m) d = std::max(d - 1, (INT)-1);

                e = i;
                while (e + 1 < f) {
                    INT j = (e + f) / 2;

                    INT lcpej;

                    if (j == n)
                        lcpej = 0;
                    else
                        lcpej = LCP[rmq(e + 1, j)];

                    if (lcpej < m)
                        f = j;
                    else
                        e = j;
                }

                INT lcpef;

                if (f == n)
                    lcpef = 0;
                else
                    lcpef = LCP[rmq(e + 1, f)];

                if (lcpef >= m) f = std::min(f + 1, n);

                interval.first = d + 1;
                interval.second = f - 1;
                return interval;

            } else if ((l == N - SA[i]) || ((SA[i] + l < N) && (l != m) && (T[SA[i] + l] < p[l]))) {
                d = i;
                ld = l;
            } else {
                f = i;
                lf = l;
            }
        }
    }

    interval.first = d + 1;
    interval.second = f - 1;
    return interval;
}

double cal_utility(vector<unsigned char> &pattern, vector<double> &PS, unsigned char *sequence, INT *SA, INT *LCP, rmq_succinct_sct<> &rmq, INT n) {
    INT m = pattern.size();
    double U = 0;
    pair<INT, INT> interval = pattern_matching(pattern, sequence, SA, LCP, rmq, n, pattern.size());
    if (interval.first > interval.second) {
        return U;
    }

    if (interval.second >= interval.first) {
        INT occs = interval.second - interval.first + 1;
        for (int i = 0; i < occs; i++) {
            if (SA[interval.first + i] == 0)
                U += PS[m - 1];
            else
                U += PS[SA[interval.first + i] + m - 1] - PS[SA[interval.first + i] - 1];
        }
    }

    return U;
}

bool tuples_sorter(B const &lhs, B const &rhs) { return (lhs.r - lhs.l) > (rhs.r - rhs.l); }

unsigned int construct_tuples(unsigned char *seq, INT n, INT *SA, INT *LCP, B *b) {
    stack<INT> st;

    b[0].l = 0;
    b[0].r = 0;
    b[0].ch.clear();
    st.push(0);

    INT x = -1;

    for (INT i = 1; i < n; i++) {
        INT l = i - 1;
        while (LCP[i] < LCP[st.top()]) {
            x = st.top();
            st.pop();

            b[x].r = i - 1;

            for (INT v : b[x].ch) {
                vector<INT>().swap(b[v].ch);
                b[v].ch.push_back(b[x].lcp);
            }
            l = b[x].l;

            if (LCP[i] <= LCP[st.top()]) {
                b[st.top()].ch.push_back(x);
                x = -1;
            }
        }
        if (LCP[i] > LCP[st.top()]) {
            b[i].lcp = LCP[i];
            b[i].l = l;
            b[i].ch.clear();
            st.push(i);
            if (~x) {
                b[i].ch.push_back(x);
                x = -1;
            }
        }
    }

    for (INT v : b[0].ch) {
        vector<INT>().swap(b[v].ch);
        b[v].ch.push_back(b[0].lcp);
    }

    cout << "Tuples corresponding to internal nodes are constructed." << endl;

    sort(b, b + n, &tuples_sorter);

    cout << "Tuples are sorted." << endl;

    return (1);
}

unsigned int construct_tuples(unsigned char *seq, INT n, vector<INT> &SSA, vector<INT> &SLCP, B *b, INT each_sample_K) {
    stack<INT> st;
    st.push(0);
    INT x = -1;

    for (INT i = 1; i < n; i++) {
        INT l = i - 1;

        while (SLCP[i] < SLCP[st.top()]) {
            x = st.top();
            st.pop();

            b[x].r = i - 1;
            for (INT v : b[x].ch) {
                vector<INT>().swap(b[v].ch);
                b[v].ch.push_back(b[x].lcp);
            }
            l = b[x].l;

            if (SLCP[i] <= SLCP[st.top()]) {
                b[st.top()].ch.push_back(x);
                x = -1;
            }
        }
        if (SLCP[i] > SLCP[st.top()]) {
            b[i].lcp = SLCP[i];
            b[i].l = l;
            b[i].ch.clear();
            st.push(i);

            if (~x) {
                b[i].ch.push_back(x);
                x = -1;
            }
        }
    }

    for (INT v : b[0].ch) {
        vector<INT>().swap(b[v].ch);
        b[v].ch.push_back(b[0].lcp);
    }

    if (each_sample_K < (SSA.size() + 1)) {
        partial_sort(b, b + each_sample_K, b + SSA.size() + 1, &tuples_sorter);
    } else {
        sort(b, b + SSA.size(), &tuples_sorter);
    }

    return (1);
}

unsigned int find_longest(B *b, INT n, INT K, INT &tau, INT &L, INT &bsize, INT &w) {
    INT k = 0;
    for (INT i = 0; i < n; i++) {
        if (k >= K) break;
        if (b[i].lcp > 0) {
            INT s = b[i].r - b[i].l + 1;
            uint32_t d = b[i].lcp;
            INT pd = b[i].ch[0];
            if (s < tau) tau = s;
            if (d > L) {
                L = d;
                w = b[i].l;
            }
            k += d - pd;
        }
        bsize++;
    }
    return (1);
}

bool cmp_util_list(const Util_a *a, const Util_a *b) { return a->freq > b->freq; }

void gen_util_list(vector<Util_a *> &util_list, unordered_map<INT, Util_a *> &util_index, unsigned char *sequence, vector<INT> &SSA, B *b, INT bsize, INT n, INT K) {
    for (INT i = 0; i < bsize; i++) {
        if (b[i].lcp == 0) continue;

        uint32_t d = b[i].lcp;
        INT pd = b[i].ch[0];
        INT freq_curt = b[i].r - b[i].l + 1;

        INT start_pos = SSA[b[i].l];

        uint64_t fp_curt = static_cast<INT>(sequence[start_pos]);

        for (INT m = 1; m < pd; m++) {
            fp_curt = karp_rabin_hashing::concat(fp_curt, static_cast<uint64_t>(sequence[start_pos + m]), 1);
        }

        for (INT ell = pd + 1; ell <= d; ell++) {
            if (ell > 1) {
                fp_curt = karp_rabin_hashing::concat(fp_curt, static_cast<uint64_t>(sequence[start_pos + ell - 1]), 1);
            }

            if (util_index.find(fp_curt) != util_index.end()) {
                util_index[fp_curt]->freq += freq_curt;
                continue;
            }
            auto pu = new Util_a;
            pu->fp = fp_curt;

            pu->ell = ell;
            pu->start = start_pos;
            pu->freq = freq_curt;
            util_list.emplace_back(pu);
            util_index[pu->fp] = pu;
        }
    }
    delete[] b;

    if (util_list.size() < K) {
        return;
    }

    sort(util_list.begin(), util_list.end(), cmp_util_list);

    for (INT i = K; i < util_list.size(); ++i) {
        auto pu = util_list[i];
        util_index.erase(pu->fp);
        delete pu;
    }
    util_list.resize(K);
}

void write_freq_pattern_file(vector<Util_e> &utilities, string filename, INT *SA, unsigned char *sequence) {
    ofstream of(filename, ios::out);
    if (!of) {
        cout << "error opening file: " << filename << endl;
        exit(1);
    }

    for (auto &u : utilities) {
        INT start = SA[u.l];
        for (INT i = start; i < u.ell + start; i++) {
            of << sequence[i];
        }

        of << " " << u.r - u.l + 1 << endl;
    }
    of.close();
}

unsigned int substrings_utility(B *b, INT bsize, ankerl::unordered_dense::map<INT, double> &H) {
    vector<Util_e> utilities;
    for (INT i = 0; i < bsize; i++) {
        if (b[i].lcp == 0) continue;

        uint32_t d = b[i].lcp;
        INT pd = b[i].ch[0];

        for (INT ell = pd + 1; ell <= d; ell++) {
            Util_e u;
            u.fp = 0;
            u.ell = ell;
            u.l = b[i].l;
            u.r = b[i].r;
            u.util = 0;
            utilities.push_back(u);
        }
    }

    for (INT i = 0; i < utilities.size(); i++) H.insert({utilities[i].fp, utilities[i].util});

    return (1);
}

unsigned int substrings_utility(INT *SA, unsigned char *sequence, vector<double> PS, B *b, INT bsize, ankerl::unordered_dense::map<INT, double> &H) {
    vector<Util_e> utilities;
    for (INT i = 0; i < bsize; i++) {
        if (b[i].lcp == 0) continue;
        INT d = b[i].lcp;
        INT pd = b[i].ch[0];
        for (INT ell = pd + 1; ell <= d; ell++) {
            Util_e u;
            u.fp = 0;
            u.ell = ell;
            u.l = b[i].l;
            u.r = b[i].r;
            u.util = 0;
            utilities.push_back(u);
        }
    }

    std::chrono::steady_clock::time_point fp_construct_begin = std::chrono::steady_clock::now();
    INT n = PS.size();
    INT hash_variable = karp_rabin_hashing::init();
    INT *FP = (INT *)calloc(n, sizeof(INT));
    FP[0] = karp_rabin_hashing::concat(0, sequence[0], 1);
    for (INT i = 1; i < n; ++i) FP[i] = karp_rabin_hashing::concat(FP[i - 1], sequence[i], 1);
    std::chrono::steady_clock::time_point fp_construct_end = std::chrono::steady_clock::now();
    cout << "FP table construction time: " << std::chrono::duration_cast<std::chrono::milliseconds>(fp_construct_end - fp_construct_begin).count() << "[ms]." << std::endl;

    INT threads = omp_get_num_procs();
    cout << threads << " threads will handle " << utilities.size() << " substrings." << endl;

#pragma omp parallel for schedule(dynamic)
    for (INT i = 0; i < utilities.size(); i++) {
        INT fp;
        INT start = SA[utilities[i].l];
        if (start == 0)
            fp = FP[utilities[i].ell - 1];
        else
            fp = karp_rabin_hashing::subtract(FP[start + utilities[i].ell - 1], FP[start - 1], utilities[i].ell);
        double U = 0;
        for (INT j = utilities[i].l; j <= utilities[i].r; j++) {
            INT start = SA[j];
            if (start == 0)
                U += PS[utilities[i].ell - 1];
            else
                U += PS[start + utilities[i].ell - 1] - PS[start - 1];
        }
        utilities[i].fp = fp;
        utilities[i].util = U;
    }

    free(FP);
    for (INT i = 0; i < utilities.size(); i++) H.insert({utilities[i].fp, utilities[i].util});
    return (1);
}

void random_weight_generation(vector<double> &PS, INT n) {
    vector<double> weight_domain{0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};

    mt19937 random_engine(12345);
    uniform_int_distribution<INT> dist(0, weight_domain.size() - 1);
    INT no_of_zeros = 0;

    for (INT i = 0; i < n; ++i) {
        double d = weight_domain.at(dist(random_engine));
        PS.push_back(d);
    }

    cout << "Number of generated weights: " << n << endl;

    std::chrono::steady_clock::time_point ps_construct_begin = std::chrono::steady_clock::now();
    for (INT i = 1; i < n; ++i) PS[i] += PS[i - 1];
    std::chrono::steady_clock::time_point ps_construct_end = std::chrono::steady_clock::now();
    cout << "Prefix-sum construction time: " << std::chrono::duration_cast<std::chrono::milliseconds>(ps_construct_end - ps_construct_begin).count() << "[ms]." << std::endl;
}

void weight_from_file(char *weights_fname, vector<double> &PS, INT n) {
    char c;
    ifstream is;
    is.open(weights_fname, ios::binary);
    std::istream_iterator<double> start(is), end;

    for (auto it = start; it != end; ++it) PS.push_back(*it);
    is.close();
    if (PS.size() != n) {
        cout << " Error: Weights size do not match the text size: " << PS.size() << endl;
        exit(-1);
    }
    std::chrono::steady_clock::time_point ps_construct_begin = std::chrono::steady_clock::now();
    for (INT i = 1; i < n; ++i) PS[i] += PS[i - 1];
    std::chrono::steady_clock::time_point ps_construct_end = std::chrono::steady_clock::now();
    cout << "Prefix-sum construction time: " << std::chrono::duration_cast<std::chrono::milliseconds>(ps_construct_end - ps_construct_begin).count() << "[ms]." << std::endl;
}

void random_pattern_generation(INT pat_minlen, INT pat_maxlen, INT no_of_patterns, vector<vector<unsigned char>> &all_patterns, INT n, unsigned char *sequence) {
    mt19937 random_engine(12345);

    std::uniform_int_distribution<INT> dist2(pat_minlen, pat_maxlen);
    std::uniform_int_distribution<INT> dist3((INT)0, max((INT)0, (INT)(n - 1 - pat_maxlen)));
    for (INT i = 0; i < no_of_patterns; ++i) {
        INT r1 = dist3(random_engine);
        INT r2 = dist2(random_engine);

        vector<unsigned char> pattern;
        for (INT j = 0; j < r2; ++j) {
            if (r1 + j < n) pattern.push_back((unsigned char)sequence[r1 + j]);
        }
        all_patterns.push_back(pattern);
    }
    cout << "The number of patterns is: " << all_patterns.size() << endl;
}

void pattern_from_file(vector<vector<unsigned char>> &all_patterns, INT n, char *patterns_fname) {
    ifstream is2;
    is2.open(patterns_fname, ios::in | ios::binary);

    char c;
    vector<unsigned char> pattern;
    while (is2.read(reinterpret_cast<char *>(&c), 1)) {
        if (c == '\n') {
            if (pattern.empty()) break;
            all_patterns.push_back(pattern);
            pattern.clear();
        } else
            pattern.push_back((unsigned char)c);
    }
    is2.close();
    pattern.clear();
    cout << "The number of patterns is: " << all_patterns.size() << endl;
}