/**
    SPACESAVINGSKETCH: like SpaceSavingSketch, except that it replaces the hash map F with a count-min-sketch (as in HeavyKeeper)
    Copyright (C) 2024 Roberto Grossi and Veronica Guerrini.

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

#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <utility>
#include <vector>

#include "unordered_dense.h"

using namespace std;

class BitVectorNear {
   public:
    BitVectorNear(int64_t _K) {
        K = _K;
        pacK = _K / 8 + 1;
        B.reserve(pacK + 1);
        for (auto i = 0; i <= pacK; i++) B[i] = 0;
        B[_K >> 3] = mask[_K & 0x07];
        n1 = 1;
        nCalls = 0;
        nOps = 0;
    }

    uint8_t get(int64_t i) {
        if (i < 0 || i > K) throw runtime_error("get index out of range");
        return B[i >> 3] & mask[i & 0x07] ? 1 : 0;
    }

    void set(int64_t i) {
        if (i < 0 || i >= K) throw runtime_error("set index out of range");
        if (get(i)) return;
        B[i >> 3] |= mask[i & 0x07];
        n1++;
    }
    void reset(int64_t i) {
        if (i < 0 || i >= K) throw runtime_error("reset index out of range");
        if (!get(i)) return;
        B[i >> 3] &= ~(mask[i & 0x07]);
        n1--;
    }

    int64_t nearestOneRight(int64_t i) {
        if (i < 0 || i >= K) throw runtime_error("rank index out of range");
        nCalls++;
        nOps++;
        auto pos = i + 1;
        while ((pos >> 3) == (i >> 3)) {
            nOps++;
            if (get(pos)) return pos;
            pos++;
        }
        auto pos3 = pos >> 3;
        while (!B[pos3]) {
            pos3++;
            nOps++;
        }
        pos = pos3 << 3;
        while (!get(pos)) {
            pos++;
            nOps++;
        }
        return pos;
    }

    void printB() {
        cout << K << ": ";
        for (auto i = 0; i <= K; i++) cout << (int)(get(i));
        cout << endl << flush;
    }

    void printStat() { cout << (nOps + 0.0) / nCalls << " average ops per NearestOneRight call" << endl; }

   private:
    uint8_t mask[8] = {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80};
    vector<uint8_t> B;
    int64_t n1 = 1;
    int64_t K = 0;
    int64_t pacK = 0;
    int64_t nCalls = 0;
    int64_t nOps = 0;
};

class SpaceSavingSketch {
   public:
    SpaceSavingSketch(int64_t _M2, int64_t _K) : start(_K), K(_K), min_pos(_K), M2(_M2), gen(rd()), dis(0, std::numeric_limits<uint64_t>::max()) {
        sorted_frequency = new node[_K];
        index.clear();
        if (_M2 > MAX_NCOLS) throw runtime_error("please increase MAX_NCOLS");
        for (int64_t i = 0; i < MAX_NROWS; i++) {
            CMS_HK[i] = new snode[_M2];
            for (int64_t j = 0; j < _M2; j++) CMS_HK[i][j].C = CMS_HK[i][j].FP = 0;
        }
    }

    ~SpaceSavingSketch() {
        delete[] sorted_frequency;
        for (int64_t i = 0; i < MAX_NROWS; i++) {
            delete[] CMS_HK[i];
        }
    }

    bool existsKey(uint64_t key) { return index.find(key) != index.end(); }

    double getUtility(uint64_t key) {
        if (existsKey(key)) {
            auto fu = index.at(key);
            return fu.second;
        } else {
            throw runtime_error("Key not found in the top-K");
        }
    }

    void processKey(uint64_t key, double& utility) {
        bool key_is_in_topk = existsKey(key);
        bool success = false;
        int64_t maxv = 0;
        for (int64_t j = 0; j < MAX_NROWS; j++) {
            uint64_t Hsh = UniversalHash(key, j) % M2;
            if (CMS_HK[j][Hsh].FP == key) {
                CMS_HK[j][Hsh].C++;
                maxv = max(maxv, CMS_HK[j][Hsh].C);
            } else {
                uint64_t modval = std::numeric_limits<uint64_t>::max();
                double powval = pow(HK_b, CMS_HK[j][Hsh].C);
                if (modval > powval) modval = powval;
                if (!(myrand64() % modval)) {
                    CMS_HK[j][Hsh].C--;
                    if (CMS_HK[j][Hsh].C <= 0) {
                        CMS_HK[j][Hsh].FP = key;
                        CMS_HK[j][Hsh].C = 1;
                        maxv = max(maxv, static_cast<int64_t>(1));
                    }
                }
            }
        }
        if (!key_is_in_topk) {
            if (maxv == 1 + getMin() || getSize() < K) {
                update(key, maxv, utility);
                success = true;
            }
        } else {
            if (maxv > getFrequency(key)) {
                if (maxv > 1 + getFrequency(key)) {
                    throw runtime_error("too large frequency");
                }
                update(key, maxv, utility);
                success = true;
            }
        }
    }

    void printTopK() {
        auto k = getSize();
        for (auto i = 0; i < k; i++) {
            auto key = get(i);
            auto fu = index.at(key);
            cout << key << " " << sorted_frequency[fu.first].second << " " << fu.second << endl;
        }
        cout << endl;
    }

    void printStat() { start.printStat(); }

   private:
    int64_t K;
    int64_t min_pos;

    ankerl::unordered_dense::map<uint64_t, pair<int64_t, double>> index;

    struct node {
        uint64_t first;
        int64_t second;
    };
    node* sorted_frequency;
    BitVectorNear start;

    uint64_t get(int64_t i) { return sorted_frequency[K - i - 1].first; }

    int64_t getMin() {
        if (min_pos == K) return 0;
        return sorted_frequency[min_pos].second;
    }

    int64_t getSize() { return K - min_pos; }

    int64_t getSizeHash() { return index.size(); }

    int64_t getFrequency(uint64_t key) {
        if (existsKey(key)) {
            auto pu = index.at(key);
            return sorted_frequency[pu.first].second;

        } else {
            throw runtime_error("Key not found in the top-K");
        }
    }

    void update(uint64_t key, int64_t frequency, double& utility) {
        auto p = 0;
        if (!existsKey(key)) {
            if (min_pos > 0) {
                p = --min_pos;
            } else {
                if (frequency != getMin() + 1) return;
                p = 0;
                index.erase(sorted_frequency[0].first);
            }
            index[key].first = p;
            index[key].second = utility;
            sorted_frequency[min_pos].first = key;
            sorted_frequency[min_pos].second = frequency - 1;
        } else {
            p = index[key].first;
            utility = index[key].second;
        }
        increment(p);
    }

    void increment(int64_t p) {
        sorted_frequency[p].second++;
        if (p < K - 1 && sorted_frequency[p].second > sorted_frequency[p + 1].second) {
            auto q = start.nearestOneRight(p) - 1;
            if (p < q) {
                index[sorted_frequency[p].first].first = q;
                index[sorted_frequency[q].first].first = p;
                auto tmp = sorted_frequency[p];
                sorted_frequency[p] = sorted_frequency[q];
                sorted_frequency[q] = tmp;
            }
            p = q;
        }
        start.set(p);
        if (p < K - 1 && sorted_frequency[p].second == sorted_frequency[p + 1].second) start.reset(p + 1);
    }

    static constexpr int64_t MAX_NROWS = 2;
    static constexpr int64_t MAX_NCOLS = 268435456;
    static constexpr double HK_b = 1.08;

    int64_t M2;
    struct snode {
        int64_t C;
        uint64_t FP;
    };
    snode* CMS_HK[MAX_NROWS];

    std::random_device rd;
    std::mt19937_64 gen;
    std::uniform_int_distribution<uint64_t> dis;

    uint64_t myrand64() { return dis(gen); }

    uint64_t UniversalHash(const uint64_t x, int64_t j) {
        if (j == 0) {
            constexpr uint64_t a1 = 0x65d200ce55b19ad7UL;
            constexpr uint64_t a2 = 0x4f2162926e40c299UL;
            return (static_cast<uint64_t>((x ^ (x >> 33U)) * a1) & 0xffffffff00000000UL) | (static_cast<uint64_t>((x ^ (x >> 33U)) * a2) >> 32U);
        } else {
            constexpr uint64_t a3 = 0x68b665e6872bd1f5UL;
            constexpr uint64_t a4 = 0xb6cfcf9d79b51db3UL;
            return (static_cast<uint64_t>((x ^ (x >> 33U)) * a3) & 0xffffffff00000000UL) | (static_cast<uint64_t>((x ^ (x >> 33U)) * a4) >> 32U);
        }
    }
};
