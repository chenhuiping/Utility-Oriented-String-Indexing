/**
    SPACESAVINGPLUS: see below
    Copyright (C) 2024 Roberto Grossi.

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

class SpaceSavingPlus {
   public:
    SpaceSavingPlus(int64_t _K) : start(_K), K(_K), min_pos(_K) {
        sorted_frequency = new node[_K];
        index.clear();
        F.clear();
    }
    ~SpaceSavingPlus() { delete[] sorted_frequency; }

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
        auto f = increaseFrequency(key);
        update(key, f, utility);
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
    ankerl::unordered_dense::map<uint64_t, int64_t> F;

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

    int64_t increaseFrequency(uint64_t key) {
        F[key]++;
        return F[key];
    }

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
};
