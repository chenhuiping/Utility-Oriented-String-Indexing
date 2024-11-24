/**
    MINHEAPDICT: see below
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

// ---> strings replaced with Karp-Rabin fingerprints uint64_t <---

#pragma once

#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include "bitvectornear.hpp"
#include "unordered_dense.h"  

using namespace std;

/*
    Min-heap that is also a hash table, using bitvectornear to handle priorities
    Optimized to use along with heavykeeper in place of ssummary.h
*/
class MinHeapDictionary
{ 
public:
    MinHeapDictionary(int64_t _K) : start(_K), K(_K), min_pos(_K)
    {
        sorted_priority = new node[_K];
        index.clear();
    }
    ~MinHeapDictionary(){
        delete[] sorted_priority;
    }

    // Check if the string exists in the dictionary
    bool exists(const uint64_t key) const
    {
        return index.find(key) != index.end();
    }

    // Get priority of a string, if it exists
    int64_t getPriority(const uint64_t key) const
    {
        if (exists(key))
        {
            return sorted_priority[index.at(key)].second;
        }
        else
        {
            throw runtime_error("Key not found");
        }
    }

    // Change the priority of an existing string
    void update(const uint64_t key, int64_t priority)
    {
        auto p = 0;
        if (!exists(key)){
            if (min_pos > 0){ // room for key
                // by construction, priority is the smallest (with ties)
                index[key] = --min_pos; // we fill sorted_priority from right to left
                sorted_priority[min_pos].first = key;
                sorted_priority[min_pos].second = priority - 1; // priority will be incremented with increment(p)
                p = min_pos;
            } else {
                // replace the pair with smallest priority
                index.erase(sorted_priority[0].first);
                index[key] = 0;
                sorted_priority[0].first = key;
                p = 0;
            }
        }
        else 
        {
            p = index[key];
        }
        increment(p);
    }

    void extractTop(ankerl::unordered_dense::set<uint64_t> &A)
    {
        A.clear();
        auto k = getSize();
        for (auto i = 0; i < k; i++)
            A.insert(get(i));
    }

    uint64_t get(int64_t i)
    {
        return sorted_priority[K - i - 1].first;
    }

    int64_t getMin()
    {
        if (min_pos == K)
            return 0;
        // throw runtime_error("Key not found"); // Check if heap is empty after cleanup
        return sorted_priority[min_pos].second;
    }

    int64_t getSize()
    {
        return K - min_pos;
    }

    int64_t getSize2()
    {
        return index.size();
    }
    
    void printStat(){
        start.printStat();
    }

private:
    struct node{
        uint64_t first;
        int64_t second;
    };
    node * sorted_priority; // (string, its estimated frequency so far) sorted by priority in non-incresing order, filled from right to left
    BitVectorNear start; // start of each segment of equal priorities in sorted_priority: start.B[i] = 1 iff i = 0 or sorted_priority[i-1].second != sorted_priority[i].second
    ankerl::unordered_dense::map<uint64_t,int64_t> index;  // map each string (at most K) to its position in sorted_priority
    int64_t K;
    int64_t min_pos;

    void printPriority(){
        for (auto i = min_pos; i < K; i++)
            cout << sorted_priority[i].second << " ";
        cout << endl << flush;
    }

    void increment(int64_t p){ // sorted_priority[p].second++ and update all the involved data structures
        // cout << "p = " << p << endl << flush;
        sorted_priority[p].second++; // priority
        if (p < K-1 && sorted_priority[p].second > sorted_priority[p+1].second ){
            auto q = start.nearestOneRight(p) - 1; // last position of equal priority pair
            // cout << "q = " << q << endl << flush;
            if (p < q){
                // swap pairs at p and q (last position with same priority before incrementing)
                index[sorted_priority[p].first] = q;
                index[sorted_priority[q].first] = p;
                auto tmp = sorted_priority[p];
                sorted_priority[p] = sorted_priority[q];
                sorted_priority[q] = tmp;
            }
            p = q;
        }
        start.set(p); // the segment starts here
        if (p < K-1 && sorted_priority[p].second == sorted_priority[p+1].second)
            start.reset(p + 1); // and, in this case, it is joined with the next segment
        // start.printB();
        // printPriority();
    }
};
