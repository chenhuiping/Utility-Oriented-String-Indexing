/**
    HEAVYKEEPERPLUS
    Branch from https://github.com/papergitkeeper/heavy-keeper-project
    Tong Yang, Haowei Zhang, Jinyang Li, Junzhi Gong, Steve Uhlig, Shigang Chen, Xiaoming Li:
    HeavyKeeper: An Accurate Algorithm for Finding Top-k Elephant Flows. IEEE/ACM Trans. Netw. 27(5): 1845-1858 (2019)
 
    Mods
    1) Karp-Rabin fingerprints for strings
    2) universal hash for the count-max sketch
    3) use space-efficient heapmindict.hpp (and bitvectornear.hpp) instead of ssummary.h
    4) extends the approach to (overlapping) substrings of the text
    5) rewritten the code

    Copyright (C) 2024 Roberto Grossi and Veronica Guerrini

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

#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <random>
#include <limits>

#include "krfp.h"
#include "unordered_dense.h"  
#include "minheapdict.hpp"

using namespace std;
using namespace karp_rabin_hashing;

constexpr int64_t MAX_NROWS = 2;           
constexpr int64_t MAX_NCOLS = 268435456;   
constexpr double HK_b = 1.08;             // base of the exponential probability to decrease a counter (original heavykepper is 1.08)
constexpr double HK_c = 1.08;

class HeavyKeeperPlus
{
public:

    HeavyKeeperPlus(int64_t _M2, int64_t _K) : heapDict(_K), M2(_M2), K(_K), gen(rd()), dis(0, std::numeric_limits<uint64_t>::max())
    {
        karp_rabin_hashing::init(); 
        if (_M2 > MAX_NCOLS)
            throw runtime_error("please increase MAX_NCOLS");
        for (int64_t i = 0; i < MAX_NROWS; i++){
            CMS_HK[i] = new node[_M2];
            for (int64_t j = 0; j < _M2; j++)
                CMS_HK[i][j].C = CMS_HK[i][j].FP = 0;
        }
    }
    

    void extractTopK(unsigned char * text, int64_t n, int64_t K, int64_t L, vector<pair<int64_t,uint32_t>> &topK) 
    {
		
		ankerl::unordered_dense::map<uint64_t,pair<int64_t,uint32_t>> fp_to_s; // maps fingerprints to their strings

        auto FP = new uint64_t[n];
        auto alphabet_size = 0;
        {
            ankerl::unordered_dense::set<unsigned char> distinct_chars;
            distinct_chars.insert(text[0]);
            FP[0] = karp_rabin_hashing::concat(0, text[0], 1);					// fingerprint of the first letter
            for(int64_t i = 1; i < n; ++i) {
                distinct_chars.insert(text[i]);
                FP[i] = karp_rabin_hashing::concat(FP[i-1], text[i], 1);	// all other fingerprints
            }
            alphabet_size = distinct_chars.size();
            distinct_chars.clear();
        }
			
        for (int64_t i = 0; i <= n-L; i++)
        {

            bool success = false;
            uint64_t itsfp; 
            for (uint32_t len = 1; len <= L; len++) {
                itsfp = KR_FP(FP, i, len);
                success = insert(itsfp);
                if (success){
                    fp_to_s[itsfp] = make_pair(i, len);
                }
            }
				
            if (!heapDict.exists(itsfp))
                continue;
			
            for (uint32_t len = L+1; i+len-1 < n; len++){
                auto old_itsfp = itsfp; 
                itsfp = KR_FP(FP, i, len);
                if (heapDict.exists(itsfp)){
                    success = insert(itsfp);
                    continue;
                } else {  
                    success = insert(itsfp);
                    if (success){
                        fp_to_s[itsfp] = make_pair(i, len);
                    } else {
                        break;
                    }
					if (heapDict.getSize() < K) continue;
                }
                uint64_t modval = std::numeric_limits<uint64_t>::max();
                double powval = pow(HK_c, len);
                if (modval > powval) modval = powval;
                if ( !(myrand64() % modval ) ) break; 
            }
        }

		int64_t maxlen = 0;
        for (int64_t i = 0; i < heapDict.getSize(); i++){
            auto x = fp_to_s[heapDict.get(i)];
            if (maxlen < x.second) maxlen = x.second;
			topK.push_back(x);
        }
		cout << "The longest top-K substring is of length = " << maxlen << endl;

        delete[] FP;

    }

private:
    // sketching table
    int64_t M2, K; // M2 <= MAX_NCOLS
    struct node
    {
        int64_t C;
        uint64_t FP;
    };
    node * CMS_HK[MAX_NROWS];

    MinHeapDictionary heapDict;

    // Random device to seed the generator
    std::random_device rd;
    // 64-bit Mersenne Twister random number generator
    std::mt19937_64 gen;
    // Uniform distribution for uint64_t
    std::uniform_int_distribution<uint64_t> dis;

    uint64_t myrand64(){
        return dis(gen);
    }

    uint64_t UniversalHash(const uint64_t x, int64_t j)
    { 
        if (j == 0)
        {
            constexpr uint64_t a1 = 0x65d200ce55b19ad7UL; 
            constexpr uint64_t a2 = 0x4f2162926e40c299UL; 
            return (static_cast<uint64_t>((x ^ (x >> 33U)) * a1) & 0xffffffff00000000UL) | (static_cast<uint64_t>((x ^ (x >> 33U)) * a2) >> 32U); // see https://en.wikipedia.org/wiki/Universal_hashing
        }
        else
        {
            constexpr uint64_t a3 = 0x68b665e6872bd1f5UL; 
            constexpr uint64_t a4 = 0xb6cfcf9d79b51db3UL; 
            return (static_cast<uint64_t>((x ^ (x >> 33U)) * a3) & 0xffffffff00000000UL) | (static_cast<uint64_t>((x ^ (x >> 33U)) * a4) >> 32U);
        }
    }

    uint64_t KR_FP(uint64_t * FP, int64_t i, size_t len)
    {
        return (i == 0)? FP[len - 1] : karp_rabin_hashing::subtract(FP[i + len - 1], FP[i - 1], len);
    }


    bool insert(uint64_t x)
    {
        bool x_is_in_topk = heapDict.exists(x);
        bool success = false;

        int64_t maxv = 0;
        for (int64_t j = 0; j < MAX_NROWS; j++) 
        {
            uint64_t Hsh = UniversalHash(x, j) % M2; 
            if (CMS_HK[j][Hsh].FP == x)             
            {
               CMS_HK[j][Hsh].C++;
               maxv = max(maxv, CMS_HK[j][Hsh].C);
            }
            else 
            {
                uint64_t modval = std::numeric_limits<uint64_t>::max();
                double powval = pow(HK_b, CMS_HK[j][Hsh].C);
                if (modval > powval) modval = powval;
                if ( !(myrand64() % modval ) )
                {
                    CMS_HK[j][Hsh].C--;
                    if (CMS_HK[j][Hsh].C <= 0)
                    {
                        CMS_HK[j][Hsh].FP = x;
                        CMS_HK[j][Hsh].C = 1;
                        maxv = max(maxv, static_cast<int64_t>(1));
                    }
                }
            }
        }
        if (!x_is_in_topk) 
        {
            if (maxv == 1 + heapDict.getMin() || heapDict.getSize() < K) 
            {
                heapDict.update(x, maxv);
                success = true;
            }
        }
        else
        {
            if (maxv > heapDict.getPriority(x)) 
            {
                if (maxv > 1 + heapDict.getPriority(x)){
                    throw runtime_error("large priority");
                }
                heapDict.update(x, maxv);
                success = true;
            }
        }
        return success;
    }
};
