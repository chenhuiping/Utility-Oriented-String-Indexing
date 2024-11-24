/**
    BITVECTORNEAR: see below
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

#include <vector>

using namespace std;

/*
    Simple bitvector where supported operations are changing a bit value, and finding the nearest 1
    Optimized to use along with heavykeeper in place of ssummary.h
*/
class BitVectorNear
{
    public:
        BitVectorNear(int64_t _K){
            K = _K;
            pacK = _K / 8 + 1;
            B.reserve(pacK + 1); // packed bits!
            for (auto i = 0; i <= pacK; i++)
                B[i] = 0;
            B[_K >> 3] =  mask[_K & 0x07]; // set Kth bit to 1 as a guard
            n1 = 1;
            nCalls = 0;
            nOps = 0;
        }

        uint8_t get(int64_t i){
            if (i < 0 || i > K)
                throw runtime_error("get index out of range");
            return B[i >> 3] & mask[i & 0x07] ? 1 : 0;
        }

        void set(int64_t i){
            if (i < 0 || i >= K)
                throw runtime_error("set index out of range");
            if (get(i)) return;
            B[i >> 3] |=  mask[i & 0x07]; // set bit to 1
            n1++;
        }
        void reset(int64_t i){
            if (i < 0 || i >= K)
                throw runtime_error("reset index out of range");
            if (!get(i)) return;
            B[i >> 3] &=  ~(mask[i & 0x07]); // set bit to 0
            n1--;
        }

        int64_t nearestOneRight(int64_t i){
            if (i < 0 || i >= K)
                throw runtime_error("rank index out of range");
            // nearest 1 always exists, it is K in the worst case
            nCalls++;
            nOps++;
            auto pos = i+1;
            while ( (pos >> 3) == (i >> 3) ){
                nOps++;
                if (get(pos)) return pos;
                pos++;
            }
            auto pos3 = pos >> 3;
            while ( !B[pos3] ) {pos3++; nOps++;} // 8 bit at a time; todo: maybe speed up this...
            pos = pos3 << 3;
            while ( !get(pos) ) {pos++; nOps++;}
            return pos;
        }

        void printB(){
            cout << K << ": ";
            for (auto i = 0; i <= K; i++)
                cout << (int)(get(i));
            cout << endl << flush;
        }

        void printStat(){
            cout << (nOps + 0.0) / nCalls << " average ops per NearestOneRight call" << endl;
        }
    private:
        uint8_t mask[8] = {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80}; // to single out a bit from a byte
        vector<uint8_t> B; // packed bitvector of pacK bytes for K+1 bits, where the last bit is 1
        int64_t n1 = 1;  // number of 1s >= 1
        int64_t K = 0;
        int64_t pacK = 0;
        // stats
        int64_t nCalls = 0;
        int64_t nOps = 0;
};
