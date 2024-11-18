/*
 *  This file is part of rk-lce.
 *  Copyright (c) by
 *  Nicola Prezza <nicolapr@gmail.com>
 *
 *   rk-lce is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   rk-lce is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details (<http://www.gnu.org/licenses/>).
 */

/*
 * rk_lce.hpp
 *
 *  Created on: Jun 4, 2016
 *      Author: Nicola Prezza
 *
 *  Encodes a text with prefix-inclusive Rabin-Karp hash function.
 *  This structure supports random access to the text and LCE queries
 *
 *  Space: n*ceil(log_2 sigma) bits + O(1) words
 *  Supports:
 *  	- access to the text in in O(1) time per block of 128/log_2 sigma characters
 *  	- LCP between any two text suffixes in O(log n) time
 *
 *  The class moreover can return a function to lexicographically-compare in O(log n) time any two
 *  text suffixes (useful for suffix-sorting in-place any subset of text positions)
 *
 *  Note: there is a probability of getting a wrong LCP result due to hash collisions.
 *  however, this probability is less than 2^-120 for real-case texts
 *
 */

#ifndef INTERNAL_RK_LCE_HPP_
#define INTERNAL_RK_LCE_HPP_

#include "includes.hpp"
#include "rk_lce_bin.hpp"

namespace rklce {

class rk_lce {
   public:
    static constexpr uint16_t w = 127;

    rk_lce(string filename) {
        char_to_uint = vector<uint8_t>(256);
        uint_to_char = vector<char>(256);

        n = 0;
        int sigma_temp = 0;
        {
            vector<bool> mapped(256, false);

            ifstream ifs(filename);

            char c;

            while (ifs.get(c)) {
                if (not mapped[uint8_t(c)]) {
                    mapped[uint8_t(c)] = true;

                    char_to_uint[uint8_t(c)] = sigma_temp;
                    uint_to_char[sigma_temp] = c;

                    sigma_temp++;
                }

                n++;
            }
        }

        assert(sigma_temp > 0);

        sigma = 1;
        log2_sigma = 0;

        while (sigma < sigma_temp) {
            sigma *= 2;
            log2_sigma++;
        }

        log2_sigma = log2_sigma == 0 ? 1 : log2_sigma;

        pad = w - (n * log2_sigma) % w;

        auto binary_text = vector<bool>(pad, false);

        {
            ifstream ifs(filename);

            char c;

            while (ifs.get(c)) {
                auto bc = char_to_uint[uint8_t(c)];

                for (int j = 0; j < log2_sigma; ++j) {
                    bool b = (bc >> (log2_sigma - j - 1)) & uint8_t(1);
                    binary_text.push_back(b);
                }
            }
        }

        bin_lce = rk_lce_bin(binary_text);
    }

    char operator[](uint64_t i) {
        auto ib = i * log2_sigma + pad;

        auto C = bin_lce(ib, log2_sigma);

        C = C >> (128 - log2_sigma);

        assert(C < uint_to_char.size());

        return uint_to_char[C];
    }

    uint64_t LCE(uint64_t i, uint64_t j) {
        auto ib = i * log2_sigma + pad;
        auto jb = j * log2_sigma + pad;

        return bin_lce.LCE(ib, jb) / log2_sigma;
    }

    uint64_t LCE_naive(uint64_t i, uint64_t j) {
        if (i == j) return n - i;

        uint64_t lce = 0;

        while (i + lce < n and j + lce < n and operator[](i + lce) == operator[](j + lce)) lce++;

        return lce;
    }

    std::function<bool(uint64_t, uint64_t)> lex_less_than() {
        return [&](uint64_t i, uint64_t j) {
            if (i == j) return false;

            uint64_t min = i < j ? j : i;
            auto lce = LCE(i, j);
            if (lce == n - 1 - min) return i != min;

            assert(i + lce < n - 1);
            assert(j + lce < n - 1);

            auto ic = operator[](i + lce);
            auto jc = operator[](j + lce);

            return ic < jc;
        };
    }

    uint64_t bit_size() { return bin_lce.bit_size() + sizeof(this) * 8 + char_to_uint.size() * 8 + uint_to_char.size() * 8; }

    uint64_t length() { return n; }
    uint64_t size() { return n; }

    uint16_t alphabet_size() { return sigma; }

   private:
    vector<uint8_t> char_to_uint;
    vector<char> uint_to_char;

    uint64_t n;

    uint64_t pad;

    uint16_t sigma;

    uint16_t log2_sigma;

    rk_lce_bin bin_lce;
};

}  // namespace rklce

#endif
