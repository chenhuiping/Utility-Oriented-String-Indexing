// Compile with g++ -std=c++20 -O3 utils.cpp kr_greg.cpp -o kr_greg
// run with ./kr_greg

// Added a subtract function

/**
 * @file    karp_rabin_hashing.cpp
 * @section LICENCE
 *
 * This file is part of Lazy-AVLG v0.1.0
 * See: https://github.com/dominikkempa/lz77-to-slp
 *
 * Copyright (C) 2021
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#include "krfp.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "utils.h"

using namespace std;

namespace karp_rabin_hashing {

std::uint64_t hash_variable;
std::uint64_t mersenne_prime_exponent;

inline std::uint64_t mul_mod_mersenne(const std::uint64_t a, const std::uint64_t b, const std::uint64_t k) {
    const std::uint64_t p = ((std::uint64_t)1 << k) - 1;
    __extension__ const unsigned __int128 ab = (unsigned __int128)a * (unsigned __int128)b;
    std::uint64_t lo = (std::uint64_t)ab;
    const std::uint64_t hi = (ab >> (uint64_t)64);
    lo = (lo & p) + ((lo >> k) + (hi << ((uint64_t)64 - k)));
    lo = (lo & p) + (lo >> k);
    return lo == p ? (uint64_t)0 : lo;
}

inline std::uint64_t mod_mersenne(std::uint64_t a, const std::uint64_t k) {
    std::uint64_t p = ((std::uint64_t)1 << k) - (uint64_t)1;
    if (k < (uint64_t)32) {
        const std::uint64_t threshold = ((std::uint64_t)1 << (k << (uint64_t)1));
        if (a <= threshold) {
            a = (a & p) + (a >> k);
            a = (a & p) + (a >> k);
            return a == p ? (uint64_t)0 : a;
        } else
            return a % p;
    } else {
        a = (a & p) + (a >> k);
        a = (a & p) + (a >> k);
        return a == p ? (uint64_t)0 : a;
    }
}

std::uint64_t rand_mod_mersenne(const std::uint64_t k) {
    const std::uint64_t p = ((std::uint64_t)1 << k) - (uint64_t)1;
    return utils::random_int<std::uint64_t>((std::uint64_t)0, (std::uint64_t(p - (uint64_t)1)));
}

std::uint64_t pow_mod_mersenne(const std::uint64_t a, std::uint64_t n, const std::uint64_t k) {
    std::uint64_t pow = mod_mersenne(a, k);
    std::uint64_t ret = mod_mersenne((uint64_t)1, k);
    while (n > (uint64_t)0) {
        if (n & (uint64_t)1) ret = mul_mod_mersenne(ret, pow, k);
        pow = mul_mod_mersenne(pow, pow, k);
        n >>= (uint64_t)1;
    }
    return ret;
}

std::uint64_t concat(const std::uint64_t left_hash, const std::uint64_t right_hash, const std::uint64_t right_len) {
    const std::uint64_t pow = pow_mod_mersenne(hash_variable, right_len, mersenne_prime_exponent);
    const std::uint64_t tmp = mul_mod_mersenne(left_hash, pow, mersenne_prime_exponent);
    const std::uint64_t ret = mod_mersenne(tmp + right_hash, mersenne_prime_exponent);
    return ret;
}

std::uint64_t subtract(const std::uint64_t long_hash, const std::uint64_t short_hash, const std::uint64_t right_len) {
    const std::uint64_t pow = pow_mod_mersenne(hash_variable, right_len, mersenne_prime_exponent);
    const std::uint64_t tmp = mul_mod_mersenne(short_hash, pow, mersenne_prime_exponent);
    const std::uint64_t p = ((std::uint64_t)1 << mersenne_prime_exponent) - 1;
    return (long_hash >= tmp) ? (long_hash - tmp) : ((long_hash + p) - tmp);
}

std::uint64_t subtract_fast(const std::uint64_t long_hash, const std::uint64_t short_hash, const std::uint64_t power) {
    const std::uint64_t pow = power;
    const std::uint64_t tmp = mul_mod_mersenne(short_hash, pow, mersenne_prime_exponent);
    const std::uint64_t p = ((std::uint64_t)1 << mersenne_prime_exponent) - 1;
    return (long_hash >= tmp) ? (long_hash - tmp) : ((long_hash + p) - tmp);
}

std::uint64_t init() {
    mersenne_prime_exponent = 61;
    hash_variable = rand_mod_mersenne(mersenne_prime_exponent);
    return hash_variable;
}
}  // namespace karp_rabin_hashing
