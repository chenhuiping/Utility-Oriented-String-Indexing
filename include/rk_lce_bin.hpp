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
 * rk_lce_bin.hpp
 *
 *  Created on: Jun 4, 2016
 *      Author: Nicola Prezza
 *
 *  Encodes a binary text with suffix-inclusive Rabin-Karp hash function. This structure supports access to the text and
 *  fast computation of the RK function of any text suffix
 *
 *  Space: n bits (n = size of original text)
 *  Supports:
 *  	- access to the text in in O(1) time per character
 *  	- LCE between any two text suffixes in O(log n) time
 *
 */

#ifndef INTERNAL_RK_LCE_BIN_HPP_
#define INTERNAL_RK_LCE_BIN_HPP_

#include "bitv.hpp"
#include "includes.hpp"

namespace rklce {

class rk_lce_bin {
   public:
    static constexpr uint64_t w = 127;

    static constexpr uint128 q = (uint128(1) << w) - 1;

    rk_lce_bin() {}

    rk_lce_bin(vector<bool> &input_bitvector) {
        n = input_bitvector.size();

        assert(n % w == 0);

        auto n_bl = n / w;

        assert(n_bl > 0);

        vector<uint128> P_vec;

        {
            auto B_vec = vector<uint128>(n_bl, 0);

            uint64_t i = 0;
            for (auto b : input_bitvector) {
                B_vec[i / w] |= (uint128(b) << (w - (i % w + 1)));
                i++;
            }

            assert(B_vec[0] != q);

            uint64_t number_full_blocks = 0;
            for (auto bl : B_vec) number_full_blocks += bl == q;

            {
                vector<bool> Q_vec;

                i = 0;
                for (auto bl : B_vec) Q_vec.push_back(bl == q);

                Q1 = bitv(Q_vec);

                assert(Q1.rank(Q1.size()) == number_full_blocks);
            }

            for (uint64_t i = 1; i < B_vec.size(); ++i) {
                B_vec[i] = (B_vec[i] + mul_pow2<w>(B_vec[i - 1], w)) % q;
            }

            P_vec = vector<uint128>(Q1.rank(Q1.size(), 0));

            uint64_t i_P_vec = 0;
            for (uint64_t j = 0; j < Q1.size(); ++j) {
                if (not Q1[j]) {
                    assert(i_P_vec < P_vec.size());
                    P_vec[i_P_vec++] = B_vec[j];
                }
            }
        }

        P = packed_vector_127(P_vec);
    }

    uint64_t bit_size() { return P.bit_size() + Q1.bit_size() + sizeof(this) * 8; }

    bool operator[](uint64_t i) {
        assert(i < n);

        auto ib = i / w;

        return (B(ib) >> (w - (i % w + 1))) & uint128(1);
    }

    uint128 operator()(uint64_t i, uint64_t len = 128) {
        assert(len <= 128);

        auto ib = i / w;

        uint128 L = B(ib) << (i % w + 1);
        uint128 R = ib == Q1.size() - 1 ? 0 : B(ib + 1) >> (w - (i % w) - 1);

        auto block = L | R;

        uint128 MASK = len == 128 ? ~0 : ((uint128(1) << len) - 1) << (128 - len);

        return block & MASK;
    }

    uint64_t LCE(uint64_t i, uint64_t j) {
        assert(i < n);
        assert(j < n);

        if (i == j) return n - i;

        if (i == n or j == n) return 0;

        auto i_len = n - i;
        auto j_len = n - j;

        auto i_block = operator()(i);
        auto j_block = operator()(j);

        if (i_block != j_block or i_len <= 128 or j_len <= 128) {
            auto lce = clz_u128(i_block ^ j_block);

            auto min = std::min(i_len, j_len);

            lce = lce > min ? min : lce;

            return lce;
        }

        auto lce = LCE_binary(i, j);
        assert(lce == LCE_naive(i, j));

        return lce;
    }

    bool equals(uint64_t i, uint64_t j, uint64_t l, uint128 i_fp = q, uint128 j_fp = q) {
        assert(i + l - 1 < n);
        assert(j + l - 1 < n);

        uint128 fpi = RK(i, i + l - 1, i_fp);
        uint128 fpj = RK(j, j + l - 1, j_fp);

        return fpi == fpj;
    }

    uint64_t LCE_naive(uint64_t i, uint64_t j) {
        if (i == j) return n - i;

        uint64_t lce = 0;

        while (i + lce < n and j + lce < n and operator[](i + lce) == operator[](j + lce)) lce++;

        return lce;
    }

    uint64_t number_of_blocks() { return Q1.size(); }

    uint64_t block_size() { return w; }

    uint64_t length() { return n; }
    uint64_t size() { return n; }

   private:
    uint128 RK(uint64_t i) {
        auto j = i / w;

        uint128 L = j == 0 ? 0 : mul_pow2<w>(P1(j - 1), i - j * w + 1);
        uint128 R = B(j) >> (w - (i % w + 1));

        return (L + R) % q;
    }

    uint128 RK(uint64_t i, uint64_t j, uint128 rki = q) {
        assert(j >= i);

        rki = rki != q ? rki : (i == 0 ? 0 : RK(i - 1));

        rki = mul_pow2<w>(rki, j - i + 1);

        uint128 rkj = RK(j);

        return sub<w>(rkj, rki);
    }

    uint128 P1(uint64_t i) {
        if (Q1.rank(Q1.size()) == 0) return P[i];

        auto t = Q1.predecessor_0(i);

        assert(Q1.rank(t, 0) < P.size());
        assert(t <= i);

        return mul_pow2<w>(P[Q1.rank(t, 0)], w * (i - t));
    }

    uint128 B(uint64_t i) {
        assert(i < Q1.size());

        assert(not Q1[0]);

        if (i == 0) return P[0];

        if (Q1[i]) return q;

        return sub<w>(P1(i), mul_pow2<w>(P1(i - 1), w));
    }

    uint64_t LCE_binary(uint64_t i, uint64_t j) {
        assert(i != j);

        auto i_len = n - i;
        auto j_len = n - j;

        auto sc = suffix_comparator(this, i, j);

        uint64_t k = 1;

        while (k < sc.size() && not sc[k]) k *= 2;
        if (k >= sc.size()) k = sc.size();

        auto lce = uint64_t(std::upper_bound(sc.begin(), sc.begin() + k, false)) - 1;

        return lce;
    }

    class suffix_comparator {
        class sc_iterator : public std::iterator<random_access_iterator_tag, bool> {
            friend class suffix_comparator;

            suffix_comparator *_sci = nullptr;
            uint64_t _index = 0;

            sc_iterator(suffix_comparator *v, uint64_t index) : _sci(v), _index(index) {}

           public:
            sc_iterator() = default;
            sc_iterator(sc_iterator const &) = default;

            operator uint64_t() { return _index; }

            sc_iterator &operator=(sc_iterator const &) = default;

            bool operator*() const { return (*_sci)[_index]; }

            sc_iterator &operator++() {
                ++_index;
                return *this;
            }

            bool operator==(sc_iterator it) const { return _index == it._index; }

            bool operator!=(sc_iterator it) const { return _index != it._index; }

            sc_iterator operator++(int) {
                sc_iterator it(*this);
                ++_index;
                return it;
            }

            sc_iterator &operator--() {
                --_index;
                return *this;
            }

            sc_iterator operator--(int) {
                sc_iterator it(*this);
                --_index;
                return it;
            }

            sc_iterator &operator+=(uint64_t n) {
                _index += n;
                return *this;
            }

            sc_iterator operator+(uint64_t n) const {
                sc_iterator it(*this);
                it += n;
                return it;
            }

            friend sc_iterator operator+(uint64_t n, sc_iterator it) { return it + n; }

            sc_iterator &operator-=(uint64_t n) {
                _index -= n;
                return *this;
            }

            sc_iterator operator-(uint64_t n) const {
                sc_iterator it(*this);
                it -= n;
                return it;
            }

            friend sc_iterator operator-(uint64_t n, sc_iterator it) { return it - n; }

            uint64_t operator-(sc_iterator it) { return uint64_t(_index) - uint64_t(it._index); }

            bool operator[](uint64_t i) const { return (*_sci)[_index + i]; }

            bool operator<(sc_iterator it) const { return _index < it._index; }

            bool operator<=(sc_iterator it) const { return _index <= it._index; }

            bool operator>(sc_iterator it) const { return _index > it._index; }

            bool operator>=(sc_iterator it) const { return _index >= it._index; }
        };

       public:
        suffix_comparator(rk_lce_bin *T, uint64_t i, uint64_t j) {
            assert(i != j);
            assert(i < T->size());
            assert(j < T->size());

            if (j < i) {
                auto t = i;
                i = j;
                j = t;
            }

            this->T = T;
            this->i = i;
            this->j = j;

            i_fp = i == 0 ? 0 : T->RK(i - 1);
            j_fp = j == 0 ? 0 : T->RK(j - 1);

            n = (T->size() - j) + 2;
        }

        uint64_t size() { return n; }

        bool operator[](uint64_t t) {
            assert(t < n);

            if (t == n - 1) return true;
            if (t == 0) return false;

            return not T->equals(i, j, t, i_fp, j_fp);
        }

        sc_iterator begin() { return sc_iterator(this, 0); }
        sc_iterator end() { return sc_iterator(this, n); }

       private:
        rk_lce_bin *T;

        uint64_t n;

        uint64_t i;
        uint64_t j;

        uint128 i_fp;
        uint128 j_fp;
    };

    packed_vector_127 P;

    bitv Q1;

    uint64_t n = 0;
};

}  // namespace rklce

#endif
