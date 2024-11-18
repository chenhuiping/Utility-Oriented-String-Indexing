/**
    USI: Useful String Indexing
    Copyright (C) 2024 Solon P. Pissis.

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

#ifndef __INCLUDE_DEF_H
#define __INCLUDE_DEF_H

#include <sdsl/bit_vectors.hpp>
#define ALLOC_SIZE 1048576
#define DEL '$'
#define DEL_STR "$"

#define DNA "ACGTN"
#define PROT "ARNDCQEGHILKMFPSTWYV"
#define IUPAC "ACGTUWSMKRYBDHVN"

using namespace sdsl;
using namespace std;

#ifdef _USE_64
typedef int64_t INT;
#endif

#ifdef _USE_32
typedef int32_t INT;
#endif

struct TSwitch {
    char* alphabet;
    char* input_filename;
    char* output_filename;
    unsigned int k;
    unsigned int K;
    unsigned int r;
    unsigned int total_length;
};

int decode_switches(int argc, char* argv[], struct TSwitch* sw);
void usage(void);

#endif