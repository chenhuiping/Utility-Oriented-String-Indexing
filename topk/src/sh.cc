/**
	SH: Substring HeavyKeeper
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
    along with this program.  If not, see <http://www.gnu.org/licenses/>.#include <iostream>
**/
#include <string>
#include <chrono>
#include <ctime>
#include <vector>
#include <cmath>
#include <omp.h>
#include <sys/time.h>
#include <tuple>
#include "defs.h"
#include "krfp.h"
#include "utils.h"
#include "unordered_dense.h"

#include "heavykeeperplus.hpp"

using namespace std;

void output_top_k(unsigned char * sequence, vector<pair<int64_t,uint32_t>> &topK, string filename) {
    ofstream of(filename, ios::out);
    if (!of) {
        cout << "Error opening file: " << filename << endl;
        exit(1);
    }

    for (auto &u : topK) {
		for(uint32_t j=0; j<u.second;j++)
			of << sequence[u.first+j]; 
		of << endl;
    }
    of.close();
}

int main(int argc, char* argv[])
{
	if(argc<3)
	{	
		cout << "Usage: " << argv[0] << " text-file " << " K " << endl;
		return -1;
	}
	
	INT K = atoll(argv[2]);
	int64_t L = 1;
	
	ifstream seq(argv[1], ios::in | ios::binary);
    unsigned char * input_seq_char = NULL;
	unsigned char * sequence;
	INT n = 0;
    char c;
    while (seq.get(c))
    {
        if(n == 0 || n % ALLOC_SIZE)	input_seq_char = ( unsigned char * ) realloc ( input_seq_char,   ( n + ALLOC_SIZE ) * sizeof ( unsigned char ) );
        input_seq_char[n] = c;
        n++;
    }
    seq.close();
	
	if( n == 0 || n % ALLOC_SIZE)    input_seq_char = ( unsigned char * ) realloc ( input_seq_char,   ( n + ALLOC_SIZE ) * sizeof ( unsigned char ) );
	if( c == '\n')	input_seq_char[--n] = 255;	
	else		input_seq_char[n] = 255;	
	//Temporary
	sequence = input_seq_char;
    sequence[++n]='\0';
	
	cout << "Text is of length n = " << n - 1 << "." << endl;
	cout << "K = " << K << endl;

	std::chrono::steady_clock::time_point  stream_begin = std::chrono::steady_clock::now();

	vector<pair<int64_t,uint32_t>> topK; 
	auto cms_hk_ncols =(256*K<268435456?256*K:268435456); // temporary: min between 256*K and 2^28
	HeavyKeeperPlus HeavySketch(cms_hk_ncols, K);
	
	HeavySketch.extractTopK(sequence, n, K, L, topK);

	std::chrono::steady_clock::time_point  stream_end = std::chrono::steady_clock::now();
	cout<<"Top-K computation time: "<< std::chrono::duration_cast<std::chrono::milliseconds>(stream_end - stream_begin).count() << "[ms]." << std::endl;	
	output_top_k(sequence,topK, string(argv[1]) + "_stream_freq_K" + to_string(K));
	free (sequence);
	return 0;
}
