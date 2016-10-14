//FileName: swindow.cpp
//Brief: Pile up a matrix that column is window width
//A-00 T-11 C-01 G-10
//This is a mex file for matlab
//Input: 1. the fasta file name. 2. the width of window
//Output: 1. the vector of slidwindow. 2. the vector of every sequence's size

#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include<ppl.h>
#include<mex.h>

using namespace std;
using namespace concurrency;

int ReadFasta(string filename, vector<string>& sequences, int wwidth, vector<int>& ibsequence, vector<int>& evseqsize) {
    ifstream infile;
    infile.open(filename);											//open the file
    //string numsequence, sequence;
    string sequence;												//every sequence in the input file are in the variable
    int indexsequence = 0;											//for recording the every sequence's begin index 
    if(infile) {
        //while(infile >> numsequence >> sequence) {
        while(getline(infile, sequence)) {							//read every line in the file
            if(sequence.size() != 0 && sequence[0] != '>') {		//skip the blank line or the '>' line
            //if(sequence[0] != '>') {
                sequences.push_back(sequence);						//put the sequence in the vector
                ibsequence.push_back(indexsequence);				//put the begin index in the vector
                int num = (sequence.size() - wwidth + 1) * 2;		//calculate every sequence will generate how many slidwindows
                indexsequence += num;								//next begin index is indexnow + num
                evseqsize.push_back(num);
            }
        }
    }
    return indexsequence;											//return the size of slidwindow
}

void PushBackExactStrand(unsigned long long& exactstrand, char c) {
    exactstrand <<= 2;
    switch(c) {
        case 'A': break;
        case 'T': exactstrand |= 3ULL; break;
        case 'C': exactstrand |= 1ULL; break;
        case 'G': exactstrand |= 2ULL; break;
    }
}

void PushFrontComplementaryStrand(unsigned long long& complementarystrand, char c, int wwidth) {
    complementarystrand >>= 2;
    switch(c) {
        case 'A': complementarystrand |= 3ULL << (wwidth-1)*2; break;
        case 'T': break;
        case 'C': complementarystrand |= 2ULL << (wwidth-1)*2; break;
        case 'G': complementarystrand |= 1ULL << (wwidth-1)*2; break;
    }
}

void GetSlidWindowMatrix(vector<string>& sequences, vector<unsigned long long>& swmatrix, int wwidth, vector<int>& ibsequence, unsigned long long amendflag) {
	int szsequences = sequences.size();
    parallel_for(0, szsequences, [&] (int k) {
        string sequence = sequences[k];
        unsigned long long exactstrand = 0;
		unsigned long long complementarystrand = 0;
        int num = ibsequence[k];
        int i;
        for(i = 0; i < wwidth; i++) {
            PushBackExactStrand(exactstrand, sequence[i]);
            PushFrontComplementaryStrand(complementarystrand, sequence[i], wwidth);
        }
        exactstrand &= amendflag;
        complementarystrand &= amendflag;
        swmatrix[num++] = exactstrand;
        swmatrix[num++] = complementarystrand;

        while(i < sequence.size()) {
            PushBackExactStrand(exactstrand, sequence[i]);
            PushFrontComplementaryStrand(complementarystrand, sequence[i], wwidth);
            exactstrand &= amendflag;
            complementarystrand &= amendflag;
            swmatrix[num++] = exactstrand;
            swmatrix[num++] = complementarystrand;

            i++;
        }
    });
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	char* infilename;							//input file's name;
    infilename = mxArrayToString(prhs[0]);
    int wwidth = mxGetScalar(prhs[1]);			//slidwindow's width
	vector<string> sequences;					//sequences in input file
    vector<int> evseqsize;						//every sequence's size (for the slidwindow)
    vector<int> ibsequence;						//every sequence's begin index (for the slidwindow)
	
    unsigned long long amendflag = 0;			//because we use the two bit represent one DNA character
    amendflag = ~amendflag;						//in order to just keep the last 2*wwidth bit there is a amendflag.
    amendflag = amendflag << wwidth*2;			//so that every number need to operate & with amendflag
    amendflag = ~amendflag;
	
	int szswmatrix = ReadFasta(infilename, sequences, wwidth, ibsequence, evseqsize);	//the size of slidwindow for all the sequences
	vector<unsigned long long> swmatrix(szswmatrix, 0);									//the output vector of slidwindow
	GetSlidWindowMatrix(sequences, swmatrix, wwidth, ibsequence, amendflag);

	double *op1, *op2;
	int mswmatrix = swmatrix.size();
	int mevseqsize = evseqsize.size();
	plhs[0] = mxCreateDoubleMatrix(mswmatrix, 1, mxREAL);
	op1 = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(mevseqsize, 1, mxREAL);
	op2 = mxGetPr(plhs[1]);
	for(int i = 0; i < mswmatrix; i++) op1[i] = swmatrix[i];
	for(int i = 0; i < mevseqsize; i++) op2[i] = evseqsize[i];
}

