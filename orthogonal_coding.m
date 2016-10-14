function [kmer_coding_m_o,  seq_length ] = orthogonal_coding( filename, motif_length )
% Perform orthogonal encoding of the input sequences set.
% Input Arguments:
% filename: character vector specifying a file name or a path and file name. The referenced file should be a FASTA-formatted file. If you specify only a file name, that file must be on the MATLAB search path or in the MATLAB Current Folder.
% motif_length : the size of the sliding window.
% Output Arguments:
% kmer_coding_m_o: a single matrix which seconcatenates the coding matrices of all input sequence together.
% seq_length: a array whose i-th entry records the size of the coding
% matrix of the i-th input matrix. Therefore, size(kmer_coding_m_o,1) = sum(seq_length).

[kmer_coding, seq_length] = swindow(filename,motif_length);
kmer_coding_m_a = zeros(length(kmer_coding),motif_length);
for i = motif_length:-1:1
    kmer_coding_m_a(:,i) = mod(kmer_coding, 4) + (1 +(i-1)*4) ;
    kmer_coding = floor(kmer_coding/4);
end

n = size(kmer_coding,1);
kmer_coding_m_a = kmer_coding_m_a';

kmer_coding_m_a = kmer_coding_m_a + repmat([0:motif_length*4:motif_length*4*(n-1)],motif_length,1);

row_id = [1:size(kmer_coding_m_a,1)];
row_id = repmat(row_id,motif_length,1);
row_id = reshape(row_id,[],1);

kmer_coding_m_o = zeros(motif_length*4, size(kmer_coding_m_a,2));
kmer_coding_m_o(kmer_coding_m_a) = 1;

kmer_coding_m_o = kmer_coding_m_o';