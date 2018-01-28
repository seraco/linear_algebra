function columnvec = get_column(A,m,n)
% This function extracts a particular column from a given matrix
% A - The given matrix
% m - The column number of A that we wish to extract
% n - The size of the column vector 
for i=1:n
columnvec(i) = A(i,m);
end
