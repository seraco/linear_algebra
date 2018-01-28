function rowvec = get_row(A,m,n)
% This function extracts a particular row from a given matrix
% A - The given matrix
% m - The row number of A that we wish to extract
% n - The size of the row vector 
for i=1:n
rowvec(i) = A(m,i);
end
