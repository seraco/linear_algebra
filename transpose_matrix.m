function B = transpose_matrix(A,n,m)
% The function returns the transpose of matrix A.
% A - input matrix
% n - the size of the matrix's rows
% m - the size of the matrix's columns

for i=1:n
for j=1:m
B(j,i)=A(i,j);
end
end