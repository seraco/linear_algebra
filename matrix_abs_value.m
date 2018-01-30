function B = matrix_abs_value(A,n)
% The function converts all entries of a matrix to absolute value
% A - the input matrix
% n - the size of the matrix A

for i=1:n
for j=1:n
B(i,j)=absolute_value(A(i,j));
end
end