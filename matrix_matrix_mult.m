function C = matrix_matrix_mult(A,B,n)
% This function computes matrix multiplication between two matrices
% A - The first given matrix
% B - The second given matrix
% n - The size of the matrices
for i=1:n
for j=1:n
C(i,j) = inner_product(get_row(A,i,n),get_column(B,j,n),n); %Matrix-matrix multiplication is nothing but dot products. I get a row vector from A and a column vector from B and then I compute dot product between them.
end
end
