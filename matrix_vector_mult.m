function c = matrix_vector_mult(A,u,n)
% This function computes matrix multiplication between two matrices
% A - The given matrix
% u - The given vector
% n - The size of the matrix and vector
for i=1:n
c(i) = inner_product(get_row(A,i,n),u,n);
end