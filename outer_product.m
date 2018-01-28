function C = outer_product(a,b,n)
% This function computes the outer product between two vectors and results a matrix
% a - the first vector
% b - the second vector
% C - the resulting matrix
for i=1:n
for j=1:n
C(i,j) = (a(i)*b(j)); % formula for outer product according to lecture 1, page 9.
end
end
