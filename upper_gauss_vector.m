function gv = upper_gauss_vector(A,m,n)
% The function constructs the Gauss vector of a given step of Gauss elimination method.
% A - the coefficient matrix
% m - the step number of Gauss elimination
% n - the size of the matrix A
for i=1:m %lines 6-8: initialisation of the Gauss vector
gv(i) = 0.0;
end
for j=m+1:n %computation of Gauss vector components according to the formula given in the lecture 3, page 10.
gv(j) = A(j,m)/A(m,m);
end
