function gv = lower_gauss_vector(A,m,n)
% The function constructs the Gauss vector of a given step of Gauss elimination method.
% A - the coefficient matrix
% m - the step number of Gauss elimination
% n - the size of the matrix A
a = n+1-m;
for i=n:-1:a %lines 6-8: initialisation of the Gauss vector
gv(i) = 0.0;
end
for j=a-1:-1:1 %computation of Gauss vector components according to the formula given in the lecture 3, page 10.
gv(j) = A(j,a)/A(a,a);
end
