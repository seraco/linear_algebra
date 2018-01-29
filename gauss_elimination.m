function x = gauss_elimination(A,b,n)
% The function solves a system of equations with Gauss elimination method (GEM).
% A - the coefficient matrix
% b - RHS of equation Ax=b
% n - the size of the matrix A

[~,L,U,c] = partial_pivoting(A,b,n);
y=forward_substitution(L,c,n);
x=backward_substitution(U,y,n);