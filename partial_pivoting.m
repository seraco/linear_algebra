function [P,L,U] = partial_pivoting(A,n)
% The function computes the upper triangularisation of a given square matrix A using the Gauss elimination method (GEM).
% A - the coefficient matrix
% b - RHS of equation Ax=b
% n - the size of the matrix A

for i=1:n
for j=1:n
if (i==j) 
L(i,j)=1.0;
I(i,j)=1.0;
P(i,j)=1.0;
else
L(i,j)=0.0;
I(i,j)=0.0;
P(i,j)=0.0;
end
U(i,j)=A(i,j);
end
end

for k=1:n-1
U_col=get_column(U,k,n);
U_col=get_subvector(U_col,k,n);
max_i=vec_max_elem(U_col,n+1-k)+k-1;
perm=permutation_matrix(max_i,k,n);
P=matrix_matrix_mult(perm,P,n);
U=swap_rows(U,max_i,k,k,n,n);
L=swap_rows(L,max_i,k,1,k-1,n);
gv=gauss_vector(U,k,n);
evec=e_vector(k,n);
LE=outer_product(gv,evec,n);
J=matrix_subtraction(I,LE,n);
U=matrix_matrix_mult(J,U,n);
L=matrix_addition(L,LE,n);
end