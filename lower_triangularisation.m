function L = lower_triangularisation(A,n)
% The function computes the lower triangularisation of a given square matrix A using the Gauss elimination method (GEM).
% A - the coefficient matrix
% n - the size of the matrix A
for i=1:n %lines 5-15: initialisation of the identity (I) and lower-triangular (L)
for j=1:n
if (i==j) 
I(i,j)=1.0;
L(i,j) = A(i,j);
else
I(i,j) = 0.0;
L(i,j) = A(i,j);
end
end
end

for k=1:n-1 % the main loop of the GEM starts here
gv = lower_gauss_vector(L,k,n); % the first step in GEM is to construct Gauss vector
evec = lower_e_vector(k,n); %second, I construct my e vector for the corresponding stage of the GEM
LE = outer_product(gv,evec,n); %Here, I take outer product between the Gauss and e vectors for the given stage
J = matrix_subtraction(I,LE,n); %Here, I construct my Gauss transformation matrix by subtracting the previously computed outer product from the identity matrix
L = matrix_matrix_mult(J,L,n); %Here, I compute the updated upper-triangular matrix through matrix multiplication with the Gauss transformation matrix.
end % the main loop of the GEM ends here