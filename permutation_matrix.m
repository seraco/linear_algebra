function P = permutation_matrix(i1,i2,n)
% The function returns the permutation matrix of rows i1-i2.
% i1 - 1st row index to be swapped
% i2 - 2nd row index to be swapped
% n - the size of the matrix to be returned

for i=1:n
for j=1:n
if (i==j) 
I(i,j)=1.0;
else
I(i,j) = 0.0;
end
end
end

e1=upper_e_vector(i1,n);
e2=upper_e_vector(i2,n);
e2me1=vector_subtraction(e2,e1,n);
e2me1T=transpose_matrix(e2me1,1,n);
out_prd=outer_product(e2me1T,e2me1,n);

for i=1:n
for j=1:n
P(i,j)=I(i,j)-out_prd(i,j);
end
end
