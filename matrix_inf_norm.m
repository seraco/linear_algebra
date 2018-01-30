function norm = matrix_inf_norm(A,n)
% The function computes infinite norm of matrix A
% A - the input matrix
% n - the size of the matrix A

norm=0.0;
for i=1:n
row_sum=0.0;
for j=1:n
row_sum=row_sum+absolute_value(A(i,j));
end
if(row_sum>norm)
norm=row_sum;
end
end