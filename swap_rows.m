function C = swap_rows(A,k,l,o,p,n)
% The function swaps the rows k-l of matrix A from column o to p
% A - the input matrix
% k - 1st index of row to be swapped
% l - 2nd index of row to be swapped
% o - initial index of columns
% p - final index of columns

for i=1:n
for j=1:n
C(i,j)=A(i,j);
end
end

for j=o:p
C(k,j)=A(l,j);
C(l,j)=A(k,j);
end