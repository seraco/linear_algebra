function c = inner_product(a,b,n)
% This function computes the dot product between two vectors.
% a - The first vector
% b - The second vector
% n - The size of the vectors a and b.
c = 0.0;
for i=1:n
c = c + (a(i)*b(i));
end
