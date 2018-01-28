function x = forward_substitution(L,b,n)
% The function solves a system of equations with forward substitution.
% L - the coefficient lower matrix
% b - the RHS of the system of equations Ux=b
% n - the size of the matrix U
for i=1:n
x(i)=0.0;
end
for i=1:n
for j=1:i-1
x(i)=x(i)+L(i,j)*x(j);
end
x(i)=(b(i)-x(i))/L(i,i);
end
