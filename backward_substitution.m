function x = backward_substitution(U,b,n)
% The function solves a system of equations with backward substitution.
% U - the coefficient upper matrix
% b - the RHS of the system of equations Ux=b
% n - the size of the matrix U
for i=1:n
x(i)=0.0;
end
for i=n:-1:1
for j=n:-1:i+1
x(i)=x(i)+U(i,j)*x(j);
end
x(i)=(b(i)-x(i))/U(i,i);
end
