function rho = growth_factor(A,n)
U=upper_triangularisation(A,n);
max_A=0.0;
max_U=0.0;
for i=1:n
for j=1:n
abs_A=absolute_value(A(i,j));
abs_U=absolute_value(U(i,j));
if(abs_A>max_A)
max_A=abs_A;
end
if(abs_U>max_U)
max_U=abs_U;
end
end
end
rho=max_U/max_A;