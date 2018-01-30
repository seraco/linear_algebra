function A = build_inital_A(n)
for i=1:n
for j=1:n
if(j==n)
A(i,j)=1.0;
elseif(i==j)
A(i,j)=1.0;
elseif(i>j)
A(i,j)=-1.0;
else
A(i,j)=0.0;
end
end
end