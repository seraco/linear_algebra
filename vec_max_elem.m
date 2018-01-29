function k = vec_max_elem(u,n)
% This function finds the index of the maximum element
% u - The input vector
% n - The size of the vector

k=1;
for i=1:n
if(u(i)>u(k))
k=i;
end
end