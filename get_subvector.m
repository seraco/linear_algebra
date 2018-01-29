function vec = get_subvector(u,k,l)
% This function extracts the elements of a vector from k to l indices
% u - The given vector
% k - The initial index
% l - The final index
for i=k:l
vec(i+1-k) = u(i);
end