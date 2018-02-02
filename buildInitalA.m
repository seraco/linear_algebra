function A = buildInitalA(n)
% The function builds an initial matrix of the form
%       A=[
%            1  0  0  0  0  1
%           -1  1  0  0  0  1
%           -1 -1  1  0  0  1
%           -1 -1 -1  1  0  1
%           -1 -1 -1 -1  1  1
%           -1 -1 -1 -1 -1  1
%       ]
% n - the size of the matrix A

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