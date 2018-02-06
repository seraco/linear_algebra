function B = Matrix1D(n)
% The function builds a matrix of the form
%       B = [
%            2 -1  0  0  0  0
%           -1  2 -1  0  0  0
%            0 -1  2 -1  0  0
%            0  0 -1  2 -1  0
%            0  0  0 -1  2 -1
%            0  0  0  0 -1  2
%       ]
% n - the size of the matrix B

iIndices = zeros(1,n+2*(n-1));
jIndices = zeros(1,n+2*(n-1));
values = zeros(1,n+2*(n-1));
k = 1;

for i=1:n
    for j=1:n
        if(i == j)
            iIndices(k) = i;
            jIndices(k) = j;
            values(k) = 2.0;
            k = k+1;
        elseif(i==j+1 || i==j-1)
            iIndices(k) = i;
            jIndices(k) = j;
            values(k) = -1.0;
            k = k+1;
        end
    end
end

B = sparse(iIndices,jIndices,values);
