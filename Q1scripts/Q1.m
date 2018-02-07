%% Script to solve the first question.
%   The question is about Gaussian elimination with partial pivoting PA=LU.
%   This script computes P,L,U matrices and solves a system of equation
%   using this information. Then, it also studies the influence of the
%   matrix size with growth factor. Finally, the variations of two
%   parameters with respect to the matrix size is also assessed.

% Version:  0.0.1
% Author:   S.Ramon

oldpath = path;
path('../LinearAlgebra',oldpath);

%% iii
A = buildInitalA(5);
b = [2 1 0 -1 -3];

sys = CLinearSystem(A,b);
[P,L,U] = sys.luWithPartialPivoting();
x = sys.directSolution('lu');

fprintf('iii part solved. The answers can be checked.\n');
pause;

%% iv
rho = zeros(3,1);
k = 0;

for i=[8,16,32]
    k = k+1;
    A = buildInitalA(i);
    b = ones(i,1);
    sys = CLinearSystem(A,b);
    rho(k) = sys.growthFactor();
end

fprintf('iv part solved. The answers can be checked.\n');
pause;

%% v
parameter1 = zeros(3,1);
parameter2 = zeros(3,1);
k = 0;

for i=[8,16,32]
    k = k+1;
    A = buildInitalA(i);
    b = ones(i,1);
    sys = CLinearSystem(A,b);
    [P,L,U] = sys.luWithPartialPivoting();
    L = L.absoluteValue();
    U = U.absoluteValue();
    LU = L*U;
    parameter1(k) = LU.infiniteNorm();
    A = CMatrix(A);
    parameter2(k) = parameter1(k)/A.infiniteNorm();
end

fprintf('v part solved. The answers can be checked.\n');
pause;