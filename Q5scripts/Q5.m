%% Script to solve Question 5.
%   The question is about a steady heat conduction problem in a L shaped
%   region. This script computes the matrix of coefficients in the problem,
%   its eigenvalues and bandwidth. Finally, it solves the problem and plots
%   the solution.

% Version:  0.0.1
% Author:   S.Ramon

addpath('../LinearAlgebra','../PhysicalProblems');

%% 1)
problem = CSteadyBidimensionalHeat(10,10,1,1,1);
A = -problem.coefficientsMatrix();
eigenvalues = eig(A);
[lower,upper] = bandwidth(A);
bwdth = max(lower,upper);
sizeOfA = size(A);

%% 2)
problem = CSteadyBidimensionalHeat(5,5,1,1,1);
A = -problem.coefficientsMatrix();
spy(A);

%% 3)
Nx = 3;
Ny = 3;
a = 1;
q = 1;
k = 1;
problem = CSteadyBidimensionalHeat(Nx,Ny,a,q,k);
solution = problem.solve();
leftSolution = solution.data(1:(Nx-1)*(2*Ny-1));
rightSolution = solution.data((Nx-1)*(2*Ny-1)+1:end);
leftSolution = reshape(leftSolution,[],Nx-1);
rightSolution = reshape(rightSolution,[],Nx);
container = zeros(2*Ny+1,2*Nx+1);
container(2:end-1,2:Nx) = leftSolution;
container(Ny+2:end-1,Nx+1:end-1) = rightSolution;
container(Ny+1:end-1,Nx+1:end-1);