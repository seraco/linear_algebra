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
A = -problem.coefficientsMatrixLshapeRegion();
eigenvalues = eig(A);
[lower,upper] = bandwidth(A);
bwdth = max(lower,upper);
sizeOfA = size(A);

%% 2)
problem = CSteadyBidimensionalHeat(5,5,1,1,1);
A = -problem.coefficientsMatrixLshapeRegion();
spy(A);

%% 3)
for i=[5,10,20]
    Nx = i;
    Ny = i;
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
    container = flipud(container);
    deltaX = a/Nx;
    deltaY = a/Ny;
    xVector = 0:deltaX:2*a;
    yVector = 0:deltaY:2*a;
    figure
    hold on
    [C,h] = contourf(xVector,yVector,container);
    clabel(C,h);
    xlabel('x');
    ylabel('y');
    title(strcat('Nx=',num2str(i),' Ny=',num2str(i)));
    hold off
end