%% Script to solve Question 3.
%   The question is about a chemical concentration problem along a porous
%   medium. The medium with initial concentration Co is exposed to an
%   environment which has concentration 0. This script calculates the
%   matrix of coefficients B. Finally, it compares the numerical solution 
%   with the analytical for different discretisations.

% Version:  0.0.1
% Author:   S.Ramon

addpath('../LinearAlgebra','../PhysicalProblems');

%% 3)
n = 5;
B = CUnidimensionalDiffusion.derxxSecondOrder(5);
B = B.data*-1;

%% 4)
N = 40;
Ntime = 6;
D = 1;
t = 0.1;
L = 1;
Co = 1;
BCleft = 0;
BCright = 0;
problem = CUnidimensionalDiffusion(N,Ntime,D,t,L,Co,BCleft,BCright);
numerical = problem.solve();
deltaX = L/N;
deltaT = t/Ntime;
xVector = (0:deltaX:L)';
analytical = zeros(N+1,1);
tVector = zeros(Ntime,1);
jTime = 0;
for j=1:Ntime
    jTime = jTime+deltaT;
    tVector(j) = jTime;
    for i=1:N+1
        analytical(i,j) = problem.analyticalSolution(xVector(i),jTime);
    end
end
for i=1:Ntime
    figure
    hold on
    plot(xVector,analytical(:,i));
    plot(xVector,numerical(i).data);
    legend('analytical','numerical');
    xlabel('x');
    ylabel('C');
    title(strcat('t=',num2str(tVector(i))));
    hold off
end