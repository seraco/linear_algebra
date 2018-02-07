classdef CUnidimensionalDiffusion
    % Class to describe the 1D diffusion problem.
    %   The problem is given by the diffusion equation which is discretised
    %   in the unidimensional domain with N+1 grid points. Boundary
    %   conditions have to be specified at the left and right edges.

    % Version:  0.0.1
    % Author:   S.Ramon

    properties
        nGridPoints
        nTimeSteps
        lSize
        totalTime
        diffusivity
        initialConcentration
        leftBoundaryCondition
        rightBoundaryCondition
    end

    methods
        function obj = CUnidimensionalDiffusion(Nx,Nt,L,t,D,Co,lBC,rBC)
            obj.nGridPoints = Nx+1;
            obj.nTimeSteps = Nt;
            obj.lSize = L;
            obj.totalTime = t;
            obj.diffusivity = D;
            obj.initialConcentration = Co;
            obj.leftBoundaryCondition = lBC;
            obj.rightBoundaryCondition = rBC;
        end
        function res = initialCondition(obj)
            res = zeros(obj.nGridPoints,1);
            for i=1:obj.nGridPoints
                res(i,1) = obj.initialConcentration;
            end
            res = CVector(res);
        end
        function res = boundaryConditions(obj,concentration)
            res = concentration;
            res.data(1) = obj.leftBoundaryCondition;
            res.data(obj.nGridPoints) = obj.rightBoundaryCondition;
        end
        function res = solve(obj)
            deltaSpace = obj.lSize/(obj.nGridPoints-1);
            deltaTime = obj.totalTime/obj.nTimeSteps;
            sigmaParameter = obj.diffusivity*deltaTime/deltaSpace^2;
            identity = CMatrix.identityMatrix(obj.nGridPoints);
            derxx = CUnidimensionalDiffusion...
                    .derxxSecondOrder(obj.nGridPoints);
            B = derxx.data*-1;
            b = identity.data-sigmaParameter/2*B;
            A = identity.data+sigmaParameter/2*B;
            b = CMatrix(b);
            A = CMatrix(A);
            concentrationN = initialCondition(obj);
            concentrationN = boundaryConditions(obj,concentrationN);
            for i=1:obj.nTimeSteps
                concentrationN = CMatrix(concentrationN.data);
                rhs = b*concentrationN;
                sys = CLinearSystem(A.data,rhs.data);
                concentrationN1 = sys.tridiagonalDirectSolution();
                concentrationN1 = boundaryConditions(obj,concentrationN1);
                concentrationN = concentrationN1;
                res(i) = concentrationN;
            end
        end
        function res = analyticalSolution(obj,x,t)
            Co = obj.initialConcentration;
            n = 100;
            D = obj.diffusivity;
            L = obj.lSize;
            res = 0.0;
            for i=1:n
                factor = (2*i-1)*pi;
                res = res+1/factor*exp(-(D*factor^2*t/L^2))*sin(factor*x/L);
            end
            res = 4*Co*res;
        end
        function compareSolutions(obj)
        end
    end

    methods(Static)
        function res = derxxSecondOrder(n)
            % The function builds a matrix of the form
            %       B = [
            %           -2  1  0  0  0  0
            %            1 -2  1  0  0  0
            %            0  1 -2  1  0  0
            %            0  0  1 -2  1  0
            %            0  0  0  1 -2  1
            %            0  0  0  0  1 -2
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
                        values(k) = -2.0;
                        k = k+1;
                    elseif(i==j+1 || i==j-1)
                        iIndices(k) = i;
                        jIndices(k) = j;
                        values(k) = 1.0;
                        k = k+1;
                    end
                end
            end
            res = sparse(iIndices,jIndices,values);
            res = CMatrix(res);
        end
    end

end
