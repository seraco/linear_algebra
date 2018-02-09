classdef CUnidimensionalDiffusion
    % Class:    CUnidimensionalDiffusion
    %
    % Summary:  Class to describe the 1D diffusion problem. The problem is given
    %           by the diffusion equation which is discretised in the
    %           unidimensional domain with N+1 grid points. Boundary conditions
    %           have to be specified at the left and right edges.
    %
    % Methods:  CUnidimensionalDiffusion
    %               Constructor of the class.
    %           initialCondition
    %               Propagates the initial condition to start the simulation.
    %           boundaryConditions
    %               Enforces boundary conditions at the left and right edges.
    %           solve
    %               Solves the problem.
    %           analyticalSolution
    %               Analytical solution of the problem.
    %           derxxSecondOrder
    %               Returns a matrix of coefficients used to build the A matrix
    %               of the system of equations.

    % Version:  0.0.1
    % Author:   S.Ramon
    % License:  MIT

    properties
        nGridPoints             % number of grid points
        nTimeSteps              % number of time steps
        lSize                   % longitude of the domain
        totalTime               % total time of simulation
        diffusivity             % diffusivity constant
        initialConcentration    % initial concentration to start the simulation
        leftBoundaryCondition   % boundary condition on the left edge
        rightBoundaryCondition  % boundary condition on the right edge
    end

    methods
        function obj = CUnidimensionalDiffusion(Nx,Nt,L,t,D,Co,lBC,rBC)
            % Summary:  Constructor. Initialises the properties of the object.
            %
            % Args:     Nx
            %               Number of divisions of the domain in the x
            %               direction.
            %           Nt
            %               Number of time steps for the simulation.
            %           L
            %               Longitude of the domain.
            %           t
            %               Total time of simulation.
            %           D
            %               Diffusivity constant.
            %           Co
            %               Initial concentration.
            %           lBC
            %               Left boundary condition.
            %           rBC
            %               Right boundary condition.
            %
            % Returns:  obj
            %               CUnidimensionalDiffusion object with initilized
            %               properties.
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
            % Summary:  Propagates the initial condition to return a vector with
            %           the initial concentration values to start the
            %           simulation.
            %
            % Returns:  res
            %               Vector with initial concentration values.
            res = zeros(obj.nGridPoints,1);
            for i=1:obj.nGridPoints
                res(i,1) = obj.initialConcentration;
            end
            res = CVector(res);
        end
        function res = boundaryConditions(obj,concentration)
            % Summary:  Enforces the boundary conditions at the edges of the
            %           domain.
            %
            % Args:     concentration
            %               Vector of concentration values to be able to return
            %               a copy of this vector, but with the enforced
            %               boundary conditions.
            %
            % Returns:  res
            %               Concentration values with enforced boundary
            %               conditions.
            res = concentration;
            res.data(1) = obj.leftBoundaryCondition;
            res.data(obj.nGridPoints) = obj.rightBoundaryCondition;
        end
        function res = solve(obj)
            % Summary:  Solves the sistem along the several time steps. The
            %           system of equations encountered in every time step is
            %           solved using a direct solution optimized for tridiagonal
            %           matrices. It is defined in the class CLinearSystem.
            %
            % Returns:  res
            %               Vector of concentration values at the end of the
            %               simulation.
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
            % Summary:  Solves the problem with an analytical solution at the
            %           location x and at time t.
            %
            % Args:     x
            %               Space x coordinate at which the solution is
            %               computed.
            %           t
            %               Time in which the solution is computed.
            %
            % Returns:  res
            %               Concentration value.
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
    end

    methods(Static)
        function res = derxxSecondOrder(n)
            % Summary:  The function builds a matrix of the following form:
            %           B = [-2  1  0  0  0  0
            %                 1 -2  1  0  0  0
            %                 0  1 -2  1  0  0
            %                 0  0  1 -2  1  0
            %                 0  0  0  1 -2  1
            %                 0  0  0  0  1 -2]
            %
            % Args:     n
            %               Size of the matrix B.
            %
            % Returns:  res
            %               The matrix of the form described in B, but of size
            %               n.
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
