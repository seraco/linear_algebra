classdef CSteadyBidimensionalHeat
    % Class:    CSteadyBidimensionalHeat
    %
    % Summary:  Class to describe the 2D steady heat problem. The problem is
    %           given by the heat equation which is discretised in the
    %           bidimensional domain with Nx+1 grid points in x direction and
    %           Ny+1 grid points in the y direction. Boundary conditions have to
    %           be specified at the edges of the domain.
    %
    % Methods:  CSteadyBidimensionalHeat
    %               Constructor of the class.
    %           coefficientsMatrixLshapeRegion
    %               Builds the matrix of coefficients in the LHS of the system
    %               to solve.
    %           solve
    %               Solves the heat problem.

    % Version:  0.0.1
    % Author:   S.Ramon
    % License:  MIT

    properties
        nHoritzontalDivisions   % number of division in the x direction
        nVerticalDivisions      % number of division in the y direction
        lSize                   % characteristic longitude
        appliedHeat             % applied heat
        conductivity            % thermal conductivity
    end

    methods
        function obj = CSteadyBidimensionalHeat(Nx,Ny,L,q,k)
            % Summary:  Constructor. Initialises the properties of the object.
            %
            % Args:     Nx
            %               Number of divisions of the domain in the x
            %               direction.
            %           Ny
            %               Number of divisions of the domain in the y
            %               direction.
            %           L
            %               Characteristic longitude of the domain.
            %           q
            %               Applied heat to the medium.
            %           k
            %               Thermal conductivity.
            %
            % Returns:  obj
            %               CUnidimensionalDiffusion object with initilized
            %               properties.
            obj.nHoritzontalDivisions = Nx;
            obj.nVerticalDivisions = Ny;
            obj.lSize = L;
            obj.appliedHeat = q;
            obj.conductivity = k;
        end
        function res = coefficientsMatrixLshapeRegion(obj)
            % Summary:  Computes the matrix of coefficients that is encountered
            %           in the LHS of the system of equations to solve but
            %           multiplied by -1. This matrix is computed for a L shaped
            %           domain with the following shape:
            %            ----------
            %           |  1    6  |
            %           |  2    7  |
            %           |  3    8  |__________
            %           |  4    9  11  13  15 |
            %           |  5    10 12  14  16 |
            %           |_____________________|
            %           The size of the left rectangle is (Nx+1)x(2Ny+1) and the
            %           size of the right rectangle is (Nx+1)x(Ny+1). The
            %           numbers in the figure are the nodes used for the
            %           discretisation. Note that the numbering goes from top to
            %           bottom and left to right. Besides, the edges of the
            %           domain are not introduced as nodes because their
            %           temperature is kept as zero for this problem.
            %
            % Returns:  res
            %               Matrix of coefficients.
            deltaXsquared = (obj.lSize/obj.nHoritzontalDivisions)^2;
            deltaYsquared = (obj.lSize/obj.nVerticalDivisions)^2;
            nHorLeft = obj.nHoritzontalDivisions-1;
            nVerLeft = 2*obj.nVerticalDivisions-1;
            nHorRight = obj.nHoritzontalDivisions;
            nVerRight = obj.nVerticalDivisions-1;
            c1 = 1/deltaXsquared;
            c2 = 1/deltaYsquared;
            c3 = -2*(c1+c2);
            res = sparse(nHorLeft*nVerLeft+nHorRight*nVerRight);
            % loop for the coefficients in the left rectangle
            for j=1:nHorLeft
                for i=1:nVerLeft
                    k = (j-1)*(nVerLeft)+i;
                    res(k,k) = c3;
                    if(i ~= 1)
                        res(k,k-1) = c2;
                    end
                    if(i ~= nVerLeft)
                        res(k,k+1) = c2;
                    end
                    if(j ~= 1)
                        res(k,k-nVerLeft) = c1;
                    end
                    if(j ~= nHorLeft)
                        res(k,k+nVerLeft) = c1;
                    end
                    if((j == nHorLeft) && (i>nVerLeft-nVerRight))
                        res(k,k+nVerRight) = c1;
                    end
                end
            end
            kTotalLeft = k;
            % loop for the coefficients in the right rectangle
            for j=1:nHorRight
                for i=1:nVerRight
                    k = kTotalLeft+(j-1)*(nVerRight)+i;
                    res(k,k) = c3;
                    res(k,k-nVerRight) = c1;
                    if(i ~= 1)
                        res(k,k-1) = c2;
                    end
                    if(i ~= nVerRight)
                        res(k,k+1) = c2;
                    end
                    if(j ~= nHorRight)
                        res(k,k+nVerRight) = c1;
                    end
                end
            end
        end
        function res = solve(obj)
            % Summary:  Solves the problem of heat conduction using the Jacobi
            %           method iterative solver.
            %
            % Returns:  res
            %               Temperature values at the different nodes in the
            %               domain.
            A = -coefficientsMatrixLshapeRegion(obj);
            b = obj.appliedHeat/obj.conductivity*ones(size(A,1),1);
            sys = CLinearSystem(A,b);
            res = sys.iterativeSolution('jacobi',0.000001);
        end
    end

end
