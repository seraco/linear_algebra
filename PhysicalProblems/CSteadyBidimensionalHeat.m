classdef CSteadyBidimensionalHeat
    % Class to describe the 2D steady heat problem.
    %   The problem is given by the heat equation which is discretised
    %   in the bidimensional domain with Nx+1 grid points in x direction
    %   and Ny+1 grid points in the y direction. Boundary conditions have
    %   to be specified at the edges of the domain.

    % Version:  0.0.1
    % Author:   S.Ramon

    properties
        nHoritzontalDivisions
        nVerticalDivisions
        lSize
        appliedHeat
    end

    methods
        function obj = CSteadyBidimensionalHeat(Nx,Ny,L,q)
            obj.nHoritzontalDivisions = Nx;
            obj.nVerticalDivisions = Ny;
            obj.lSize = L;
            obj.appliedHeat = q;
        end
        function res = coefficientsMatrix(obj)
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
                end
            end
            kTotalLeft = k;
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
            A = coefficientsMatrix(obj);
            b = obj.appliedHeat*ones(size(A,1),1);
            sys = CLinearSystem(A,b);
            res = sys.iterativeSolution('jacobi',0.000001);
        end
    end

end
