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
            res = sparse(nHorLeft*nVerLeft+nHorRight*nVerRight);
            for j=1:nHorLeft
                for i=1:nVerLeft
                    k = (j-1)*(nVerLeft)+i;
                    res(k,k) = -2/deltaXsquared-2/deltaYsquared;
                    if(i ~= 1)
                        res(k,k-1) = 1/deltaYsquared;
                    end
                    if(i ~= nVerLeft)
                        res(k,k+1) = 1/deltaYsquared;
                    end
                    if(j ~= 1)
                        res(k,k-nVerLeft) = 1/deltaXsquared;
                    end
                    if(j ~= nHorLeft)
                        res(k,k+nVerLeft) = 1/deltaXsquared;
                    end
                end
            end
            kTotalLeft = k;
            for j=1:nHorRight
                for i=1:nVerRight
                    k = kTotalLeft+(j-1)*(nVerRight)+i;
                    res(k,k) = -2/deltaXsquared-2/deltaYsquared;
                    res(k,k-nVerRight) = 1/deltaXsquared;
                    if(i ~= 1)
                        res(k,k-1) = 1/deltaYsquared;
                    end
                    if(i ~= nVerRight)
                        res(k,k+1) = 1/deltaYsquared;
                    end
                    if(j ~= nHorRight)
                        res(k,k+nVerRight) = 1/deltaXsquared;
                    end
                end
            end
        end
        function res = solve(obj)
            A = coefficientsMatrix(obj);
            b = obj.appliedHeat*ones(size(A,1),1);
            res = b;
        end
    end

end
