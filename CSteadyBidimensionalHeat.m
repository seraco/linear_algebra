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
        end
    end

end
