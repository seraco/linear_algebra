classdef CLinearSystem
    % Class defined for handling linear systems of equations.
    %   This implementation is intended for performing calculations with
    %   linear systems.
    
    % Version:  0.0.1
    % Author:   S.Ramon
    
    properties
        A
        b
    end
    
    methods
        function obj = CLinearSystem(A,b)
            obj.A = CMatrix(A);
            obj.b = CVector(b);
            isSystemValid(obj);
        end
        function isSystemValid(obj)
            isSquareMatrix = obj.A.mRows==obj.A.nColumns;
            isCompleteVector = obj.A.mRows==obj.b.nElements;
            if(~isSquareMatrix)
                msg = 'Matrix of coefficients should be square matrix.';
                error(msg);
            end
            if(~isCompleteVector)
                msg = ['Matrix of coefficients and RHS vector should '...
                       'be of equivalent size.'];
                error(msg);
            end
        end
        function res = backwardSubstitution(obj)
            if(~obj.A.isUpperTriangular)
                msg = ['For backward substitution the matrix of '...
                       'coefficients should be upper triangular'];
                error(msg);
            end
            res = zeros(size(obj.b.data));
            for i=obj.b.nElements:-1:1
                for j=obj.b.nElements:-1:i+1
                    res(i)=res(i)+obj.A.data(i,j)*res(j);
                end
                res(i)=(obj.b.data(i)-res(i))/obj.A.data(i,i);
            end
        end
        function res = forwardSubstitution(obj)
            if(~obj.A.isLowerTriangular)
                msg = ['For forward substitution the matrix of '...
                       'coefficients should be lower triangular'];
                error(msg);
            end
            res = zeros(size(obj.b.data));
            for i=1:obj.b.nElements
                for j=1:i-1
                    res(i)=res(i)+obj.A.data(i,j)*res(j);
                end
                res(i)=(obj.b.data(i)-res(i))/obj.A.data(i,i);
            end
        end
        function [P,L,U] = luWithPartialPivoting(obj)
            L = CMatrix.identityMatrix(obj.A.mRows);
            I = CMatrix.identityMatrix(obj.A.mRows);
            P = CMatrix.identityMatrix(obj.A.mRows);
            U = obj.A;
            for k=1:obj.A.mRows-1
                columnU = U.getColumn(k);
                columnU = columnU.getSubvector(k,obj.A.mRows);
                iMax = columnU.maxElementIndex()+k-1;
                permutation = ...
                    CMatrix.permutationMatrix(iMax,k,obj.A.mRows);
                P = permutation*P;
                U = U.swapSubrows(iMax,k,k,obj.A.mRows);
                L = L.swapSubrows(iMax,k,1,k-1);
                gv = U.gaussVector(k);
                e = I.getColumn(k);
                LE = gv*e;
                J = I-LE;
                U = J*U;
                L = L+LE;
            end
        end
    end
    
end

