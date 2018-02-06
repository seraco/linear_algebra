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
                       'coefficients should be upper triangular.'];
                error(msg);
            end
            res = zeros(size(obj.b.data));
            for i=obj.b.nElements:-1:1
                for j=obj.b.nElements:-1:i+1
                    res(i)=res(i)+obj.A.data(i,j)*res(j);
                end
                res(i)=(obj.b.data(i)-res(i))/obj.A.data(i,i);
            end
            res = CVector(res);
        end
        function res = tridiagonalBackwardSubstitution(obj)
            if(~obj.A.isUpperTriangular)
                msg = ['For backward substitution the matrix of '...
                       'coefficients should be upper triangular.'];
                error(msg);
            end
            nElements = obj.b.nElements;
            res = zeros(size(obj.b.data));
            res(nElements) = obj.b.data(nElements)/...
                                   obj.A.data(nElements,nElements);
            for i=nElements-1:-1:1
                res(i) = res(i)+obj.A.data(i,i+1)*res(i+1);
                res(i) = (obj.b.data(i)-res(i))/obj.A.data(i,i);
            end
            res = CVector(res);
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
                    res(i) = res(i)+obj.A.data(i,j)*res(j);
                end
                res(i) = (obj.b.data(i)-res(i))/obj.A.data(i,i);
            end
            res = CVector(res);
        end
        function res = tridiagonalForwardSubstitution(obj)
            if(~obj.A.isLowerTriangular)
                msg = ['For forward substitution the matrix of '...
                       'coefficients should be lower triangular'];
                error(msg);
            end
            res = zeros(size(obj.b.data));
            res(1) = obj.b.data(1)/obj.A.data(1,1);
            for i=2:obj.b.nElements
                res(i) = res(i)+obj.A.data(i,i-1)*res(i-1);
                res(i) = (obj.b.data(i)-res(i))/obj.A.data(i,i);
            end
            res = CVector(res);
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
                gv = CLinearSystem.gaussVector(U,k);
                e = I.getColumn(k);
                LE = gv*e;
                J = I-LE;
                U = J*U;
                L = L+LE;
            end
        end
        function res = growthFactor(obj)
            [~,~,U] = luWithPartialPivoting(obj);
            maxA=0.0;
            maxU=0.0;
            for i=1:obj.A.mRows
                for j=1:obj.A.mRows
                    absA = abs(obj.A.data(i,j));
                    absU = abs(U.data(i,j));
                    if(absA > maxA)
                        maxA = absA;
                    end
                    if(absU > maxU)
                        maxU = absU;
                    end
                end
            end
            res = maxU/maxA;
        end
        function res = directSolution(obj,type)
            rhs = CMatrix(obj.b.data);
            switch type
                case 'lu'
                    [P,L,U] = luWithPartialPivoting(obj);
                    c = P*rhs;
                case 'cholesky'
                    [L,U] = choleskyDecomposition(obj);
                    c = rhs;
                otherwise
                    error('Unknown decomposition type.');
            end
            firstSystem = CLinearSystem(L.data,c.data);
            y = firstSystem.forwardSubstitution();
            secondSystem = CLinearSystem(U.data,y.data);
            res = secondSystem.backwardSubstitution();
        end
        function [L,U] = choleskyDecomposition(obj)
            msg = 'Matrix of coefficients should be symmetric.';
            if(~obj.A.isSymmetric)
                error(msg);
            end
            matrix = obj.A.data;
            mRows = obj.A.mRows;
            nColumns = obj.A.nColumns;
            U = CMatrix(zeros(mRows,nColumns));
            for j=1:nColumns
                s = 0.0;
                for k=1:j-1
                    kU = U.getColumn(k);
                    jU = U.getColumn(j);
                    kU = kU.getSubvector(1,k-1);
                    jU = jU.getSubvector(1,k-1);
                    t = matrix(k,j)-kU.*jU;
                    t = t/U.data(k,k);
                    U.data(k,j) = t;
                    s = s+t*t;
                end
                s = matrix(j,j)-s;
                if s<0
                    error(['Matrix of coefficients should be definite '...
                           'positive.']);
                end
                U.data(j,j) = sqrt(s);
            end
            L = U.transpose;
        end
        function [L,U] = tridiagonalCholesky(obj)
            msg = 'Matrix of coefficients should be symmetric.';
            if(~obj.A.isSymmetric)
                error(msg);
            end
            msg1 = 'Matrix of coefficients should be tridiagonal.';
            if(~obj.A.isTridiagonal)
                error(msg1);
            end
            matrix = obj.A.data;
            mRows = obj.A.mRows;
            iIndices = zeros(1,mRows+mRows-1);
            jIndices = zeros(1,mRows+mRows-1);
            values = zeros(1,mRows+mRows-1);
            iIndices(1) = 1;
            jIndices(1) = 1;
            values(1) = sqrt(matrix(1,1));
            for i=2:mRows
                iIndices(2*i-2) = i;
                jIndices(2*i-2) = i-1;
                values(2*i-2) = matrix(i,i-1)/values(2*i-3);
                iIndices(2*i-1) = i;
                jIndices(2*i-1) = i;
                values(2*i-1) = sqrt(matrix(i,i)-values(2*i-2)^2);
            end
            values = real(values);
            L = sparse(iIndices,jIndices,values);
            L = CMatrix(L);
            U = L.transpose;
        end
        function res = tridiagonalDirectSolution(obj)
            [L,U] = tridiagonalCholesky(obj);
            firstSystem = CLinearSystem(L.data,obj.b.data);
            y = firstSystem.tridiagonalForwardSubstitution();
            secondSystem = CLinearSystem(U.data,y.data);
            res = secondSystem.tridiagonalBackwardSubstitution();
        end
    end
    
    methods(Static)
        function res = gaussVector(matrix,k)
            res = zeros(matrix.mRows-k,1);
            for i=k+1:matrix.mRows
                res(i) = matrix.data(i,k)/matrix.data(k,k);
            end
            res = CVector(res);
        end
    end
    
end

