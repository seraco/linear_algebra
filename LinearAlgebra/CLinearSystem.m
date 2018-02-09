classdef CLinearSystem
    % Class:    CLinearSystem
    %
    % Summary:  Class defined for handling linear systems of equations.
    %           This implementation is intended for performing calculations with
    %           linear systems.
    %
    % Methods:  CLinearSystem
    %               Constructor of the class.
    %           isSystemValid
    %               Checks if system of equation is valid.
    %           backwardSubstitution
    %               Backward substitution algorithm.
    %           tridiagonalBackwardSubstitution
    %               Backward substitution improved for tridiagonal matrices.
    %           forwardSubstitution
    %               Forward substitution algorithm.
    %           tridiagonalForwardSubstitution
    %               Forward substitution improved for tridiagonal matrices.
    %           luWithPartialPivoting
    %               LU decomposition process with partial pivoting.
    %           growthFactor
    %               Calculates the growth factor of A matrix.
    %           directSolution
    %               Solves the system with a direct method.
    %           choleskyDecomposition
    %               Cholesky decomposition process.
    %           tridiagonalCholesky
    %               Cholesky decomposition improved for tridiagonal matrices.
    %           tridiagonalDirectSolution
    %               Direct solution improved for tridiagonal matrices.
    %           iterativeSolution
    %               Solves the system with an iterative method.
    %           jacobiMethod
    %               Jacobi method iterative algorithm.
    %           gaussVector
    %               Computes the Gauss vector for Gaussian elimination.

    % Version:  0.0.1
    % Author:   S.Ramon
    % License:  MIT

    properties
        A   % matrix of coefficients
        b   % RHS vector
    end

    methods
        function obj = CLinearSystem(A,b)
            % Summary:  Constructor. Initialises the properties of the object
            %           and checks that the system of equations is valid.
            % Args:     A
            %               Matrix of coefficients.
            %           b
            %               RHS vector.
            % Returns:  obj
            %               CLinearSystem object with initilized properties.
            obj.A = CMatrix(A);
            obj.b = CVector(b);
            isSystemValid(obj);
        end
        function isSystemValid(obj)
            % Summary:  Method that checks if a system of equations is valid by
            %           evaluating if the matrix of coefficients is square and
            %           if the size of the RHS vector agrees with the size that
            %           matrix.
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
            % Summary:  Backward substitution algorithm. If the matrix of
            %           coefficients of the system is upper triangular, the
            %           method calculates each unknown of the system from the
            %           bottom to the top.
            % Returns:  res
            %               Vector of unknowns.
            if(~obj.A.isUpperTriangular)
                msg = ['For backward substitution the matrix of '...
                       'coefficients should be upper triangular.'];
                error(msg);
            end
            res = zeros(size(obj.b.data));
            for i=obj.b.nElements:-1:1
                for j=obj.b.nElements:-1:i+1
                    res(i) = res(i)+obj.A.data(i,j)*res(j);
                end
                res(i) = (obj.b.data(i)-res(i))/obj.A.data(i,i);
            end
            res = CVector(res);
        end
        function res = tridiagonalBackwardSubstitution(obj)
            % Summary:  The same algorithm as in backwardSubstitution() but
            %           modified to be more efficient for tridiagonal matrices.
            % Returns:  res
            %               Vector of unknowns.
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
            % Summary:  Forward substitution algorithm. If the matrix of
            %           coefficients of the system is lower triangular, the
            %           method calculates each unknown of the system from the
            %           top to the bottom.
            % Returns:  res
            %               Vector of unknowns.
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
            % Summary:  The same algorithm as in forwardSubstitution() but
            %           modified to be more efficient for tridiagonal matrices.
            % Returns:  res
            %               Vector of unknowns.
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
            % Summary:  LU decomposition with partial pivoting. The procedure
            %           of partial pivoting avoids roundoff errors and dividing
            %           by zero, which is absolutely necessary when using a
            %           computer to perform the LU factorisation method.
            % Returns:  P
            %               Permutation matrix.
            %           L
            %               Lower triangular matrix from LU decomposition.
            %           U
            %               Upper triangular matrix from LU decomposition.
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
            % Summary:  The growth factor is calculated here. The maximum
            %           magnitude is found for both A and U matrices and their
            %           quotient gives the growth factor.
            % Returns:  res
            %               Growth factor.
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
            % Summary:  This method gives the direct solution of the system of
            %           equations by combining a factorization method (LU
            %           decomposition, Cholesky factorization) with forward and
            %           backward substitution methods to get the vector of
            %           unknowns.
            % Args:     type
            %               Type of decomposition of A for the direct solution.
            % Returns:  res
            %               Vector of unknowns.
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
            % Summary:  Cholesky decomposition is a more efficient way to
            %           factorize A matrix in LU for symmetric positive definite
            %           matrices.
            % Returns:  L
            %               Lower triangular matrix from the decomposition.
            %           U
            %               Upper triangular matrix from the decomposition.
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
            % Summary:  Cholesky decomposition improved for tridiagonal
            %           matrices. The procedure takes advantage of the fact that
            %           L will be lower banded matrix of bandwidth 1. The
            %           algorithm can be checked at
            %           http://www.cs.cornell.edu/courses/cs4220/2014sp/CVLBook/chap7.pdf
            %           in chapter 7.3.2.
            % Returns:  L
            %               Lower triangular matrix from the decomposition.
            %           U
            %               Upper triangular matrix from the decomposition.
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
            % Summary:  Direct solution improved for tridiagonal matrices with
            %           Cholesky factorization.
            % Returns:  res
            %               Vector of unknowns.
            [L,U] = tridiagonalCholesky(obj);
            firstSystem = CLinearSystem(L.data,obj.b.data);
            y = firstSystem.tridiagonalForwardSubstitution();
            secondSystem = CLinearSystem(U.data,y.data);
            res = secondSystem.tridiagonalBackwardSubstitution();
        end
        function res = iterativeSolution(obj,type,tolerance)
            % Summary:  This method gives the iterative solution of the system
            %           of equations.
            % Returns:  res
            %               Vector of unknowns.
            switch type
                case 'jacobi'
                    res = jacobiMethod(obj,tolerance);
                otherwise
                    error('Unknown decomposition type.');
            end
        end
        function res = jacobiMethod(obj,tolerance)
            % Summary:  The Jacobi method is implemented by decomposing A matrix
            %           in its lower, diagonal and upper parts, and using sparse
            %           operations to operate more efficiently in matrices where
            %           most of the entries are zero.
            % Args:     tolerance
            %               Tolerance margin for convergence.
            % Returns:  res
            %               Vector of unknowns.
            convergenceError = tolerance;
            initGuess = CVector(obj.b.data);
            Dinv = obj.A.getDiagonalInverse;
            L = obj.A.getLower;
            U = obj.A.getUpper;
            LplusU = L+U;
            rhs = CMatrix(obj.b.data);
            xK = CMatrix(initGuess.data);
            res = CMatrix(zeros(obj.b.nElements,1));
            iter = 0;
            while(true)
                iter = iter+1;
                factor = rhs-LplusU*xK;
                for i=1:obj.b.nElements
                    res.data(i,1) = Dinv.data(i,i)*factor.data(i,1);
                end
                err = res-xK;
                maxErrorMagnitude = err.findMaximumMagnitude;
                if(maxErrorMagnitude < convergenceError)
                    break;
                end
                xK = res;
            end
            disp(iter);
            res = CVector(res.data);
        end
    end

    methods(Static)
        function res = gaussVector(matrix,k)
            % Summary:  Returns the Gauss vector of matrix for the k step of
            %           Gaussian elimination.
            % Args:     matrix
            %               Matrix at k step of Gaussian elimination.
            %           k
            %               Step of Gaussian elimination.
            % Returns:  res
            %               Gauss vector.
            res = zeros(matrix.mRows-k,1);
            for i=k+1:matrix.mRows
                res(i) = matrix.data(i,k)/matrix.data(k,k);
            end
            res = CVector(res);
        end
    end

end
