classdef CMatrix
    % Own implementation for matrices.
    %   The class is helpful for defining custom algorithms to handle
    %   calculations with matrices.
    
    % Version:  0.0.1
    % Author:   S.Ramon
    
    properties
        mRows
        nColumns
        data
    end
    
    methods
        function obj = CMatrix(matrix)
            obj.data = matrix;
            obj.mRows = size(matrix,1);
            obj.nColumns = size(matrix,2);
        end
        function res = plus(obj1,obj2)
            areEqualSize(obj1,obj2);
            res = zeros(obj1.mRows,obj1.nColumns);
            for i=1:obj1.mRows
                for j=1:obj1.nColumns
                    res(i,j) = obj1.data(i,j)+obj2.data(i,j);
                end
            end
            res = CMatrix(res);
        end
        function res = minus(obj1,obj2)
            areEqualSize(obj1,obj2);
            res = zeros(obj1.mRows,obj1.nColumns);
            for i=1:obj1.mRows
                for j=1:obj1.nColumns
                    res(i,j) = obj1.data(i,j)-obj2.data(i,j);
                end
            end
            res = CMatrix(res);
        end
        function res = mtimes(obj1,obj2)
            isMultiplicationPossible(obj1,obj2);
            res = zeros(obj1.mRows,obj2.nColumns);
            for i=1:obj1.mRows
                for j=1:obj1.nColumns
                    for k=1:obj2.nColumns
                        res(i,k) = res(i,k)+obj1.data(i,j)*obj2.data(j,k);
                    end
                end
            end
            res = CMatrix(res);
        end
        function areEqualSize(obj1,obj2)
            areRowsDifferent = obj1.mRows~=obj2.mRows;
            areColumnsDifferent = obj1.nColumns~=obj2.nColumns;
            if(areRowsDifferent || areColumnsDifferent)
                msg = ['The number of rows and columns should be the '...
                       'same in both matrices.'];
                error(msg);
            end
        end
        function isMultiplicationPossible(obj1,obj2)
            if(obj1.nColumns ~= obj2.mRows)
                msg = ['The number of columns in the first matrix '...
                       'should be equal to the number of rows in '...
                       'the second.'];
                error(msg);
            end
        end
        function res = infiniteNorm(obj)
            res = 0.0;
            for i=1:obj.mRows
                rowSum = 0.0;
                for j=1:obj.nColumns
                    rowSum = rowSum+abs(obj.data(i,j));
                end
                if(rowSum > res)
                    res = rowSum;
                end
            end
        end
        function res = absoluteValue(obj)
            res = zeros(obj.mRows,obj.nColumns);
            for i=1:obj.mRows
                for j=1:obj.nColumns
                    res(i,j) = abs(obj.data(i,j));
                end
            end
            res = CMatrix(res);
        end
        function res = isUpperTriangular(obj)
            for i=2:obj.mRows
                for j=1:i-1
                    if(obj.data(i,j) ~= 0)
                        res = false;
                        return;
                    end
                end
            end
            res = true;
        end
        function res = isLowerTriangular(obj)
            for i=1:obj.mRows-1
                for j=i+1:obj.nColumns
                    if(obj.data(i,j) ~= 0)
                        res = false;
                        return;
                    end
                end
            end
            res = true;
        end
        function res = isSymmetric(obj)
            if (obj.mRows ~= obj.nColumns)
                res = false;
                return;
            end
            for i=1:obj.mRows
                for j=i:obj.nColumns
                    if(obj.data(i,j) ~= obj.data(j,i))
                        res = false;
                        return;
                    end
                end
            end
            res = true;
        end
        function res = isTridiagonal(obj)
            if (obj.mRows ~= obj.nColumns)
                res = false;
                return;
            end
            for i=1:obj.mRows
                for j=1:i-2
                    if(obj.data(i,j) ~= 0)
                        res = false;
                        return;
                    end
                end
                for j=i+2:obj.nColumns
                    if(obj.data(i,j) ~= 0)
                        res = false;
                        return;
                    end
                end
            end
            res = true;
        end
        function res = transpose(obj)
            res = zeros(obj.nColumns,obj.mRows);
            for i=1:obj.mRows
                for j=1:obj.nColumns
                    res(j,i) = obj.data(i,j);
                end
            end
            res = CMatrix(res);
        end
        function res = getRow(obj,i)
            res = zeros(1,obj.nColumns);
            for j=1:obj.nColumns
                res(1,j) = obj.data(i,j);
            end
            res = CVector(res);
        end
        function res = getColumn(obj,j)
            res = zeros(obj.mRows,1);
            for i=1:obj.nColumns
                res(i,1) = obj.data(i,j);
            end
            res = CVector(res);
        end
        function res = swapSubrows(obj,i1,i2,j1,j2)
            res = obj.data;
            for j=j1:j2
                res(i1,j)=obj.data(i2,j);
                res(i2,j)=obj.data(i1,j);
            end
            res = CMatrix(res);
        end
        function res = getDiagonal(obj)
            if (obj.mRows ~= obj.nColumns)
                error('Number of rows and columns should be equal.')
            end
            iIndices = zeros(obj.mRows,1);
            jIndices = zeros(obj.mRows,1);
            values = zeros(obj.mRows,1);
            for i=1:obj.mRows
                iIndices(i) = i;
                jIndices(i) = i;
                values(i) = obj.data(i,i);
            end
            res = sparse(iIndices,jIndices,values);
            res = CMatrix(res);
        end
        function res = getDiagonalInverse(obj)
            if (obj.mRows ~= obj.nColumns)
                error('Number of rows and columns should be equal.')
            end
            iIndices = zeros(obj.mRows,1);
            jIndices = zeros(obj.mRows,1);
            values = zeros(obj.mRows,1);
            for i=1:obj.mRows
                iIndices(i) = i;
                jIndices(i) = i;
                values(i) = 1/obj.data(i,i);
            end
            res = sparse(iIndices,jIndices,values);
            res = CMatrix(res);
        end
        function res = getLower(obj)
            if (obj.mRows ~= obj.nColumns)
                error('Number of rows and columns should be equal.')
            end
            iIndices = zeros(obj.mRows,1);
            jIndices = zeros(obj.mRows,1);
            values = zeros(obj.mRows,1);
            k = 0;
            for i=1:obj.mRows
                for j=1:i-1
                    k = k+1;
                    iIndices(k) = i;
                    jIndices(k) = j;
                    values(k) = obj.data(i,j);
                end
            end
            iIndices(k+1) = obj.nColumns;
            jIndices(k+1) = obj.nColumns;
            values(k+1) = 0;
            res = sparse(iIndices,jIndices,values);
            res = CMatrix(res);
        end
        function res = getUpper(obj)
            if (obj.mRows ~= obj.nColumns)
                error('Number of rows and columns should be equal.')
            end
            iIndices = zeros(obj.mRows,1);
            jIndices = zeros(obj.mRows,1);
            values = zeros(obj.mRows,1);
            k = 0;
            for j=1:obj.nColumns
                for i=1:j-1
                    k = k+1;
                    iIndices(k) = i;
                    jIndices(k) = j;
                    values(k) = obj.data(i,j);
                end
            end
            iIndices(k+1) = obj.nColumns;
            jIndices(k+1) = obj.nColumns;
            values(k+1) = 0;
            res = sparse(iIndices,jIndices,values);
            res = CMatrix(res);
        end
        function res = findMaximumMagnitude(obj)
            res = 0;
            for i=1:obj.mRows
                for j=1:obj.nColumns
                    absElement = abs(obj.data(i,j));
                    if(absElement > res)
                        res = absElement;
                    end
                end
            end
        end
    end
    
    methods(Static)
        function res = identityMatrix(n)
            res = zeros(n,n);
            for i=1:n
                for j=1:n
                    if(i == j)
                        res(i,j) = 1.0;
                    else
                        res(i,j) = 0.0;
                    end
                end
            end
            res = CMatrix(res);
        end
        function res = permutationMatrix(i1,i2,n)
            identity = CMatrix.identityMatrix(n);
            e1 = identity.getColumn(i1);
            e2 = identity.getColumn(i2);
            e2me1 = e2-e1;
            outerProduct = e2me1*e2me1;
            res = zeros(n,n);
            for i=1:n
                for j=1:n
                    res(i,j) = identity.data(i,j)-outerProduct.data(i,j);
                end
            end
            res = CMatrix(res);
        end
    end
    
end
