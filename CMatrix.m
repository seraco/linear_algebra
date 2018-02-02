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
        function res = transpose(obj)
            res = zeros(obj.nColumns,obj.mRows);
            for i=1:obj.mRows
                for j=1:obj.nColumns
                    res(j,i) = obj.data(i,j);
                end
            end
            res = CMatrix(res);
        end
    end
    
end

