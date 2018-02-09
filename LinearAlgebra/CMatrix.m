classdef CMatrix
    % Class:    CMatrix
    %
    % Summary:  Own implementation for matrices. The class is helpful for
    %           defining custom algorithms to handle calculations with matrices.
    %
    % Methods:  CMatrix
    %               Constructor of the class.
    %           plus
    %               Overloading of + operator.
    %           minus
    %               Overloading of - operator.
    %           mtimes
    %               Overloading of * operator.
    %           areEqualSize
    %               Checks that two CMatrix objects are of equal size.
    %           isMultiplicationPossible
    %               Checks that multiplication between two CMatrix objects is
    %               possible.
    %           infiniteNorm
    %               Calculates the infinite norm.
    %           absoluteValue
    %               Returns the matrix with all its entries with absolute value.
    %           isUpperTriangular
    %               Checks if the matrix is upper triangular.
    %           isLowerTriangular
    %               Checks if the matrix is lower triangular.
    %           isSymmetric
    %               Checks if the matrix is symmetric.
    %           isTridiagonal
    %               Checks if the matrix is tridiagonal.
    %           transpose
    %               Returns the transpose of the matrix.
    %           getRow
    %               Returns the specified row of the matrix.
    %           getColumn
    %               Returns the specified column of the matrix.
    %           swapSubrows
    %               Interchanges two subrows of the matrix.
    %           getDiagonal
    %               Returns the diagonal of the matrix.
    %           getDiagonalInverse
    %               Returns the inverse of the matrix returned by getDiagonal().
    %           getLower
    %               Returns the lower part of the matrix.
    %           getUpper
    %               Returns the upper part of the matrix.
    %           findMaximumMagnitude
    %               Returns the entry with maximum magnitude.
    %           identityMatrix
    %               Returns the identity matrix of specified size.
    %           permutationMatrix
    %               Returns a permutation matrix to swap rows or columns.

    % Version:  0.0.1
    % Author:   S.Ramon
    % License:  MIT

    properties
        mRows       % number of rows of the matrix
        nColumns    % number of columns of the matrix
        data        % stores the matrix data
    end

    methods
        function obj = CMatrix(matrix)
            % Summary:  Constructor. Initialises the properties of the object.
            %
            % Args:     matrix
            %               Matrix that is stored as data.
            %
            % Returns:  obj
            %               CMatrix object with initilized properties.
            obj.data = matrix;
            obj.mRows = size(matrix,1);
            obj.nColumns = size(matrix,2);
        end
        function res = plus(obj1,obj2)
            % Summary:  Overloading of + operator.
            %
            % Args:     obj1
            %               CMatrix object at on the left of operator +.
            %           obj2
            %               CMatrix object at on the right of operator +.
            %
            % Returns:  res
            %               CMatrix object with summed data entries.
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
            % Summary:  Overloading of - operator.
            %
            % Args:     obj1
            %               CMatrix object at on the left of operator -.
            %           obj2
            %               CMatrix object at on the right of operator -.
            %
            % Returns:  res
            %               CMatrix object with substracted data entries.
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
            % Summary:  Overloading of * operator.
            %
            % Args:     obj1
            %               CMatrix object at on the left of operator *.
            %           obj2
            %               CMatrix object at on the right of operator *.
            %
            % Returns:  res
            %               CMatrix object obtained after multiplication.
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
            % Summary:  Checks that the data of two CMatrix objects is of equal
            %           size. If the number of rows or columns are different,
            %           displays an error message.
            %
            % Args:     obj1
            %               First CMatrix object.
            %           obj2
            %               Second CMatrix object.
            areRowsDifferent = obj1.mRows~=obj2.mRows;
            areColumnsDifferent = obj1.nColumns~=obj2.nColumns;
            if(areRowsDifferent || areColumnsDifferent)
                msg = ['The number of rows and columns should be the '...
                       'same in both matrices.'];
                error(msg);
            end
        end
        function isMultiplicationPossible(obj1,obj2)
            % Summary:  Checks that the multiplication between two CMatrix
            %           objects is possible. If the number of columns from the
            %           first matrix is equal to the number of rows from the
            %           second matrix, then it is possible. If not, it displays
            %           an error message.
            %
            % Args:     obj1
            %               First CMatrix object.
            %           obj2
            %               Second CMatrix object.
            if(obj1.nColumns ~= obj2.mRows)
                msg = ['The number of columns in the first matrix '...
                       'should be equal to the number of rows in '...
                       'the second.'];
                error(msg);
            end
        end
        function res = infiniteNorm(obj)
            % Summary:  Calculates the infinite norm of the matrix, this is the
            %           maximum row sum of the magnitude of every entry.
            %
            % Returns:  res
            %               Infinite norm.
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
            % Summary:  Computes the absolute value of the matrix. By looping
            %           along each entry it returns a new matrix with each entry
            %           converted to its absolute value.
            %
            % Returns:  res
            %               Matrix with entries in absolute value.
            res = zeros(obj.mRows,obj.nColumns);
            for i=1:obj.mRows
                for j=1:obj.nColumns
                    res(i,j) = abs(obj.data(i,j));
                end
            end
            res = CMatrix(res);
        end
        function res = isUpperTriangular(obj)
            % Summary:  Checks that the matrix is upper triangular by verifying
            %           that every entry above the main diagonal is zero.
            %
            % Returns:  res
            %               Boolean to say if it is upper triangular or not.
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
            % Summary:  Checks that the matrix is lower triangular by verifying
            %           that every entry below the main diagonal is zero.
            %
            % Returns:  res
            %               Boolean to say if it is lower triangular or not.
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
            % Summary:  Checks that the matrix is symmetric by verifying that
            %           every entry at (i,j) is the same as the entry at (j,i).
            %
            % Returns:  res
            %               Boolean to say if it is symmetric or not.
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
            % Summary:  Checks that the matrix is tridiagonal by verifying that
            %           every entry not within the band is zero.
            %
            % Returns:  res
            %               Boolean to say if it is tridiagonal or not.
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
            % Summary:  Computes the transpose of the matrix. It exchanges the
            %           entries at (i,j) with the entries at (j,i).
            %
            % Returns:  res
            %               The transposed matrix.
            res = zeros(obj.nColumns,obj.mRows);
            for i=1:obj.mRows
                for j=1:obj.nColumns
                    res(j,i) = obj.data(i,j);
                end
            end
            res = CMatrix(res);
        end
        function res = getRow(obj,i)
            % Summary:  Gets the row i as a CVector object.
            %
            % Args:     i
            %               The row number to get.
            %
            % Returns:  res
            %               CVector object containing the specified row.
            res = zeros(1,obj.nColumns);
            for j=1:obj.nColumns
                res(1,j) = obj.data(i,j);
            end
            res = CVector(res);
        end
        function res = getColumn(obj,j)
            % Summary:  Gets the column j as a CVector object.
            %
            % Args:     j
            %               The row number to get.
            %
            % Returns:  res
            %               CVector object containing the specified column.
            res = zeros(obj.mRows,1);
            for i=1:obj.nColumns
                res(i,1) = obj.data(i,j);
            end
            res = CVector(res);
        end
        function res = swapSubrows(obj,i1,i2,j1,j2)
            % Summary:  Swaps the subrows of the matrix. A subrow is considered
            %           the subvector originated from taking the entries from
            %           column j1 to column j2.
            %
            % Args:     i1
            %               The first row number to swap.
            %           i2
            %               The second row number to swap.
            %           j1
            %               The initial column of the subvector.
            %           j2
            %               The final column of the subvector.
            %
            % Returns:  res
            %               CMatrix object with swapped subrows.
            res = obj.data;
            for j=j1:j2
                res(i1,j)=obj.data(i2,j);
                res(i2,j)=obj.data(i1,j);
            end
            res = CMatrix(res);
        end
        function res = getDiagonal(obj)
            % Summary:  Returns the main diagonal of the matrix. Just with one
            %           for loop it is possible to extract the entries in the
            %           diagonal and build a new matrix with them. Sparse matrix
            %           is used to speed computations in matrices with many zero
            %           entries.
            %
            % Returns:  res
            %               Diagonal matrix with the main diagonal of the
            %               matrix.
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
            % Summary:  Returns the inverse of the matrix obtained with
            %           getDiagonal().
            %
            % Returns:  res
            %               Inverse diagonal matrix.
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
            % Summary:  Returns the lower part of the matrix by looping through
            %           its lower part. It returns a new sparse matrix to speed
            %           computations.
            %
            % Returns:  res
            %               The lower part of the matrix.
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
            % Summary:  Returns the upper part of the matrix by looping through
            %           its upper part. It returns a new sparse matrix to speed
            %           computations.
            %
            % Returns:  res
            %               The upper part of the matrix.
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
            % Summary:  Finds the entry with maximum magnitude in the matrix by
            %           looping through all the elements.
            %
            % Returns:  res
            %               The maximum magnitude.
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
            % Summary:  Returns the identity matrix of size n.
            %
            % Args:     n
            %               Size of the matrix.
            % Returns:  res
            %               The identity matrix.
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
            % Summary:  Returns the permutation matrix of size n, for swapping
            %           i1 with i2 rows or columns.
            %
            % Args:     i1
            %               First index to swap.
            %           i2
            %               Second index to swap.
            %           n
            %               Size of the matrix.
            % Returns:  res
            %               The permutation matrix.
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
