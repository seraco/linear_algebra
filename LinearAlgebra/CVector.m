classdef CVector
    % Class:    CVector
    %
    % Summary:  Own implementation for vectors. The class is helpful for
    %           defining custom algorithms to handle calculations with vectors.
    %
    % Methods:  CVector
    %               Constructor of the class.
    %           plus
    %               Overloading of + operator.
    %           minus
    %               Overloading of - operator.
    %           times
    %               Overloading of .* operator with inner product.
    %           mtimes
    %               Overloading of * operator with outer product.
    %           areEqualSize
    %               Checks that two vectors are of equal size.
    %           getSubvector
    %               Gets a subvector from the whole vector.
    %           maxElementIndex
    %               Returns the index of the maximum magnitude element.
    %           maxElementMagnitude
    %               Returns the maximum magnitude element.

    % Version:  0.0.1
    % Author:   S.Ramon
    % License:  MIT

    properties
        nElements   % number of elements
        data        % stores the vector data
    end

    methods
        function obj = CVector(vector)
            % Summary:  Constructor. Initialises the properties of the object
            %           and returns the vector as a column vector.
            %
            % Args:     vector
            %               Vector that is stored as data.
            %
            % Returns:  obj
            %               CVector object with initilized properties.
            vectorSize = size(vector);
            if(vectorSize(1) > vectorSize(2))
                obj.nElements = vectorSize(1);
                data = zeros(obj.nElements,1);
                for i=1:obj.nElements
                    data(i,1) = vector(i,1);
                end
            else
                obj.nElements = vectorSize(2);
                data = zeros(obj.nElements,1);
                for i=1:obj.nElements
                    data(i,1) = vector(1,i);
                end
            end
            obj.data = data;
        end
        function res = plus(obj1,obj2)
            % Summary:  Overloading of + operator.
            %
            % Args:     obj1
            %               CVector object on the left of operator +.
            %           obj2
            %               CVector object on the right of operator +.
            %
            % Returns:  res
            %               CVector object with summed data entries.
            areEqualSize(obj1,obj2);
            res = zeros(obj1.nElements,1);
            for i=1:obj1.nElements
                res(i,1) = obj1.data(i,1)+obj2.data(i,1);
            end
            res = CVector(res);
        end
        function res = minus(obj1,obj2)
            % Summary:  Overloading of - operator.
            %
            % Args:     obj1
            %               CVector object on the left of operator -.
            %           obj2
            %               CVector object on the right of operator -.
            %
            % Returns:  res
            %               CVector object with substracted data entries.
            areEqualSize(obj1,obj2);
            res = zeros(obj1.nElements,1);
            for i=1:obj1.nElements
                res(i,1) = obj1.data(i,1)-obj2.data(i,1);
            end
            res = CVector(res);
        end
        function res = times(obj1,obj2)
            % Summary:  Overloading of .* operator to implement the inner
            %           product.
            %
            % Args:     obj1
            %               CVector object on the left of operator .*.
            %           obj2
            %               CVector object on the right of operator .*.
            %
            % Returns:  res
            %               Inner product of obj1 with obj2.
            areEqualSize(obj1,obj2);
            res = 0.0;
            for i=1:obj1.nElements
                res = res+obj1.data(i,1)*obj2.data(i,1);
            end
        end
        function res = mtimes(obj1,obj2)
            % Summary:  Overloading of * operator to implement the outer
            %           product.
            %
            % Args:     obj1
            %               CVector object on the left of operator *.
            %           obj2
            %               CVector object on the right of operator *.
            %
            % Returns:  res
            %               CMatrix object with the outer product of obj1 with
            %               obj2.
            areEqualSize(obj1,obj2);
            res = zeros(obj1.nElements,obj1.nElements);
            for i=1:obj1.nElements
                for j=1:obj1.nElements
                    res(i,j) = obj1.data(i,1)*obj2.data(j,1);
                end
            end
            res = CMatrix(res);
        end
        function areEqualSize(obj1,obj2)
            % Summary:  Checks that the data of two CVector objects is of equal
            %           size. If the number of elements is different, displays
            %           an error message.
            %
            % Args:     obj1
            %               First CVector object.
            %           obj2
            %               Second CVector object.
            if(obj1.nElements ~= obj2.nElements)
                msg = ['The number of elements should be the '...
                       'same in both vectors.'];
                error(msg);
            end
        end
        function res = getSubvector(obj,i1,i2)
            % Summary:  Returns a subvector of a vector, from the element i1 to
            %           the element i2.
            %
            % Args:     i1
            %               Initial element index of the subvector.
            %           i2
            %               Final element index of the subvector.
            %
            % Returns:  res
            %               CVector object with the data of the subvector.
            if(i2 < i1)
                res = CVector([]);
                return
            end
            res = zeros(i2-i1,1);
            for i=i1:i2
                res(i+1-i1) = obj.data(i,1);
            end
            res = CVector(res);
        end
        function res = maxElementIndex(obj)
            % Summary:  Returns the index of the maximum magnitude element by
            %           looping through all of the elements of the vector.
            %
            % Returns:  res
            %               The index of the maximum magnitude element.
            res = 1;
            for i=1:obj.nElements
                absElement = abs(obj.data(i));
                if(absElement > obj.data(res))
                    res = i;
                end
            end
        end
        function res = maxElementMagnitude(obj)
            % Summary:  Returns the maximum magnitude element by looping through
            %           all of the elements of the vector.
            %
            % Returns:  res
            %               The maximum magnitude element.
            res = 1;
            for i=1:obj.nElements
                absElement = abs(obj.data(i));
                if(absElement > obj.data(res))
                    res = i;
                end
            end
            res = obj.data(i);
        end
    end

end
