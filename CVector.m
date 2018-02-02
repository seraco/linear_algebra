classdef CVector
    % Own implementation for vectors.
    %   The class is helpful for defining custom algorithms to handle
    %   calculations with vectors.
    
    % Version:  0.0.1
    % Author:   S.Ramon
    
    properties
        nElements
        data
    end
    
    methods
        function obj = CVector(vector)
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
            areEqualSize(obj1,obj2);
            res = zeros(obj1.nElements,1);
            for i=1:obj1.nElements
                res(i,1) = obj1.data(i,1)+obj2.data(i,1);
            end
            res = CVector(res);
        end
        function res = minus(obj1,obj2)
            areEqualSize(obj1,obj2);
            res = zeros(obj1.nElements,1);
            for i=1:obj1.nElements
                res(i,1) = obj1.data(i,1)-obj2.data(i,1);
            end
            res = CVector(res);
        end
        function res = times(obj1,obj2)
            areEqualSize(obj1,obj2);
            res = 0.0;
            for i=1:obj1.nElements
                res = res+obj1.data(i,1)*obj2.data(i,1);
            end
        end
        function res = mtimes(obj1,obj2)
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
            if(obj1.nElements ~= obj2.nElements)
                msg = ['The number of elements should be the '...
                       'same in both vectors.'];
                error(msg);
            end
        end
        function res = getSubvector(obj,i1,i2)
            res = zeros(i2-i1,1);
            for i=i1:i2
                res(i+1-i1) = obj.data(i,1);
            end
            res = CVector(res);
        end
        function res = maxElementIndex(obj)
            res = 1;
            for i=1:obj.nElements
                if(obj.data(i) > obj.data(res))
                    res = i;
                end
            end
        end
    end
    
end

