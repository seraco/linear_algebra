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
            obj.b = CMatrix(b);
        end
        function res = backwardSubstitution(obj)
            if(~obj.A.isUpperTriangular)
                msg = ['For backward substitution the matrix of '...
                       'coefficients should be upper triangular'];
                error(msg);
            end
            res = zeros(size(obj.b.data));
            for i=obj.b.mRows:-1:1
                for j=obj.b.mRows:-1:i+1
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
            for i=1:obj.b.mRows
                for j=1:i-1
                    res(i)=res(i)+obj.A.data(i,j)*res(j);
                end
                res(i)=(obj.b.data(i)-res(i))/obj.A.data(i,i);
            end
        end
    end
    
end

