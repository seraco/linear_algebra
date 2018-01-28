function evec = upper_e_vector(m,n)
% This function constructs the e evector for the given Gauss elimination step
% m - the step number of Gauss elimination method
% n - the size of the matrix A
for i=1:n %lines 5-7: initialisation of the e-vector
evec(i) = 0.0;
end
evec(m)=1.0; %since e-vector denotes to the m_th column of identity matrix, make the m_th entry in e-vector as 1.0
