function A = generateSPDmatrix(n)
% genearte a matrix which is symmetric positive definite
A = rand(n,n); 
A = 0.5*(A+A'); 
A = A*A';
A = A + n*eye(n);
end
