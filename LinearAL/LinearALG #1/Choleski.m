clc
clear all
close all;
%In linear algebra, the Cholesky decomposition or Cholesky factorization
%is a decomposition of a Hermitian, positive-definite matrix into the product of 
%a lower triangular matrix and its conjugate transpose, which is useful for efficient numerical solutions, e.g., Monte Carlo simulations.
n=10;
M=randn(n,n)
A=M*M';
[n,n] = size( A );
L = zeros( n, n );
for i=1:n
    L(i, i) = sqrt( A(i, i) - L(i, :)*L(i, :)' );
 
    for j=(i + 1):n
        L(j, i) = ( A(j, i) - L(i, :)*L(j, :)' )/L(i, i);
    end
end
L
