clc
clear all
close all;
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
