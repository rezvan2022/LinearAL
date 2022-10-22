clc
clear all
close all
m=100; n=20;
A=randn(m,n)
e = zeros(m,n);
u = zeros(n);
 
for i=1:n
    u(i,i) = norm(A(:,i));
    e(:,i) = A(:,i)/u(i,i);
    u(i,i+1:n) = e(:,i)'*A(:,i+1:n);
    A(:,i+1:n) = A(:,i+1:n) - e(:,i)*u(i,i+1:n);
end
