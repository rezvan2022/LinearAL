clc
clear all
close all
%iterative matrix inversion,
n=3;
A=randn(n,n);
E=inf; 
err=10^-5;
landa=max(eig(A*A'));
c=(2/landa);  
alpha=c-eps
x=alpha*A';
while E>err
 x=x*(2*eye(n)-A*x); 
 E=norm(A*x-eye(n)); 
end 
Inv_Mat=x
inv(A)
  