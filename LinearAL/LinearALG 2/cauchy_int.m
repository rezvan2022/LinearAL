clc
clear all
close all
%Cauchy Inegral
n=3;
A=full(gallery('tridiag',n,1,4,1));
m=size(A);
fun=@(z)sinh(z)*inv(z*eye(m)-A);
r=15;
sinh_cauchy=quadv(@(t)fun(exp(1i*t)*r)*1i*r*exp(1i*t),0,2*pi)/(2i*pi);
sinh_mat=(expm(A)-expm(-A))*0.5;
err=norm(sinh_cauchy-sinh_mat)/norm(sinh_mat)