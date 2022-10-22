clc
clear all
close all
mm=2;
nn=3;
A=rand(mm,nn);
tol=10^-3;
[U_out,S_out,V_out] = svd_decomp(A,tol);
A_est=U_out*S_out*V_out';
err=norm(A_est-A)/norm(A)
