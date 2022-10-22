clc;
close all;
clear all;
A = hilb(5);
n=length(A);
D=diag(diag(A))
L=tril(A,-1)
U=triu(A,1)
M_jc=D;
M_gs=L+D;
A_jc=M_jc^(-1)*A;
A_gs=M_gs^(-1)*A;
w=1.35;
M_sor=(D+w*L)/w;
A_sor=M_sor^(-1)*A;
cond_A_sor=cond(A_sor);
x=[cond(A);cond(A_jc);cond(A_gs);cond(A_sor)]

% bar(x)