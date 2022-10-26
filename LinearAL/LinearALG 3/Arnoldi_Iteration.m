clc
clear all
close all
n=10;
%The Arnoldi iteration uses the modified Gram–Schmidt process to produce a sequence of orthonormal vectors, 
%q1, q2, q3, …, called the Arnoldi vectors, such that for every n, the vectors q1, …, qn span the Krylov subspace
%K_n
%The idea of the Arnoldi iteration as an eigenvalue algorithm is to compute the eigenvalues in the Krylov subspace
%  A=hilb(n);
 A=generateSPDmatrix(n);
for i=1:n
    b(i,1)=1;
end
res=(A\b-b);
v1=res/norm(res);
%v1=zeros(n)
m=5;
% Run m steps of basic Arnoldi iteration
n = length(A);
V = zeros(n,m+1);
H = zeros(m+1,m);
V(:,1) = v1;
for j = 1:m
r = A*V(:,j);
for k = 1:j
H(k,j) = V(:,k)'*r;
r = r-V(:,k)*H(k,j);
end
H(j+1,j) = norm(r);
V(:,j+1) = r/H(j+1,j);
end
