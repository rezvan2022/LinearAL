clc
clear all
close all
tic
n=50;
A=generateSPDmatrix(n);
[m,m]=size(A);
b=zeros(m,1);
for i=1:m
    b(i,1)=1;
end
n = size(A,1);
D = diag(diag(A));
L = tril(-A,-1);
U = triu(-A,1);
it_nu=10000;
I=eye(n);
w=0.1;
Tj=-(D)^(-1)*(L+U);
T_wj=w*Tj-((1-w)*I);
c_wj=w*(D)^(-1)*b;
Ln_wj=eig(T_wj);
r_wj= max(abs (Ln_wj));

% Tj = inv(D)*(L+U);
% cj = inv(D)*b;

err_jc = 1e-05;
k = 1;
x = zeros(n,1);			

while k <= it_nu
   x(:,k+1) = T_wj*x(:,k) + c_wj;
   if norm(x(:,k+1)-x(:,k)) < err_jc
       disp('Convergence has been achived providing these results')
      
      disp('Interation enough to converge= ');disp(k); disp('x = ');disp(x(:,k+1));
      break
   end
   k = k+1;
end

if norm(x(:,k+1)- x(:,k)) > err_jc || k > it_nu
   disp('No. of iteration is not enough')
    disp('Results for current No. of iteration:')
   disp(err_jc);
   disp(x');
end
toc
