clc
clear all
close all
tic
n=40;
B=full(gallery('tridiag',n,-1,4,-1));
I=eye(n);
A=full(blktridiag(B,-I,-I,n));
[m,m]=size(A)
b=zeros(m,1);
for i=1:m
    b(i,1)=1;
end

n = size(A,1);
D = diag(diag(A));
L = tril(-A,-1);
U = triu(-A,1);
it_nu=10000;
Tg = inv(D-L)*U; 
cg = inv(D-L)*b;

err_gs = 1e-03;
k = 1;
x = zeros(n,1);			

while k <= it_nu
   x(:,k+1) = Tg*x(:,k) + cg;
   if norm(x(:,k+1)-x(:,k)) < err_gs
        disp('Convergence has been achived providing these results')
      
      disp('Interation enough to converge= ');disp(k); disp('x = ');disp(x(:,k+1));
      break
   end
   k = k+1;
end

if norm(x(:,k+1)- x(:,k)) > err_gs || k > it_nu
   disp('No. of iteration is not enough')
    disp('Results for current No. of iteration:')
   disp(err_jc);
   disp(x');
end
toc