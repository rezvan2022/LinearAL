clc
clear all
close all
tic
n=10;
A1=randn(n,n);
A2=A1';
A=A1*A2;
b=zeros(n,1);
for i=1:n
    b(i,1)=1;
end
it_nu=1000000;
n = size(A,1);
D = diag(diag(A));
L = tril(-A,-1);
U = triu(-A,1);
b=zeros(n,1);
for i=1:n
    b(i,1)=1;
end
Tj = inv(D)*(L+U);				
rho_Tj = max(abs(eig(Tj)));     
w=0.1;
disp('w =');disp(w);
B1=inv(D+w*U)*(-w*L+(1-w)*D);% BACKWARD
B2=inv(D+w*L)*(-w*U+(1-w)*L);%Forward
cw=w*(2-w)*inv(D+w*U)*D*inv(D+w*L)*b;
size(cw)
Tw=B1*B2;
size(Tw)
err_sor = 1e-02;
k = 1;
x = zeros(n,1);                
while k <= it_nu
   x(:,k+1) = Tw*x(:,k) + cw;
   if norm(x(:,k+1)-x(:,k)) < err_sor
       disp('Convergence has been achived providing these results')
      
      disp('Interation enough to converge= ');disp(k); disp('x = ');disp(x(:,k+1));
      break
   end
   k = k+1;
end

if norm(x(:,k+1)- x(:,k)) > err_sor || k > it_nu
    disp('No. of iteration is not enough')
    disp('Results for current No. of iteration:')
   disp(err_jc);
   disp(x');
end
   