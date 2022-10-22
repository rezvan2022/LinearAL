clc
clear all
close all
tic
n=10;
h=1/(n+1);
ta=0.1;
I_n=eye(n^2);
V_n=(1/h^2)*full(gallery('tridiag',n,-1,2,-1));
K=kron(eye(n),V_n)+kron(V_n,eye(n));
A=K+((3-sqrt(3))/ta)*I_n;
for i=1:n^2
    for j=1:n^2
        B(i,j)=j*(1-i)/(ta*(j+1)^2);
    end
end
m = size(A,1);
D = diag(diag(A));
L = tril(-A,-1);
U = triu(-A,1);
it_nu=2000;
Tg = inv(D-L)*U; 
cg = inv(D-L)*B;

err_gs = 1e-05;
k = 1;
x = zeros(m,m);			

while k <= it_nu
   x(:,:,k+1) = Tg*x(:,:,k) + cg;
   if norm(x(:,:,k+1)-x(:,:,k)) < err_gs
        disp('Convergence has been achived providing these results')
      
      disp('Interation enough to converge= ');disp(k); disp('x = ');disp(x(:,:,k+1));
      break
   end
   k = k+1;
end
if norm(x(:,:,k+1)- x(:,:,k)) > err_gs || k > it_nu
   disp('No. of iteration is not enough')
   disp('Results for current No. of iteration:')
   disp(err_jc);
   disp(x');
end
err=norm(A*x(:,:,end)-B)/norm(B)
% surf(x(:,:,end))
toc;
