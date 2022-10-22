clc
clear all
close all
tic
n=10;
A=zeros(n,n)
for i=1:n
    for j=1:n
        if i==j
            A(i,j)=10;
        
        end
    end
end
for i=2:n
    A(i-1,i)=-1;
end
for i=2:n
    A(i,i-1)=-1;
end
for i=3:n
    A(i,i-2)=-1;
end
for i=3:n
    A(i-2,i)=-1;
end
it_nu=1000;
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
w = 2./(1+sqrt(1-rho_Tj^2));	
disp('w =');disp(w);
Tw = inv(D-w*L)*((1-w)*D+w*U);		
cw = w*inv(D-w*L)*b;			
err_sor = 1e-05;
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
   