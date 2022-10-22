tic
clc
clear all
close all
%minimum residual methods
n=10000;
% A=hilb(n);
 A=generateSPDmatrix(n);
 A=sparse(A);
for i=1:n
    b(i,1)=1;
end
x0=A\b;
it_nu=1000;
tol=10^-3;
   iter = 0;                                         % initialization
   flag = 0;
   res_norm = zeros(it_nu+1,1);
   b_norm = norm( b );
   if  ( b_norm == 0.0 ), b_norm = 1.0; end
   tol = tol*b_norm;
   x = x0;
   r = ( b-A*x );
   res_norm(1) = norm( r );
   if ( res_norm(1) < tol ) return, end
   [n,n] = size(A);                                  % initialize workspace
   V(1:n,1:it_nu+1) = zeros(n,it_nu+1);
   H(1:it_nu+1,1:it_nu) = zeros(it_nu+1,it_nu);
   e1    = zeros(it_nu+1,1);
   e1(1) = 1.0;
   s = norm( r )*e1;
   V(:,1) = r / norm( r );
   for i = 1:it_nu                              % begin iteration

       w = A*V(:,i);                                 % basis using Gram-Schmidt
       for k = 1:i
           H(k,i)= w'*V(:,k);
           w = w - H(k,i)*V(:,k);
       end
       H(i+1,i) = norm( w );
       V(:,i+1) = w / H(i+1,i);
       
       y = H(1:i+1,1:i)\s(1:i+1);                    % solve projected system 
       res_norm(i+1) = norm( H(1:i+1,1:i)*y - s(1:i+1) );
       if ( res_norm(i+1) <= tol )                  % update approximation
          res_norm = res_norm(1:i+1);
	  break;
       end
   end
   x = x + V(:,1:i)*y;
   toc