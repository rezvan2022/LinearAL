clc
clear all
close all
n=5;
A=randn(n,n);
  L=zeros(n,n);
  U=zeros(n,n);
  for i=1:n
  % Finding L
  for k=1:i-1
  L(i,k)=A(i,k);
  for j=1:k-1
  L(i,k)= L(i,k)-L(i,j)*U(j,k);
  end
  L(i,k) = L(i,k)/U(k,k);
  end
  % Finding U
  for k=i:n
  U(i,k) = A(i,k);
  for j=1:i-1
  U(i,k)= U(i,k)-L(i,j)*U(j,k);
  end
  end
  end
  for i=1:n
  L(i,i)=1;
  end
b=zeros(n,1);
for i=1:n
    b(i,1)=1;
end
y=linsolve(L,b);
x=linsolve(U,y)
x_d=linsolve(A,b)


