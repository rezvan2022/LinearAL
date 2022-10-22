clc
clear all
close all
n=3;
A=randn(n,n);
m=length(A);
L=zeros(size(A));
U=zeros(size(A));
L(:,1)=A(:,1);
U(1,:)=A(1,:)/L(1,1);
U(1,1)=1;
for k=2:m
for j=2:m
    for i=j:m
        L(i,j)=A(i,j)-dot(L(i,1:j-1),U(1:j-1,j));
    end
    U(k,j)=(A(k,j)-dot(L(k,1:k-1),U(1:k-1,j)))/L(k,k);
end
end
