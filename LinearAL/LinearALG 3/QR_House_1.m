clc
clear all
close all
%QR HouseHolder
n=5;
A=hilb(n);
% A=[1 0 2;2 1 4;4 6 4];
    D=diag(sort(eig(A),'descend'));
% while norm(D-A)>0.01
for i=1:10
A=[1 0 2;
    2 1 4;
    4 6 4];
[q,r]=qr_h_main(A);
A=r*q;
i=i+1;
end
