%pade_approximation
syms x
f1= exp(x)
%f1=sin(x)
pade_f1=pade(f1,'Order',[2 2])
n=3;
%  A=generateSPDmatrix(n);
A=hilb(n);
% pade_f1_matrix=double(subs(pade_f1,x,A));
% f1o=expm(A);
% err1=norm(expm(A)-pade_f1_matrix)/norm(f1o)
[num,den] = numden(pade_f1)
num=coeffs(num)
den=coeffs(den)
N_p=num(1)*eye(n)+num(2)*A+num(3)*A^2;
D_p=den(1)*eye(n)+den(2)*A+den(3)*A^2;
pade_app_f1=double(N_p./D_p)
% num=polyvalm(num,A)
% den=polyvalm(den,A);
% pade_f1_matrix=n./d;
expm(A)


 