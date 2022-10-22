X=[0,0.5,1,1.5,2];
Y=sin(X);
[C,D] = Newton_divdiff(X,Y);
A=hilb(2);
A_f=polyvalm(C,A)
F = funm(A,@sin)
err= norm(A_f-F)/norm(F)