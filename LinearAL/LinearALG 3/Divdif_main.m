X=[0,0.5,1,1.5,2];
Y=sin(X);
[C,D] = Newton_divdiff(X,Y)
% Given a sequence of data points(x0,y0),... ,(xn,yn)}
%the method calculates the coefficients of the interpolation polynomial of these points in the Newton form.
%Newton_divdiff(X,Y) is a code that deployed to model this method
A=hilb(2);
A_f=polyvalm(C,A)
F = funm(A,@sin)
err= norm(A_f-F)/norm(F)
