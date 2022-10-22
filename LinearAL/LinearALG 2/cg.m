clc
clear all
close all
tic
n=10000;
% A = generateSPDmatrix(n);
% A=A*A';
err_matrix=eps*ones(n,n);
A=hilb(n)+err_matrix;
cond(A)
b=zeros(n,1);
for i=1:n
            b(i,1)=1;
end
max_it=1000;
in_gs=b;% initail guass
err_cg=10^-3;
x = in_gs;         
r = b - A*in_gs;
p = r;         
rho = r'*r;
it_nu = 0;     % Init counter for number of iterations
normb = norm(b);
while (norm(r)/normb > err_cg)   % Test break condition
	a = A*p;
	alpha = rho/(a'*p);
	x = x + alpha*p;
	r = r - alpha*a;
	rho_new = r'*r;
	p = r + rho_new/rho * p;
	rho = rho_new;
	it_nu = it_nu + 1;
	if (it_nu == max_it)         
		disp( 'Need larger Iteration Number');                  
		break
        else
    end
    

end
        disp('Final answer is:')
        disp(x)
        disp('Required Iteration Number is:')
        disp(it_nu)
toc