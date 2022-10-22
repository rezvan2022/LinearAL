function [U_out,S_out,V_out] = svd_decomp(A,tol)

if ~exist('tol','var')
   tol=eps*1024;
end
size_A=size(A);
itr_max=100*max(size_A);
itr_count=0;
U_out=eye(size_A(1));
S_out=A';
V_out=eye(size_A(2));

Err=realmax;
while Err>tol & itr_count<itr_max ;
    [q,S_out]=qr(S_out'); U_out=U_out*q;
    [q,S_out]=qr(S_out'); V_out=V_out*q;
    e=triu(S_out,1);
    E=norm(e(:));
    F=norm(diag(S_out));
    if F==0, F=1;end
    Err=E/F;
    itr_count=itr_count+1;
end
ss=diag(S_out);
S_out=zeros(size_A);
for n=1:length(ss)
    ssn=ss(n);
    S_out(n,n)=abs(ssn);
    if ssn<0
       U_out(:,n)=-U_out(:,n);
    end
end
if nargout<=1
   U_out=diag(S_out);
end
return