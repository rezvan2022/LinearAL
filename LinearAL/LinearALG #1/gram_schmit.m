for i=1:100
    for j=1:100
        A(i,j)=i-j;
    end
end
[m,n]=size(A)
u=zeros(m,n);
e=zeros(m,n);
for i=1:n
    v=A(:,i);
    for j=1:i-1
        u(j,i)=e(:,j)'*A(:,i);
        v=v-u(j,i)*e(:,j);
    end
    u(i,i)=norm(v);
    e(:,i)=v/u(i,i);
end
s=e(:,1)'*e(:,2)
