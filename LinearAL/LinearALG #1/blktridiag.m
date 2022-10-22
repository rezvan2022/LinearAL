function A = blktridiag(Amd,Asub,Asup,n)
[p,q] = size(Amd);
if (p~=q) || (p<1)
  error 'Blocks must be (non-empty) square arrays or scalars'
end
if any(p~=size(Asub)) || any(p~=size(Asup)) 
  error 'Amd, Asub, Asup are not identical in size'
end

if isempty(n) || (length(n)>1) || (n<1) || (n~=floor(n))
  error 'n must be a positive scalar integer'
end

[ind1,ind2,ind3]=ndgrid(0:p-1,0:p-1,0:n-1);
rind = 1+ind1(:)+p*ind3(:);
cind = 1+ind2(:)+p*ind3(:);
v = repmat(Amd(:),n,1);
if n>1
  % sub-diagonal
  [ind1,ind2,ind3]=ndgrid(0:p-1,0:p-1,0:n-2);
  rind = [rind;1+p+ind1(:)+p*ind3(:)];
  cind = [cind;1+ind2(:)+p*ind3(:)];
  v=[v;repmat(Asub(:),n-1,1)];
  
  [ind1,ind2,ind3]=ndgrid(0:p-1,0:p-1,0:n-2);
  rind = [rind;1+ind1(:)+p*ind3(:)];
  cind = [cind;1+p+ind2(:)+p*ind3(:)];
  v=[v;repmat(Asup(:),n-1,1)];
end
A = sparse(rind,cind,v,n*p,n*p);