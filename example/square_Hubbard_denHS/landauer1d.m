function[gg]=landauer1d(L,eta,cdw,W,en)

h=diag(ones(1,L-1),1)+diag(ones(1,L-1),-1);
h=h+diag((rand(1,L)-0.5))*2*W;
h=h+diag((-1).^(1:L))*cdw;

Gr=inv(en-h+1i*eye(L)*eta);
Ga=inv(en-h-1i*eye(L)*eta);

gg=Gr(1,L)*Ga(L,1);

end

