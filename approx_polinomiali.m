clear all
close all
% f=@(x) 1./(x.^2+1);
f=@(x) 5+x.^2+x.^3;
f=@(x) sin(x)+log(x.^8);

N=50;
noiselevel=0;
maxdeg=15;

x0=rand(N,1);
X=2.*x0-1;


noise=noiselevel.*randn(size(X));

Y=(f(X)+noise);
W=zeros(maxdeg,maxdeg); 
W2=W; 
ErrLs=zeros(maxdeg,1);

for i=1:maxdeg
    A=RegMat(X,i);
    B=pinv(A);
    W(end:-1:end-i+1,i)=B*Y;
    ErrLs(i)=norm(A*B*Y-Y,2);
end

M=100;

% tt0=rand(M,1); tt=2*tt0-1; tt=sort(tt);

t=linspace(-1,1,M);
R=zeros(maxdeg,M); R2=0*R;

for i=1:maxdeg
%     R2(i,:)= max(abs(polyval(W(:,i),tt)-f(tt)));
    R(i,:)=max(abs(polyval(W(:,i),t)-f(t)));
end

semilogy(1:maxdeg,ErrLs',"-o")
hold on
semilogy(1:maxdeg,R',"r-o")
legend("Errore LS","Errore residuo")

figure()
hold on
fplot(f, [-1,1])
for i=1:maxdeg
    plot(t,polyval(W(:,i),t),"*")
end
legend()

function A=RegMat(X,m)

    [n,r]=size(X);
    A=zeros(n,m);
    
    for i=0:m-1
        A(:,i+1)=X.^i;
    end

end