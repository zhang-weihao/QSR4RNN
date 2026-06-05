function y=RNN(u,A,B,b,C,d,Nnum,out)
    [~,k]=size(u);
    y=zeros(out,k);
    x=zeros(Nnum,k);
    x(:,1)=tanh(B*u(:,1)+b);
    y(:,1)=C*x(:,1)+d;
    for i=2:k
        x(:,i)=tanh(A*x(:,i-1)+B*u(:,i)+b);
        y(:,i)=C*x(:,i)+d;
    end
end