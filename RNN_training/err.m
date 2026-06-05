function e=err(x,u,yreal)
    Nnum = 2; in = 2; out = 2;
    [A,B,b,C,d,~,~,~,~,~,~] = unpack(x,Nnum,in,out);
    y=RNN(u,A,B,b,C,d,Nnum,out);
    e=norm(y-yreal');
end