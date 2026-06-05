function [c,ceq]=nonlinearcon(x,m,mm)
Nnum = 2; in = 2; out = 2;
a=0.1;bbb=0.008;r=0.004;te=-10;th=50;
[A,B,b,C,d,Q,S,R,P,L,Lambda,l] = unpack(x,Nnum,in,out);
E = eye(Nnum);
K = [E; A];

Coin = [eye(in), zeros(in, 2*Nnum+1);
        zeros(1, in+Nnum*2), 1];
CoN = [zeros(Nnum,in), eye(Nnum), zeros(Nnum,Nnum), zeros(Nnum,1);
           B, zeros(Nnum,Nnum), eye(Nnum), b;zeros(1,2*Nnum+in) 1];
CoP = [zeros(out,in), C,zeros(out,Nnum), d;
           eye(in), zeros(in, 2*Nnum+1);
           zeros(1,2*Nnum+in), 1];
CoL = [zeros(Nnum,in), eye(Nnum), zeros(Nnum,Nnum), zeros(Nnum,1);
           zeros(Nnum,in), zeros(Nnum,Nnum), eye(Nnum), zeros(Nnum,1)];
 
Lamb = diag(Lambda);
LMI1 = [Q - K' * P * K, zeros(Nnum), S;
        zeros(Nnum), zeros(Nnum), zeros(Nnum);
        S', zeros(Nnum), R];

LMI1 = LMI1 + blkdiag(zeros(Nnum),P);

Piin = [-l(1) l(1)/2 0;0 -l(2) l(2)/2;l(1)/2 l(2)/2 0];
PiL = [Q, S; S', R];
PiL = CoL' * PiL * CoL;

% PiN 构造
PiN = zeros(2 * Nnum, 2 * Nnum);
for i = 1:Nnum
    for j = i:Nnum
        ei = E(:, i);
        ej = E(:, j);
        if i == j
            e = ei;
        else
            e = ei - ej;
        end
        M = [-2*(e*e'),  e*e';
              e*e',     zeros(Nnum)];
        PiN = PiN + L(i,j) * M;
    end
end
PiN = [PiN zeros(2*Nnum,1);zeros(1,2*Nnum+1)]+ [-2*Lamb zeros(Nnum,Nnum) zeros(Nnum,1);zeros(Nnum,2*Nnum+1);zeros(1,Nnum) zeros(1,Nnum) 2*ones(1,Nnum)*Lamb*ones(Nnum,1)];
PiN =  CoN' * PiN * CoN;
Pi0=[0 2*a r*th 0 bbb*te;2*a 0 0 r*th bbb*te;r*th 0 0 0 0;0 r*th 0 0 0;bbb*te bbb*te 0 0 0];
Coin0=[mm 0 0 0 m(1);0 mm 0 0 m(2);0 0 1 0 0;0 0 0 1 0;0 0 0 0 1];
LMI2 = PiL + PiN +Coin'*Piin*Coin-CoP' * Coin0' * Pi0 * Coin0 * CoP;

c=[-min(eig(LMI1));max(eig(LMI2));-min(eig(P))];
ceq=[];
end