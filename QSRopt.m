function [value, time] =  QSRopt(An,B1,Bn,bn,C,d,Piin,method)
    tic;
    sizeB1=size(B1);
    Nnum=sizeB1(1);
    in=sizeB1(2);
    sizebn=size(bn);
    Lnum=sizebn(3);
    sized=size(d);
    out=sized(1);
    n=Nnum*Lnum;
    if method == 1
        A = An(:,:,1); B = B1; b=bn(:,:,1);

for i = 2 : Lnum
    A = blkdiag(A,An(:,:,i));
    B = blkdiag(B,Bn(:,:,i-1));
    b = vertcat(b,bn(:,:,i));
end
B = [B, zeros(n, Nnum)];

E = eye(n);
K = [E; A];


Coin = [eye(in), zeros(in, 2*n+1);
        zeros(1, in+n*2), 1];

CoN = [zeros(n,in), eye(n), zeros(n,n), zeros(n,1);
       B, eye(n), b];

CoL = [zeros(n,in), eye(n), zeros(n,n), zeros(n,1);
       zeros(n,in), zeros(n,n), eye(n), zeros(n,1)];

En = [zeros(Nnum,Nnum*(Lnum-1)), eye(Nnum)];
CoP = [zeros(out, in), C*En, zeros(out, n), d;
       eye(in), zeros(in, 2*n+1);
       zeros(1,2*n+in), 1];

cvx_begin sdp
    cvx_quiet(true)
    cvx_solver mosek

    variable Q(n,n) symmetric
    variable S(n,n)
    variable R(n,n) symmetric
    variable P(2*n,2*n) semidefinite
    variable Ay(out,out) semidefinite
    variable by(out,1)
    variable L(n,n) symmetric nonnegative
    variable l1(1,1) nonnegative

    LMI1 = [Q - K' * P * K, zeros(n), S;
            zeros(n), zeros(n), zeros(n);
            S', zeros(n), R];
    
    LMI1 = LMI1 + blkdiag(zeros(n),P);

    PiL = [Q, S; S', R];
    PiL = CoL' * PiL * CoL;

    % PiN 构造
    PiN = zeros(2 * n, 2 * n);
    for i = 1:n
        for j = i:n
            ei = E(:, i);
            ej = E(:, j);
            if i == j
                e = ei;
            else
                e = ei - ej;
            end
            M = [-2*(e*e'),  e*e';
                  e*e',     zeros(n)];
            PiN = PiN + L(i,j) * M;
        end
    end
    PiN =  CoN' * PiN * CoN;
    temp1 = l1 * Coin' * Piin * Coin;
    % temp2 = CoP' * diag([0 0 1]) * CoP;
    temp2 = zeros(2*n+in+1,2*n+in+1);
    temp2(2*n+in+1,2*n+in+1) = 1;
    Aby = [Ay zeros(out,in) by];
    LMI2 = [PiL + PiN + temp1 - temp2, CoP' * Aby';
            Aby * CoP, -eye(out)];
    
    maximize(log_det(Ay))
    subject to
        LMI1 >= 0
        LMI2 <= 0
        P >= 0
        for i=1:n
            for j=1:n
                L(i,j) >= 1e-8
            end
        end
cvx_end
    end
    if method == 2
        E = eye(Nnum);
LMI1=cell(Lnum,1);
LMI2=cell(Lnum,1);
Coin = [zeros(1,Nnum*(Lnum-1)+in+out-1) 1 0;zeros(out,Nnum*(Lnum-1)+in+out) 1];

cvx_begin sdp
    cvx_solver mosek
    cvx_quiet(true)
    variable blk(in+out+n-Nnum+1, in+out+n-Nnum+1) symmetric
    variable Ay(1,1) semidefinite
    variable by(1,1)
    variable l1(1,1) nonnegative
    for i = 1:Lnum
        eval(sprintf('variable L%d(%d,%d) symmetric nonnegative', i, Nnum,Nnum));
        eval(sprintf('variable Q%d(%d,%d) symmetric', i, Nnum,Nnum));
        eval(sprintf('variable S%d(%d,%d)', i, Nnum,Nnum));
        eval(sprintf('variable R%d(%d,%d) symmetric', i, Nnum,Nnum));
        eval(sprintf('variable P%d(%d,%d) semidefinite', i, 2*Nnum,2*Nnum));
        if i==1
            eval(sprintf('variable Pi%d(%d,%d) symmetric', i, in+Nnum+1, in+Nnum+1));
        elseif i==Lnum
            eval(sprintf('variable Pi%d(%d,%d) symmetric', i, out+Nnum+1, out+Nnum+1));
        else
            eval(sprintf('variable Pi%d(%d,%d) symmetric', i, 2*Nnum+1, 2*Nnum+1));
        end
    end
    Pi = cell(1,Lnum);
    L = cell(1,Lnum);
    Q = cell(1,Lnum);
    S = cell(1,Lnum);
    R = cell(1,Lnum);
    P = cell(1,Lnum);
    P11=cell(Lnum,1);
    P12=cell(Lnum,1);
    P13=cell(Lnum,1);
    P22=cell(Lnum,1);
    P23=cell(Lnum,1);
    P33=cell(Lnum,1);
    for i = 1:Lnum
        Pi{i} = eval(sprintf('Pi%d', i));
        L{i} = eval(sprintf('L%d', i));
        Q{i} = eval(sprintf('Q%d', i));
        S{i} = eval(sprintf('S%d', i));
        R{i} = eval(sprintf('R%d', i));
        P{i} = eval(sprintf('P%d', i));
    end
    for i = 1 : Lnum
        if i == 1
            P11{i}=Pi{i}(1:Nnum,1:Nnum);
            P12{i}=Pi{i}(1:Nnum,Nnum+1:Nnum+in);
            P13{i}=Pi{i}(1:Nnum,Nnum+in+1);
            P22{i}=Pi{i}(Nnum+1:Nnum+in,Nnum+1:Nnum+in);
            P23{i}=Pi{i}(Nnum+1:Nnum+in,Nnum+in+1);
            P33{i}=Pi{i}(Nnum+in+1,Nnum+in+1);
        elseif i==Lnum
            P11{i}=Pi{i}(1:out,1:out);
            P12{i}=Pi{i}(1:out,out+1:out+Nnum);
            P13{i}=Pi{i}(1:out,out+Nnum+1);
            P22{i}=Pi{i}(out+1:out+Nnum,out+1:out+Nnum);
            P23{i}=Pi{i}(out+1:out+Nnum,out+Nnum+1);
            P33{i}=Pi{i}(out+Nnum+1,out+Nnum+1);
        else
            P11{i}=Pi{i}(1:Nnum,1:Nnum);
            P12{i}=Pi{i}(1:Nnum,Nnum+1:Nnum+Nnum);
            P13{i}=Pi{i}(1:Nnum,Nnum+Nnum+1);
            P22{i}=Pi{i}(Nnum+1:Nnum+Nnum,Nnum+1:Nnum+Nnum);
            P23{i}=Pi{i}(Nnum+1:Nnum+Nnum,Nnum+Nnum+1);
            P33{i}=Pi{i}(Nnum+Nnum+1,Nnum+Nnum+1);
        end   
    end
    for i = 1 : Lnum
        if i == 1
            CoN = [zeros(Nnum,in), eye(Nnum), zeros(Nnum,Nnum), zeros(Nnum,1);
           B1, zeros(Nnum,Nnum), eye(Nnum), bn(:,:,i)];
            CoP = [zeros(Nnum,in), eye(Nnum),zeros(Nnum,Nnum), zeros(Nnum,1);
           eye(in), zeros(in, 2*Nnum+1);
           zeros(1,2*Nnum+1), 1];
            CoL = [zeros(Nnum,in), eye(Nnum), zeros(Nnum,Nnum), zeros(Nnum,1);
           zeros(Nnum,in), zeros(Nnum,Nnum), eye(Nnum), zeros(Nnum,1)];
        elseif i==Lnum
            CoN = [zeros(Nnum,Nnum), eye(Nnum), zeros(Nnum,Nnum), zeros(Nnum,1);
           Bn(:,:,i-1), zeros(Nnum,Nnum), eye(Nnum), bn(:,:,i)];
            CoP = [zeros(out,Nnum), C,zeros(out,Nnum), d;
           eye(Nnum), zeros(Nnum, 2*Nnum+1);
           zeros(1,3*Nnum), 1];
            CoL = [zeros(Nnum,Nnum), eye(Nnum), zeros(Nnum,Nnum), zeros(Nnum,1);
           zeros(Nnum,Nnum), zeros(Nnum,Nnum), eye(Nnum), zeros(Nnum,1)];
        else
            CoN = [zeros(Nnum,Nnum), eye(Nnum), zeros(Nnum,Nnum), zeros(Nnum,1);
           Bn(:,:,i-1), zeros(Nnum,Nnum), eye(Nnum), bn(:,:,i)];
            CoP = [zeros(Nnum,Nnum), eye(Nnum),zeros(Nnum,Nnum), zeros(Nnum,1);
           eye(Nnum), zeros(Nnum, 2*Nnum+1);
           zeros(1,3*Nnum), 1];
            CoL = [zeros(Nnum,Nnum), eye(Nnum), zeros(Nnum,Nnum), zeros(Nnum,1);
           zeros(Nnum,Nnum), zeros(Nnum,Nnum), eye(Nnum), zeros(Nnum,1)];
        end
        K = [E; An(:,:,i)];
        LMI1{i} = [Q{i} - K' * P{i} * K, zeros(Nnum), S{i};
                zeros(Nnum), zeros(Nnum), zeros(Nnum);
                S{i}', zeros(Nnum), R{i}];
        
        LMI1{i} = LMI1{i} + blkdiag(zeros(Nnum),P{i});
    
        PiL = [Q{i}, S{i}; S{i}', R{i}];
        PiL = CoL' * PiL * CoL;
    
        % PiN 构造
        if i==1
            PiN = zeros(2*Nnum+1+in, 2*Nnum+1+in);
        else
            PiN = zeros(3*Nnum+1, 3*Nnum+1);    
        end
        for j = 1:Nnum
            for k = i:Nnum
                ej = E(:, j);
                ek = E(:, k);
                if j == k
                    e = ej;
                else
                    e = ej - ek;
                end
                M = [-2*(e*e'),  e*e';
                      e*e',     zeros(Nnum)];
                PiN = PiN + L{i}(j,k) * CoN' * M * CoN;
            end
        end
    
        LMI2{i} = PiL + PiN - CoP'*Pi{i}*CoP;
    end
    
    index = [out, Nnum*ones(1,Lnum-1),in,1];
    endi = cumsum(index);
    starti = cumsum([1 index(1:end-1)]);
    
    blk(starti(1):endi(1),starti(1):endi(1))==P11{Lnum};
    blk(starti(1):endi(1),starti(2):endi(2))==P12{Lnum};
    blk(starti(1):endi(1),starti(Lnum+2):endi(Lnum+2))==P13{Lnum};
    blk(starti(Lnum+1):endi(Lnum+1),starti(Lnum+1):endi(Lnum+1))==P22{1};
    blk(starti(Lnum+1):endi(Lnum+1),starti(Lnum+2):endi(Lnum+2))==P23{1};
    P33temp=P33{1};
    for i=2:Lnum
        P33temp=P33temp+P33{i};
    end
    blk(starti(Lnum+2):endi(Lnum+2),starti(Lnum+2):endi(Lnum+2))==P33temp-eye(1-starti(Lnum+2)+endi(Lnum+2));
    for i=2:Lnum
        j=Lnum+1-i;
        blk(starti(i):endi(i),starti(i):endi(i))==P11{j}+P22{j+1};
        blk(starti(i):endi(i),starti(i+1):endi(i+1))==P12{j};
        blk(starti(i):endi(i),starti(Lnum+2):endi(Lnum+2))==P13{j}+P23{j+1};
    end
    for i=1:Lnum-1
        for j=i:Lnum-1
            jj=j+2;
            blk(starti(i):endi(i),starti(jj):endi(jj))==zeros(1+endi(i)-starti(i),1+endi(jj)-starti(jj));
        end
    end

    Aby=[Ay;zeros(n+in-Nnum,1);by];
    LMI3=[blk+l1 * Coin' * Piin * Coin,Aby;Aby' -1];
    maximize(Ay)
    subject to
        for i=1:Lnum
            LMI1{i} >= 0
            LMI2{i} <= 0
        end
        LMI3<=0
        for i=1:Lnum
            for j=1:Nnum
                for k=j:Nnum
                    L{i}(j,k) >= 1e-8
                end
            end
        end
cvx_end
    end
    time = toc;
    value = Ay;
end
