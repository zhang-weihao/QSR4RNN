clc;clear all;close all;

n = 100;
value_err = zeros(n,1);
t1 = zeros(n,1);
t2 = zeros(n,1);

Nnum = 5;
in = 1;
out = 1;
Piin = [-1,0;0,1];

colors = [
    0.85 0.33 0.10;  
    0.00 0.45 0.74;
    0.49 0.18 0.56;  
];
figure;

for l = 4:6
    Lnum = l;
    for i = 1:n
        An = zeros(Nnum,Nnum,Lnum);
        for k=1:Lnum
            [Q,~]=qr(randn(Nnum));
            D=0.5 * diag(2*rand(Nnum,1)-1);
            An(:,:,k)=Q*D*Q';
        end
        B1 = randn(Nnum, 1);
        Bn = randn(Nnum, Nnum, Lnum-1);
        bn = randn(Nnum, 1, Lnum);
        C = randn(out, Nnum);
        d = randn(out,1);
        [value1,time1]=QSRopt(An,B1,Bn,bn,C,d,Piin,1);
        [value2,time2]=QSRopt(An,B1,Bn,bn,C,d,Piin,2);
        value_err(i) = (value1-value2)/value2*100;
        t1(i) = time1;
        t2(i) = time2;
    end
    subplot(3,3,3*l-11)
    histogram(value_err, 50, 'FaceColor', colors(1,:));
    xlabel('$$(A_{y_1}-A_{y_2})/A_{y_2}\times 100$$', 'Interpreter', 'latex','fontsize',15,'fontname','Times');
    title(['$$\ell=$$',num2str(l)], 'Interpreter', 'latex', 'FontWeight','bold','fontsize',15,'fontname','Times');
    grid on;
    subplot(3,3,3*l-10)
    histogram(t1, 25, 'FaceColor', colors(2,:));
    xlabel('Running Time 1 (sec)', 'fontsize',15,'fontname','Times');
    title(['$$\ell=$$',num2str(l)], 'Interpreter', 'latex', 'FontWeight','bold','fontsize',15,'fontname','Times');
    grid on;
    subplot(3,3,3*l-9)
    histogram(t2, 25, 'FaceColor', colors(3,:));
    xlabel('Running Time 2 (sec)', 'fontsize',15,'fontname','Times');
    title(['$$\ell=$$',num2str(l)], 'Interpreter', 'latex', 'FontWeight','bold','fontsize',15,'fontname','Times');
    grid on;
end



