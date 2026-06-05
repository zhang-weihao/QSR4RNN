clc;clear all;
Nnum = 2; in = 2; out = 2;
index = [Nnum*Nnum Nnum*in Nnum*1 out*Nnum out*1 Nnum*Nnum Nnum*Nnum Nnum*Nnum 4*Nnum*Nnum Nnum*Nnum Nnum in];
iii=sum(index);
rng('shuffle');
s=rng;
seed0=s.Seed;
% rng(1737784809);
rng(450853087);
x0 = rand(sum(index),1);

ii=sum([Nnum*Nnum Nnum*in Nnum*1 out*Nnum out*1 Nnum*Nnum Nnum*Nnum Nnum*Nnum 4*Nnum*Nnum out*1]);
lb=[-inf*ones(ii,1);zeros(iii-ii,1)];
ub=inf*ones(iii,1);

options = optimoptions('fmincon',...
    'StepTolerance',1e-18,...
    'MaxFunctionEvaluations', 1e6,...
    'MaxIterations', 1e6,...
    'MaxFunctionEvaluations', 1e6);

[u, u_fun, ~, seed1] = generateMPRS(0, 100000, 10, 0:0.01:1, 3000, 4000,425461574);
t0=0;tf=100000;step=10;
[~,yreal1]=ode45(@(t,x) sys(t,x,u_fun),t0:step:tf,[0;0]);
yrealmean=mean(yreal1);
mm=max(max(abs(yreal1-yrealmean)));
yrealnorm=(yreal1-yrealmean)/mm;

endi = cumsum(index);
starti = cumsum([1 index(1:end-1)]);
xx = fmincon(@(x) err(x,u,yrealnorm),x0,[],[],[],[],lb,ub,[],options);

[A,B,b,C,d,~,~,~,~,~,~,~] = unpack(xx,Nnum,in,out);
CC=mm*C;
dd=mm*d+yrealmean';
[u, u_fun, ~, seed2] = generateMPRS(0, 50000, 10, 0:0.1:1, 3000, 4000,425463376);
t0=0;tf=50000;step=10;
yyy=RNN(u,A,B,b,CC,dd,Nnum,out);
[~,yreal]=ode45(@(t,x) sys(t,x,u_fun),t0:step:tf,[0;0]);
FIT = 100*(1-norm(yyy-yreal')/norm(yreal'));
disp(FIT)

t=t0:step:tf;
figure;
subplot(1,2,1)
plot(t,yyy(1,:),'b',t,yreal(:,1)','c','LineWidth',1.2);
legend('Room 1(RNN)','Room 1(Real)')
xlabel('Time(s)','fontsize',12,'fontname','Times');
ylabel('Temperature (°C)','FontSize', 12,'fontname','Times');
grid on;

subplot(1,2,2)
plot(t,yyy(2,:),'b',t,yreal(:,2)','c','LineWidth',1.2);
legend('Room 2(RNN)','Room 2(Real)')
xlabel('Time(s)','fontsize',12,'fontname','Times');
ylabel('Temperature (°C)','fontsize',12,'fontname','Times');
grid on;

a=0.1;bbb=0.008;r=0.004;te=-10;th=50;
pi=[0 2*a r*th 0 bbb*te;2*a 0 0 r*th bbb*te;r*th 0 0 0 0;0 r*th 0 0 0;bbb*te bbb*te 0 0 0];
z=zeros(1,5001);
temp=[yyy(:,1);u(:,1);1];
z(1)=temp'*pi*temp;
for i=2:5001
    temp=[yyy(:,i);u(:,i);1];
    z(i)=z(i-1)+temp'*pi*temp;
end
figure;
plot(t,z,'b','LineWidth',2);
xlabel('Time(s)','fontsize',12,'fontname','Times');
ylabel('Safety constraint value','fontsize',12,'fontname','Times');
grid on;
