function xdot = sys(t, x, u_fun)
a=0.05;b=0.008;r=0.004;te=-10;th=50;
u=u_fun(t);

xdot = zeros(2,1);

xdot(1) = 2*a*(x(2)-x(1)) + b*(te-x(2)) + r*(th-x(1))*u(1);
xdot(2) = 2*a*(x(1)-x(2)) + b*(te-x(1)) + r*(th-x(2))*u(2);