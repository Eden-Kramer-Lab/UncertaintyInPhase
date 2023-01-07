%% Fitzhugh-Nagumo Model
% dx/dt = x(x-1)(1-a1x) - y + I
% dy/dt = bx
%% chaotic regime
a1 = 10;
b = 1;
a = 0.1;
f = 0.13;
omega = 2*pi*f;

%%
% let y(1) = u; y(2) = v
rng('shuffle')
Fs = 50;
F = @(t,y)  [y(1)*(y(1)-1)*(1-a1*y(1)) - y(2)+ (a/(omega))*cos(omega*t)+randn();... % +10*randn()
            b*y(1)];
tspan = [1/Fs : 1/Fs: 500];              % time
y0 = [.2 .2];                 % u0 = 1; v0 = 2;
[t,y] = ode45(F,tspan,y0);
% plot(1/1000:1/1000:25,y(:,1))
% legend('u(t)','v(t)')