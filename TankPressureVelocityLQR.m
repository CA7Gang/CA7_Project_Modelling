clc
clear all
close all

% Define nominal system and system with modelling error, then discretize

tau = 0.000096; % Absolute value of time constant

ts = 10;
fs = 1/ts;

T = tau*ts; % Euler discretization of the time constant

A = 1;

Bp = [T T];

Bc = [T T];


C = 1;


dNom = ss(A,2*Bp,C,[],ts);
dSys = ss(A,Bp,C,[],ts);
dCon = ss(A,Bc,C,[],ts);

n = size(dSys.A,1);
m = size(dSys.B,2);
y = size(dSys.C,1);

x0 = zeros(n+y,1);


% Construct velocity-form system matrices:

Av = [dSys.A zeros(n,y); dSys.C*dSys.A eye(y,y)];
Bv = [dSys.B ; dSys.C*dSys.B];
Cv = [zeros(y,n) eye(y,y)];
Bdist = [eye(n,n);dSys.C];

VSys = ss(Av,Bv,Cv,[],ts);

% Make an LQR gain matrix

Q = Cv'*Cv; % Reference deviation cost
Qnom = C'*C;
% Q = 0.01*eye(2,2);
R = 0.1*eye(m,m); % Actuation cost

[K,P,e] = lqr(VSys,Q,R);
[Kn,Pn,en] = lqr(dSys,Qnom,R);

Kff = inv(R)*dSys.B'*inv(dSys.A-Pn*dSys.B/R*dSys.B')*Pn*dSys.B

CLSysV = ss(Av-Bv*K,Bv,Cv,[],ts);

% step(CLSysV,100)

pObs = 0.2*eig(Av-Bv*K);

L = place(Av',Cv',[0.1 0.4]);

rank(ctrb(Av,Bv))

[Av_bar, Bv_bar, Cv_bar, T, k] = ctrbf(Av,Bv,Cv);
n_uc = size(Av, 1) - sum(k); % Number of uncontrollable modes is 8 - 6 = 2
Av_uc = Av_bar(1:n_uc, 1:n_uc)

%%


x = zeros(n+1,1);
dU(1:2,1) = 0; % Control input delta
x_real = zeros(n,1);
refval(1) = 1;
uLQR(1:2,1) = 0;
t_end = 4*3600;
yLQR(:,1) = [0;0];                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
t = 0:ts:(t_end-ts);

e = -2+0.5*sin(2*pi*12/max(t)*t);
dc = [e;e];

% dc = -1*ones(2,length(t));

dist = dSys.B*dc;

plot(dist(1,:))

t = seconds(t);

T = ceil(length(t));

%%

for ii = 1:length(t)-1   
    
   dU(:,ii+1) = uLQR(:,ii)+dU(:,ii);  
   
   u(:,ii) = dU(:,ii+1)-dc(:,ii); % Just cancel the fucking disturbance out 4Head
    
   x(:,ii+1) = Av*x(:,ii)+Bv*uLQR(:,ii); 
   yLQR(ii+1) = Cv*x(:,ii+1);
 
 
   x_real(:,ii+1) = dSys.A*x_real(:,ii)+dSys.B*u(:,ii)+dSys.B*dc(:,ii);
   y_real(ii+1) = dSys.C*(x_real(:,ii+1));
   x(end,ii+1) = y_real(ii+1)-refval(ii);
   
   uLQR(:,ii+1) = -K*(x(:,ii+1));
   
   
%    if (ii > ceil(T/2))
% 
%        refval(ii+1) = 2;
%    else
%        refval(ii+1) = refval(ii);
%    end
   
    refval(ii+1) = refval(ii);
end

%%

t = hours(t);

figure(2)
subplot(2,2,1)
ylabel('Output')
xlabel('Samples')
plot(t,x(1,:))
hold on
plot(t,x(2,:))
title('Velocity-form system')
legend('Pressure change','Tracking error','interpreter','latex')
xlabel('Time [hr]')

subplot(2,2,2)
plot(t,uLQR(1,:))
hold on
plot(t,uLQR(2,:))
hold off
title('Differential control input')
xlabel('Time [hr]')

subplot(2,2,3)
ylabel('Output')
xlabel('Samples')
plot(t,y_real)
hold on
plot(t,refval,'--r')
hold off
title('Real System')
legend('Process value','Reference','interpreter','latex')
xlabel('Time [hr]')

subplot(2,2,4)
plot(t(1:end-1),u(1,:),'b')
hold on
plot(t(1:end-1),u(2,:),'r')
plot(t(1:end-1),dc(1,1:end-1),'--k')
plot(t(1:end-1),dc(2,1:end-1),'--m')
hold off
title('Full control input')
xlabel('Time [hr]')