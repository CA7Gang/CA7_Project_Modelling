clc
clear all
close all

% Define nominal system and system with modelling error, then discretize

tau = -0.000096; % Absolute value of time constant

A = 1;

Bp = [-tau -tau];

Bc = [-1 -1];

C = 1;

ts = 1;
fs = 1/ts;

nomsys = ss(A,Bp,C,[]);
consumsys = ss(A,Bc,C,[]);

% dSys = ss(A,Bp,C,[],ts);
% dCon = ss(A,Bc,C,[],ts);

dSys = c2d(nomsys,ts);
dCon = c2d(consumsys,ts);

n = size(dSys.A,1);
m = size(dSys.B,2);
y = size(dSys.C,1);

x0 = zeros(n+y,1);


% Construct velocity-form system matrices:

Av = [dSys.A zeros(n,y); dSys.C*dSys.A eye(y,y)];
Bv = [dSys.B ; dSys.C*dSys.B];
Cv = [zeros(y,n) eye(y,y)];

VSys = ss(Av,Bv,Cv,[],ts);

% Make an LQR gain matrix

Q = Cv'*Cv; % Reference deviation cost
% Q = 0.01*eye(2,2);
R = eye(m,m); % Actuation cost

[K,P,e] = lqr(VSys,Q,R);
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
t_end = 100;
yLQR(:,1) = [0;0];                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
t = 0:ts:(t_end-ts);

% e = sin(2*pi*1/max(t)*t);
% dc = [e;e];

dc = [-0.0003;-0.00002].*ones(2,length(t))*0;

T = ceil(length(t));

%%

for ii = 1:length(t)-1   
    
   dU(:,ii+1) = uLQR(:,ii)+dU(:,ii);   
    
   x(:,ii+1) = Av*x(:,ii)+Bv*uLQR(:,ii); 
   yLQR(ii+1) = Cv*x(:,ii+1);
 
 
   x_real(:,ii+1) = dSys.A*x_real(:,ii)+dSys.B*dU(:,ii+1);
   y_real(ii+1) = dSys.C*x_real(:,ii+1);
   x(end,ii+1) = y_real(ii+1)-refval(ii);
   
   uLQR(:,ii+1) = -K*x(:,ii+1);
   
%    if (ii < ceil(2*T/3)) && (ii > ceil(T/3))
%        refval(ii+1) = 2.5;
%    elseif ii > ceil(2*T/3)
%        refval(ii+1) = 5;
%    else
%        refval(ii+1) = refval(ii);
%    end
   
    refval(ii+1) = refval(ii);
end

%%


figure(1)
subplot(2,2,1)
ylabel('Output')
xlabel('Samples')
plot(t,x(1,:))
hold on
plot(t,x(2,:))
title('Velocity-form system')
legend('Pressure change','Tracking error','interpreter','latex')
subplot(2,2,2)
plot(t,uLQR(1,:))
hold on
plot(t,uLQR(2,:))
hold off
title('Differential control input')

subplot(2,2,3)
ylabel('Output')
xlabel('Samples')
plot(t,y_real)
hold on
plot(t,refval,'--r')
hold off
title('Real System')
legend('Process value','Reference')
subplot(2,2,4)
plot(t,dU(1,:))
hold on
plot(t,dU(2,:))
hold off
title('Full control input')