clc
clear all
close all

% Define nominal system and system with modelling error, then discretize

tau = 0.000096; % Absolute value of time constant

A = 1;

Bp = [-1 -1]*tau

Bc = [-1 -1]*tau;

C = 1;

ts = 0.1;
fs = 1/ts;

nomsys = ss(A,Bp,C,[]);
consumsys = ss(A,Bc,C,[]);

dSys = c2d(nomsys,ts,'tustin');
dCon = c2d(nomsys,ts,'tustin');

n = size(dSys.A,1);
m = size(dSys.B,2);
y = size(dSys.C,1);

x0 = zeros(n+y,1);
ref = 5;


% Construct velocity-form system matrices:

Av = [dSys.A zeros(n,y); dSys.C*dSys.A eye(y,y)];
Bv = [dSys.B ; zeros(y,m)];
Cv = [zeros(y,n) eye(y,y)];

VSys = ss(Av,Bv,Cv,[],ts);

% Make an LQR gain matrix

Q = Cv'*Cv; % Reference deviation cost
R = 1.*eye(m,m); % Actuation cost

[K,P,e] = lqr(VSys,Q,R);

pObs = 0.2*eig(Av-Bv*K);

L = place(Av',Cv',[0.1 0.4]);

%%


x = zeros(n+1,1);
dU(1:2,1) = 0; % Control input delta
x_real = zeros(n,1);
refval(1) = 0.1;
uLQR(1:2,1) = 0;
t_end = 10;
yLQR(:,1) = [0;0];                                                                                 5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
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
 
   
%    y_real(ii+1,1) = dSys.C*x_real(:,ii+1)+dCon.B*dc;
   x_real(:,ii+1) = dSys.A*x_real(:,ii)+dSys.B*dU(:,ii+1);
   y_real(ii+1,1) = dSys.C*x_real(:,ii+1)+dCon.B*dc(:,ii);
   x(end) = y_real(ii+1,1)-refval(ii);
   
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
plot(t,yLQR)
title('Velocity-form system')
legend('Tracking error','interpreter','latex')
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