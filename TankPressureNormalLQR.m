clc
clear all
close all

% Fast dynamics

A=[-0.3236   -0.0406   -0.1577   10.0671   -0.7623    0.0000;
   -0.1429   -0.3189    0.2176   -1.3489   -6.0984         0;
   -0.0275    0.0687   -0.3968    7.1758    1.5246         0;
    0.1089    0.0687    0.0196  -20.2792    1.5246   -0.0000;
   -0.0551   -0.0443   -0.0272    1.4425  -15.3466   -0.0999;
    0.1486    0.0817   -0.2877   11.1436    8.2419   -0.4174];

B=[0.0982    0.0000;
    0.0078         0;
    0.1147         0;
   -0.0408         0;
   -0.0087   -0.0193;
   -0.0651   -0.0715];

C = [0 0 1 0 0 0;
	0 0 -1 -1 -1 -1];

D = zeros(2,2);

FastSys = ss(A,B,C,D);
%% Define nominal system and system with modelling error, then discretize

tau = -0.000096; % Absolute value of time constant

A = 1;

Bp = [tau tau];

Bc = [-1 -1];

C = 1;

ts = 100;
fs = 1/ts;

nomsys = ss(A,Bp,C,[]);
consumsys = ss(A,Bc,C,[]);

dSys = c2d(nomsys,ts);
dCon = c2d(consumsys,ts);

n = size(dSys.A,1);
m = size(dSys.B,2);
y = size(dSys.C,1);

x0 = zeros(n+y,1);


% Construct velocity-form system matrices:

Ai = [dSys.A zeros(n,y); dSys.C eye(y,y)];
Bi = [dSys.B ; zeros(n,m)];
Ci = [C zeros(y,y)];

ISys = ss(Ai,Bi,Ci,[],ts);

% Make an LQR gain matrix

Q = [0 0; 0 C'*C];; % Reference deviation cost
% Q = 0.01*eye(2,2);
R = 0.01*eye(m,m); % Actuation cost

[K,P,e] = lqr(ISys,Q,R);
CLSysV = ss(Ai-Bi*K,Bi,Ci,[],ts);

% step(CLSysV,100)

pObs = 0.2*eig(Ai-Bi*K);

% L = place(Ai',Ci',[0.1 0.4]);

rank(ctrb(Ai,Bi))

[Ai_bar, Bi_bar, Ci_bar, T, k] = ctrbf(Ai,Bi,Ci);
n_uc = size(Ai, 1) - sum(k); % Number of uncontrollable modes
Ai_uc = Ai_bar(1:n_uc, 1:n_uc)

%%


x = zeros(n+1,1);
dU(1:2,1) = 0; % Control input delta
x_real = zeros(n,1);
refval(1) = 1;
uLQR(1:2,1) = 0;
t_end = 24*3600;
yLQR(:,1) = [0;0];                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
t = 0:ts:(t_end-ts);

% e = sin(2*pi*1/max(t)*t);
% dc = [e;e];

dc = [-0.0003;-0.00002].*ones(2,length(t))*0;

T = ceil(length(t));

%%

for ii = 1:length(t)-1   
   
   x(:,ii+1) = Ai*x(:,ii)+Bi*uLQR(:,ii); 
   yLQR(ii+1) = Ci*x(:,ii+1);
 
   x(end,ii+1) = x(end,ii)+yLQR(ii+1)-refval(ii);
   
%    if x(end,ii+1) > 3.5
%        x(end,ii+1) = 3.5;
%    elseif x(end,ii+1) < -3.5
%        x(end,ii+1) = -3.5;
%    end
       
   
   uLQR(:,ii+1) = -K*x(:,ii+1);
   
   if (ii < ceil(2*T/3)) && (ii > ceil(T/3))
       refval(ii+1) = 2.5;
   elseif ii > ceil(2*T/3)
       refval(ii+1) = 5;
   else
       refval(ii+1) = refval(ii);
   end
%    
%     refval(ii+1) = refval(ii);
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
legend('Pressure','Integrator state','interpreter','latex')
subplot(2,2,2)
plot(t,uLQR(1,:))
hold on
plot(t,uLQR(2,:))
hold off
title('Control input')

subplot(2,2,3)
ylabel('Output')
xlabel('Samples')
plot(t,yLQR)
hold on
plot(t,refval,'--r')
hold off
title('Real System')
legend('Process value','Reference')