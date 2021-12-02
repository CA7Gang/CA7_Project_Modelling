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

s = tf('s');

[num,den] = pade(4,1);

[Ad,Bd,Cd,Dd] = tf2ss(num,den);

FastSys = ss(A,B,C,D);

p = eig(FastSys.A);

L = place(FastSys.A',FastSys.C',4*p);



% Define nominal system and system with modelling error, then discretize

tau = 0.000096; % Absolute value of time constant

ts = 3*60;
fs = 1/ts;

T = tau*ts; % Euler discretization of the time constant

A = 1.00;

Bp = [T T];
 
Bc = [T T];


C = 1;


dNom = ss(A,Bp,C,[],ts);
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
% Q = 0.01*eye(2,2);
R = 10*eye(m,m); % Actuation cost

[K,P,e] = lqr(VSys,Q,R);

CLSysV = ss(Av-Bv*K,Bv,Cv,[],ts);

% step(CLSysV,100) 

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
t_end = 12*3600;
yLQR(:,1) = [0;0];                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
t = 0:ts:(t_end-ts);

e = -2+0.5*sin(2*pi*0.5/max(t)*t);
dc = [e;e];

% dc = -1*ones(2,length(t));

plot(dc(1,:))

t = seconds(t);

T = ceil(length(t));

%%

for ii = 1:length(t)-1   
    
    dU(:,ii+1) = uLQR(:,ii)+dU(:,ii);  

   u(:,ii) = dU(:,ii+1)-pinv(dSys.B)*dSys.B*dc(:,ii); % Just cancel the fucking disturbance out 4Head
   
   for jj = 1:2
       if u(jj,ii) < 0
           u(jj,ii) = 0;
       elseif u(jj,ii) > 3.5
           u(jj,ii) = 3.5;
       end
   end
    
   x(:,ii+1) = Av*x(:,ii)+Bv*uLQR(:,ii); 
   yLQR(ii+1) = Cv*x(:,ii+1);
 
 
   x_real(:,ii+1) = dSys.A*x_real(:,ii)+dNom.B*u(:,ii)+dNom.B*dc(:,ii);
%    y_real(ii+1) = dSys.C*(x_real(:,ii+1));
   y_real(ii+1) = dSys.C*(x_real(:,ii+1))-0.1;
   x(end,ii+1) = y_real(ii+1)-refval(ii);
   
   uLQR(:,ii+1) = -K*(x(:,ii+1));
   
   
   if (ii > ceil(T/2))

       refval(ii+1) = 2;
   else
       refval(ii+1) = refval(ii);
   end
   
%     refval(ii+1) = refval(ii);
end

%% Save stuff to an appropriate file

% save('NomSysData.mat','x','uLQR','u','y_real','refval','dc','t')
% save('HalfNomSysData.mat','x','uLQR','u','y_real','refval','dc','t')
% save('DoubleNomSysData.mat','x','uLQR','u','y_real','refval','dc','t')
% save('OutDistNomSysData.mat','x','uLQR','u','y_real','refval','dc','t')


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
legend('Pressure','Reference','interpreter','latex')
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

%% Transfer function representations
% 
% close all
% z = tf('z');
% 
% GLQR = K*inv(z*eye(2,2)-Av)*Bv;
% LQRtf = (GLQR);
% LQRtf.Ts = ts;
% 
% [num,den] = ss2tf(VSys.A,VSys.B,VSys.C,zeros(1,2),1);
% G1 = tf(num,den,ts);
% GOL = minreal(LQRtf*G1);
% GCL = minreal(GOL/(1+GOL))
% pzmap(GCL)

% Using data from Sheth et al, compute an upper bound on the delay in each
% scenario

Loss2 = 0.08; % Worst-case loss, 2 km distance urban
Loss8 = 0.15; % Worst-case loss, 8 km distance urban
Loss20 = 0.6; % Worst-case loss, 20 km distance urban


tkalm = 10;
n = 4;

w = 2*pi*1/(3600*12);
Akalm = [0 -w; w 0];
AkalmBlk = blkdiag(Akalm,2*Akalm);
AkalmBlk = blkdiag(AkalmBlk,4*Akalm);
AkalmBlk = blkdiag(AkalmBlk,8*Akalm);
AkalmBlk = [zeros(1,2*n);AkalmBlk];
AkalmBlk = [zeros(2*n+1,1) AkalmBlk];
Bkalm = ones(n*2+1,1);
Ckalm = [1 1 0];

for ii = 1:n-1
    Ckalm = [Ckalm [1 0]];
end

dAkalm = eye(2*n+1,2*n+1)+AkalmBlk*tkalm;

KalmSys = ss(dAkalm,Bkalm,Ckalm,[],tkalm);
