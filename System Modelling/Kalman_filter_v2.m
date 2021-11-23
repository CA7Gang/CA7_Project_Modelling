% clear all
% med vores ss system.
A = [-0.3236   -0.0406   -0.1577   10.0671   -0.7623    0.0000;
   -0.1429   -0.3189    0.2176   -1.3489   -6.0984         0;
   -0.0275    0.0687   -0.3968    7.1758    1.5246         0;
    0.1089    0.0687    0.0196  -20.2792    1.5246   -0.0000;
   -0.0551   -0.0443   -0.0272    1.4425  -15.3466   -0.0999;
    0.1486    0.0817   -0.2877   11.1436    8.2419   -0.4174;];

B = [0.0982    0.0000;
    0.0078         0;
    0.1147         0;
   -0.0408         0;
   -0.0087   -0.0193;
   -0.0651   -0.0715];

C = [0 0 1 0 0 0; 0 0 -1 -1 -1 -1];

D = zeros(size(C,1),size(B,2));

vd = 0.05*eye(6); % disturbence covarianc matrix
vn = 0.1; % measurement noise, scalar, typically bigger than vd
BF = [B vd 0*B] % augment B matrix, wtih disturbence and noise (0 noise)
Daug = [0 0 0 0 0 0 0 0 0 vn;
        0 0 0 0 0 0 0 0 0 vn]
sysC = ss(A,BF,C,Daug)
sysFullOutput = ss(A,BF,eye(6),zeros(6,size(BF,2)))

%[Kf,P,E] = lqe(A,vd,C,vd,vn)
kf = (lqr(A',C',vd,vn))'
sysKf = ss(A-kf*C,[B kf],eye(6),0*[B kf])

% uDist = randn(6,length(t));
% uNoise = randn(2,length(t));
% 
% u = ones(length(t),2)*66;
% 
% uAug = [u'; vd*vd*uDist; uNoise];

%% simulering af ulineært system

ts = 0.01;



f_p = 2/(24*10); % Assumed pump frequency



t = 0:ts:(5*60)-ts; % Time vector corresponding to 24 hours
% w = 50*(1+square(12*pi*f_p*t,50))/2+50*(1+square(12*pi*f_p*t,100))/2; % Pump waveforms
% OD =(1+sin(pi*f_p*t+t(end)/2))/2; % Valve waveforms


clear flow tankpres df d_t pt w1 w2 OD1 OD2 pressures q_lin pt_lin flowchange linpump
w1 = 66*ones(length(t),1); w2 = 66*ones(length(t),1);
% w1 = w; w2 = w; 
% OD1 = OD; OD2 = OD;
OD1 = 0.5*ones(length(t),1); OD2 = 0.5*ones(length(t),1);
% qc = [0;0];
qc = randn(2,1);
% qc = q0(1:2);
% df = [q0(3:end-1);0];
% df = [0;0;0;0]; 
df = randn(4,1);
pt =0*4*rho*g/10^5;
% d_t = 0;
d_t = 0;

w0 = 66;

LinSys = ss(A,B,C,[]);
LinSys = c2d(LinSys,ts,'Tustin');
u = [w0;w0];
[kalmf, L, P] = kalman(LinSys,1,1,0);
x = zeros(6,1) 
q0 = double(q0)

for ii = 1:length(t)
    [dqdt,pbar,pt_new] = fooSim.Model_TimeStep([w1(ii) w2(ii)],[OD1(ii) OD2(ii)],[qc(1) qc(2)],[df(1) df(2) df(3) df(4)],d_t,pt);
    qc(1) = qc(1)+dqdt(1)*ts;
    qc(2) = qc(2)+dqdt(2)*ts;
    df(1) = df(1)+dqdt(3)*ts; 
    df(2) = df(2)+dqdt(4)*ts; 
    df(3) = df(3)+dqdt(5)*ts; 
    d_t = d_t+dqdt(end)*ts;

    % Must obey basic mass conservation if no leaks are assumed
    masscon = -cumsum([df(1:3);d_t]);
    df(4) = masscon(end);

    pt = pt_new;
    flow(:,ii) = [qc;df;d_t];
    flowL(:,ii) = [qc;df(1:3);d_t];
    pressures(:,ii) = double(pbar);
    tankpres(ii) = pt;

    % Subtracts the eqpoint flows, to get the rael flow
    x_true(:,ii) = [flow(3,ii)-q0(3);flow(6,ii)-sum(q0(3:6))]';
    
if ii == 1
    %estimate states
    x(:,ii) = LinSys.B*u;
else
    %estimate states
    x(:,ii) = LinSys.A*x(:,ii-1)+LinSys.B*u;
end
    z(:,ii) = LinSys.C*x(:,ii);
    %estimate observation
    
    x(:,ii) = x(:,ii)+L*(x_true(:,ii)-z(:,ii));


%

%


%     y_simulering = [flow(3,:);flow(6,:)]';
% 
%     x_hat = sysKf, [u'; y_simulering'],t;
%     plot(t,x_hat(:,1),'k--','LineWidth',2.0)
%     hold off
%     
%     figure()
%     plot(t,flowL,'-',t,x_hat,'--','LineWidth',2)
%     legend('q3','q7','d1','d5','d9','d8')
end    
   
% for ii = 1:length(q_lin)
%     masscon = cumsum(q_lin(3:end,ii));
%     linpump(ii) = -(masscon(end));
% end

% t = seconds(t);
% t = minutes(t);
%% plot

plot(t,flow) % simuleret resulatat
legend('chord1','chord2','pump1','cunsumer1','consumer2','pump2','tankflow')

%% kalman filter

vd = 0.05*eye(6) % disturbence covarianc matrix
vn = 0.1; % measurement noise, scalar, typically bigger than vd
BF = [B vd 0*B] % augment B matrix, wtih disturbence and noise (0 noise)
Daug = [0 0 0 0 0 0 0 0 0 vn;
        0 0 0 0 0 0 0 0 0 vn]
sysC = ss(A,BF,C,Daug)
sysFullOutput = ss(A,BF,eye(6),zeros(6,size(BF,2)))

%[Kf,P,E] = lqe(A,vd,C,vd,vn)
kf = (lqr(A',C',vd,vn))'
sysKf = ss(A-kf*C,[B kf],eye(6),0*[B kf])

uDist = randn(6,length(t));
uNoise = randn(2,length(t));

% dette skal laves ordenligt!!! Hvad skal vores kontrol signal være for at
% få pumpen til at køre 66% ? er det bare 66? eller skal det være noget
% helt andet?
% u = [0;0]*t; 
% u(100:120) = 66;
% u(1500:end) = 66;

u = ones(length(t),2)*66;

uAug = [u'; vd*vd*uDist; uNoise];

% [y,t] = lsim(sysC,uAug,t);
% plot(t,y)

% [xtrue,t] = lsim(sysFullOutput,uAug,t);
% hold on
% plot(t,xtrue(:,1),'r','LineWidth',2.0)

testy = [flow(3,:);flow(6,:)]'; % ved ikke helt...

[x_hat,t] = lsim(sysKf, [u'; testy'],t);
plot(t,x_hat(:,1),'k--','LineWidth',2.0)
hold off

figure()
plot(t,flowL,'-',t,x_hat,'--','LineWidth',2)
legend('q3','q7','d1','d5','d9','d8')

% Overvejelser til at sammenligne med ulineært system:
% Har vi den rigtige C matrix?
% brug samme setup som i simuleringen i conponentUnitTest?
% Udskrift xtrue, med flows fra simuleringen i component unit test
% Skal vi have et nyt sys (sysFullOutput) hvor det andet pumpe flow også er med? i hvert fald laves om?

%% tro kopi af det eksempel som Steve Brunton laver på youtube.
m = 1;
M = 5;
L = 2;
g = -10;
d = 1;

s = -1;

A_ex = [0 1 0 0;
       0 -d/M -m*g/M 0;
       0 0 0 1;
       0 -s*d/(M*L) -s*(m+M)*g/(M*L) 0]

B_ex = [0; 1/M; 0; s*1/(M*L)]

C_ex = [1 0 0 0]

D_ex = zeros(size(C_ex,1),size(B_ex,2))

vd = 0.1*eye(4) % disturbence covarianc matrix
vn = 1; % measurement noise, scalar, typically bigger than vd
BF_ex = [B_ex vd 0*B_ex] % augment B matrix, wtih disturbence and noise (0 noise)
Daug_ex = [0 0 0 0 0 vn]
sysC_ex = ss(A_ex,BF_ex,C_ex,Daug_ex)
sysFullOutput_ex = ss(A_ex,BF_ex,eye(4),zeros(4,size(BF_ex,2)))

kf_ex = (lqr(A_ex',C_ex',vd,vn))'
sysKf_ex = ss(A_ex-kf_ex*C_ex,[B_ex kf_ex],eye(4),0*[B_ex kf_ex])

dt = 0.01;
t = dt:dt:50;

uDist = randn(4,size(t,2));
uNoise = randn(size(t));
u = 0*t;
u(100:120) = 100;
u(1500:1520) = -100;

uAug = [u; vd*vd*uDist; uNoise];

[y_ex,t_ex] = lsim(sysC_ex,uAug,t);
plot(t_ex,y_ex)

[xtrue_ex,t] = lsim(sysFullOutput_ex,uAug,t);
hold on
plot(t,xtrue_ex(:,1),'r','LineWidth',2.0)

[x_hat_ex,t] = lsim(sysKf_ex, [u; y_ex'],t);
plot(t,x_hat_ex(:,1),'k--','LineWidth',2.0)

figure()
plot(t,xtrue_ex,'-',t,x_hat_ex,'--','LineWidth',2)
