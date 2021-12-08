clear all
close all

ts = 1;
t = 0:ts:(3600*48)-ts;

T1 = 1/(1.15495*10^-5); 
T2 = 1/(2.30989*10^-5);
T3 = 1/(3.46484*10^-5);
T4 = 1/(4.63291*10^-5);

max_flow_lab = 0.3; % max consumer flow i lab

w = 2*pi*1/(T1); %computes correct w given length of t, and choosen ts, such that we have two peaks in 24 hours

n = 4; % nr. of harmonics
Akalm = [0 w; -w 0];
AkalmBlk = blkdiag(Akalm,2*Akalm);
AkalmBlk = blkdiag(AkalmBlk,3*Akalm);
AkalmBlk = blkdiag(AkalmBlk,4*Akalm);
AkalmBlk = [zeros(1,2*n);AkalmBlk];
AkalmBlk = [zeros(2*n+1,1) AkalmBlk];
Bkalm = ones(n*2+1,1);
Ckalm = [1 1 0];

for ii = 1:n-1
    Ckalm = [Ckalm [1 0]];
end


%% leak
m = 60000;              %sort of arbitrary position of leak
leak_development = 10;        %30 minutes of development from leak happens to full leakage
max_leak = 0.1;               %assuming nominal flow of 2, max leak is assumed to be 10% of nominal flow
leak_increment = max_leak/10;      %how fast leak develops

leak(1:m) = 0;                  %gradual leak vector, begins at m, ramp increase over 2000 samples, then constant
leak(m:m+leak_development) = 0:leak_increment:max_leak;
leak(m+leak_development:length(t)) = max_leak;
%% Simulation

% scalars to simulate model error (1 = no error, perfect model)
a_c1 = 1.2;   % changes the amplitude of the sineusiods with 20% 
a_c2 = 1.2;
a_c3 = 0.5;
a_c4 = 1;

k_c = 1.05; % multiplier to the time-shift of the consumer pattern
w_c = 1.12;

noise_bit = 0; % 1 = noise on, 0 = noise off

%Identify DC + 4 harmonic components:
K = 5.4274;
h1 = 2.84716;
h2 = 1.2767;
h3 = 0.411351;
h4 = 0.428664;
% h5 = 0.2731;

tot = K + h1+h2+h3+h4;

% faserne fra fft analysen i rad. Det kan også ses i det matlab dokument
% hvor modellen udledes fra dataet på consumerne.
ph1 = -5.2790;
ph2 = -5.0431;
ph3 = 5.2835;
ph4 = -2.0002;

time_shift = (31500-8825)*k_c; % den tid vi flytter toppen, for at passe med datasettets top, k_c bruges til uderligere faseforskydning -> usikkerhed
time_shift_clean = 31500-8825;

% with uncertencties
theta1 = time_shift*(-2*pi*1/(T1));
theta2 = time_shift*(-2*pi*1/(T2));
theta3 = time_shift*(-2*pi*1/(T3));
theta4 = time_shift*(-2*pi*1/(T4));

w1 = h1/tot*max_flow_lab * sin(2*pi*1/T1 * t + ph1+theta1);
w2 = h2/tot*max_flow_lab * sin(2*pi*1/T2 * t + ph2+theta2);
w3 = h3/tot*max_flow_lab * sin(2*pi*1/T3 * t + ph3+theta3);
w4 = h4/tot*max_flow_lab * sin(2*pi*1/T4 * t + ph4+theta4);

% whitout uncertenties
theta1_clean = time_shift_clean*(-2*pi*1/(T1));
theta2_clean = time_shift_clean*(-2*pi*1/(T2));
theta3_clean = time_shift_clean*(-2*pi*1/(T3));
theta4_clean = time_shift_clean*(-2*pi*1/(T4));

w1_clean = h1/tot*max_flow_lab * sin(2*pi*1/T1 * t + ph1+theta1_clean);
w2_clean = h2/tot*max_flow_lab * sin(2*pi*1/T2 * t + ph2+theta2_clean);
w3_clean = h3/tot*max_flow_lab * sin(2*pi*1/T3 * t + ph3+theta3_clean);
w4_clean = h4/tot*max_flow_lab * sin(2*pi*1/T4 * t + ph4+theta4_clean);

ConsumerPattern = (K/tot)*max_flow_lab+a_c1*w1+a_c2*w2+a_c3*w3+a_c4*w4; % ændre vores "measurement", for at se om kalmanfiltret kan følge.
ConsumerPattern_clean = (K/tot)*max_flow_lab+w1_clean+w2_clean+w3_clean+w4_clean;

figure(1)
plot(t,ConsumerPattern) % model vi går efter
hold on
plot(t,ConsumerPattern_clean) % clean consumer pattern
title('Model of consumer pattern')
xlabel('Time [sec]')
ylabel('Consumption scaled to lab')
legend('With uncertenties','clean')
xlim([0 3600*48]) % consumer for 2 dag
hold off


%%
itteration = 1;
% x(:,1) = [K/tot*max_flow_lab h1/tot*max_flow_lab*cos(ph1*2*pi*1/(T1)-pi/2) h1/tot*max_flow_lab*sin(ph1*2*pi*1/(T1)-pi/2) h2/tot*max_flow_lab*cos(ph2*2*pi*1/(T2)-pi/2) h2/tot*max_flow_lab*sin(ph2*2*pi*1/(T2)-pi/2) h3/tot*max_flow_lab*cos(ph3*2*pi*1/(T3)-pi/2) h3/tot*max_flow_lab*sin(ph3*2*pi*1/(T3)-pi/2) h4/tot*max_flow_lab*cos(ph4*2*pi*1/(T4)-pi/2) h4/tot*max_flow_lab*sin(ph4*2*pi*1/(T4)-pi/2)]';
x(:,1) = [K/tot*max_flow_lab h1/tot*max_flow_lab*sin(ph1+theta1_clean) h1/tot*max_flow_lab*cos(ph1+theta1_clean) h2/tot*max_flow_lab*sin(ph2+theta2_clean) h2/tot*max_flow_lab*cos(ph2+theta2_clean) h3/tot*max_flow_lab*sin(ph3+theta3_clean) h3/tot*max_flow_lab*cos(ph3+theta3_clean) h4/tot*max_flow_lab*cos(ph4+theta4_clean) h4/tot*max_flow_lab*sin(ph4+theta4_clean)]';
estimate(:,1) = Ckalm * x(:,1)+noise_bit*normrnd(0,0.01);    % start værdi for estimatet
meassurement(1) = ConsumerPattern(1);               % start værdi for measurement
% for jj = 1000:1000:2000
Q = eye(1)*0.1;                                     % Uncertainty on system model   
R = eye(1)*70000;                                   % Uncertainty on observation
% R = eye(1)*jj;
kf = lqr(AkalmBlk',Ckalm',Q,R);                     % kalman gain
kf = 0;
for ii = 2:ts:length(t)
     %meassurement(ii) =  ConsumerPattern(ii)+0.05*(rand(1,1))+leak(ii); % med leak
     meassurement(ii) =  ConsumerPattern(ii)+noise_bit*normrnd(0,0.01); % uden leak

    x(:,ii) = x(:,ii-1) + ts*AkalmBlk*x(:,ii-1);
    estimate(:,ii) = Ckalm * x(:,ii);
    
    x(:,ii) = x(:,ii) + kf'*(meassurement(:,ii)-estimate(:,ii));
    estimate(:,ii) = Ckalm * x(:,ii);
    
    if ii > 5000000
        delta(ii) = estimate(ii) - (sum(meassurement(ii-5:ii))/5);
    else
        delta(ii) = abs((estimate(ii) - meassurement(ii)));
    end

end

% maxleakerror = max(delta(length(t)*0.6:length(t)*0.6+leak_development));        %search for max error in interval of leak development
% maxmodelerror_beforeleak = max(delta(1:length(t)*0.6));                             %search for max error before leak
% maxmodelerror_afterleak = max(delta((length(t)*0.61+leak_development):length(t)));        %search for max error after leak
% maxmodelerror = max([maxmodelerror_beforeleak maxmodelerror_afterleak]);               %find max error of before and after search
% maxilars(itteration) = maxleakerror-maxmodelerror;                                 %computes error difference
% noise_value(itteration) = jj;                                                 %stores value of observation uncertainty at given error difference 
% itteration = itteration + 1                                                  %increment itteration counter
% end

figure()
plot(ConsumerPattern_clean, 'r')
hold on
plot(estimate, 'b')
legend('measurement','estimate')
hold off
% figure()
% plot(delta)
% figure()
% plot(maxilars)

