% DRTS Assignment 1
clc;clear;close all;


% 1: Create a DNC/RTC arrival-model of the periodic external traffic sources 
% of the in-car network (choose arrival curve types and determine parameters)
% ------------------------------------------------------------------------
% 2: Create a DNC/RTC model for your communication network (Identify network 
% elements, select curve types and parameters)
% ------------------------------------------------------------------------


BW = 1*10^3; % [bpms] Serialization speed

% Wheel:
N_W = 4;	% Amount of wheel sensors   
p_W = 20;	% Amount of bits per sensor per period
O = 0;		% Overhead
T = 10;		% [ms] Transfer period for both wheel sensor and ESP

n = 1;		% Amount of sensor packets in each dataframe (CAN packet)
k_W = p_W*n+O;
tau = 0;	% jitter

% Arrival curves:
r_W = k_W/(n*T);			% Rate
b_W = k_W * (tau/T + 1);	% Burst
arr_W = rtccurve([0, b_W, r_W]);
figure
rtcplot(arr_W,T)
title('Arrival Curve for one Wheel sensor')
xlabel('time [ms]')
ylabel('bits')

r_W_all = N_W*r_W;
b_W_all = N_W*b_W;
arr_W_all = rtccurve([0, b_W_all, r_W_all]);
figure
rtcplot(arr_W_all,T)
title('Arrival Curve for all Wheel sensors')
xlabel('time [ms]')
ylabel('bits')


% ESP
p_ESP = 8;			% Amount of bits per sensor per period

k_ESP = p_ESP+O;

% Arrival:	
r_ESP = k_ESP/(n*T);			% Rate
b_ESP = k_ESP * (tau/T + 1);	% Burst

arr_ESP = rtccurve([0, b_ESP, r_ESP]);
figure
rtcplot(arr_ESP,T)
title('Arrival Curve ESP')
xlabel('time [ms]')
ylabel('bits')

% Combination of arrival constraints:
r_All = r_W_all + r_ESP;
b_All = b_W_all + b_ESP;

figure
arr_All = rtccurve([0, b_All, r_All]);
rtcplot(arr_All, T)
title('Arrival Curve for all Method 2')
xlabel('time [ms]')
ylabel('bits')

% Service curves:

d = (n-1)*T; % Delay - 0 since we are not packing more than one sensor
serv = rtccurve([d, 0, BW]);
figure
rtcplot(serv,T)
title('Service curve')
xlabel('time [ms]')
ylabel('bits')

% Backlog and delay:
figure
rtcplot(arr_All, 0.1) % Arrival curve
hold on
rtcplot(serv) % Service curve
rtcplotv(arr_All,serv) % Backlog
rtcploth(arr_All,serv, 'b') % Delay
title('Delay and backlog from arrival and service curves')
xlabel('time [ms]')
ylabel('bits')
legend('Arrival', 'Service', 'backlog', 'delay')

delay = rtch(arr_All,serv)
backlog = rtcv(arr_All,serv)

findfigs % Pull figures into screen


% 3. Guess initial parameters for token bucket filters (for the Poisson 
% traffic sources) and include parameters in the DNC/RTC model.
% ------------------------------------------------------------------------

% - Packet size is defined in task as 1400 bits
% - Average period is set to 30 ms (although its a poisson process)

% Rear Camera parameters
T_RC = 40; % [ms]
p_RC = 1400; % Packet size - Number of bits per packet
% Multimedia parameters
T_MM = 40; % [ms]
p_MM = 1400; % Packet size - Number of bits per packet

% Token bucket filter parameters - initial guesses
M_TB = 5; % Bucket size - Burst parameter (random guess=
T_TB = 30; % [ms] Token replenishment rate - Curve rate (Bit faster than 40ms)


% 4. Determine arrival curves for outputs of token bucket filters and 
% include in the DNC/RTC model.
% ------------------------------------------------------------------------


% Arrival curves:
r_RC = 1/T_TB * p_RC;
b_RC= p_RC * M_TB;
arr_W = rtccurve([0, b_RC, r_RC]);
figure
rtcplot(arr_W,T)
title('Arrival Curve for one token bucket filter for rear camera')
xlabel('time [ms]')
ylabel('bits')

% 5. Compute max backlogs and max waiting times for all flows for the 
% (non)deterministic part of the network
% ------------------------------------------------------------------------

