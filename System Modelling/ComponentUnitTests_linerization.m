clear all
close all

g = 9.82; % Gravitational acceleration
h0 = 0; % Height of the reference node
rho = 1000; % Density of water (kg/L)
v = 1.012*10^-6; % Kinematic viscosity of water
q_mean = 1.5/3600; % Mean assumed flow per second (1.5 cubic meters / hr)

kf = 1.8;

PipeIndex = [2 3 5 6 7 8 9 10 12 13];

PumpIndex = [1 14];

ValveIndex = [4 11];

PipeLen = [5 15 15 15 5 15 5 15 20 5]'; 
PipeDiam = 10^-3.*[25 25 25 25 25 15 25 25 15 15]';

eta = 5*10^-5.*ones(10,1); % Pipe roughness factors


% Pump coefficients
a2 = -[0.0355 0.0355]';
a1 = [0.0004 0.0004]'*0;
a0 = [0.0001 0.0001]';


Kv = ones(1,2);


for ii = 1:length(PipeLen)
Pipes(ii) =  PipeComponent(PipeLen(ii),rho,PipeDiam(ii),g,eta(ii),q_mean,v,kf);
end

for ii = 1:length(a2)
Pumps(ii) = PumpComponent(a0(ii),a1(ii),a2(ii));
end

for ii = 1:length(Kv)
Valves(ii) = ValveComponent('Linear',Kv(ii));
end

for ii = 1:length(PipeIndex)
    Components{PipeIndex(ii)} = @(q,w,OD) Pipes(ii).Lambda(q)+Pipes(ii).mu(q,OD)+Pipes(ii).alpha(q,w);
    Inertias(PipeIndex(ii)) = Pipes(ii).J;
end

for ii = 1:length(PumpIndex)
    Components{PumpIndex(ii)} = @(q,w,OD) Pumps(ii).Lambda(q)+Pumps(ii).mu(q,OD)+Pumps(ii).alpha(q,w);
    Inertias(PumpIndex(ii)) = Pumps(ii).J;
end

for ii = 1:length(ValveIndex)
    Components{ValveIndex(ii)} = @(q,w,OD) Valves(ii).Lambda(q)+Valves(ii).mu(q,OD)+Valves(ii).alpha(q,w);
    Inertias(ValveIndex(ii)) = Valves(ii).J;
end

Inertias = diag(Inertias);



% H = [1 0 1 0 0 0 0;
%     -1 1 0 0 0 0 0;
%     0 0 -1 1 1 0 0;
%     0 -1 0 -1 0 1 0;
%     0 0 0 0 -1 0 -1;
%     0 0 0 0 0 -1 1];

 
H = [1	0	0	0	0	0	0	0	0	0	0	0	0	0;
-1	1	0	0	0	0	0	0	0	0	0	0	0	0;
0	-1	1	0	1	0	0	0	0	0	0	0	0	0;
0	0	-1	1	0	1	0	0	0	0	0	0	0	0;
0	0	0	-1	0	0	0	0	0	0	0	0	0	0;
0	0	0	0	-1	0	1	1	0	0	0	0	0	0;
0	0	0	0	0	-1	-1	0	1	0	0	0	0	0;
0	0	0	0	0	0	0	0	-1	1	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	-1	0	0	0;
0	0	0	0	0	0	0	-1	0	0	1	-1	0	0;
0	0	0	0	0	0	0	0	0	-1	0	1	-1	0;
0	0	0	0	0	0	0	0	0	0	0	0	1	-1;
0	0	0	0	0	0	0	0	0	0	0	0	0	1];


 
% NodeHeights = zeros(size(H,1)-1,1)-h0;
NodeHeights = [0; 0; 0; 0; 0.9; 0; 0; 3; 0.9; 0; 0; 0]-h0; 
p0 = 0;
 
 
fooGraph = GraphModel(H,[3,7],13,[1 13],[5 9],8);

% We need to know where pumps and valves are on the FULL graph before any
% reductions
PumpEdges = [1;14];
ValveEdges = [4;11];

fooSim = HydraulicNetworkSimulation(fooGraph,Components,PumpEdges,ValveEdges,Inertias,NodeHeights,p0,rho,g);

%%

 % Collect the pressure equations in one vector (chords at the top)
                Omega_T = sym(zeros(length(fooGraph.spanT),1));
                Omega_C = sym(zeros(length(fooGraph.chords),1));
                
                for ii = 1:length(fooGraph.spanT)
                    Omega_T(ii) = fooSim.Omega_T(ii);
                end
                for ii = 1:length(fooGraph.chords)
                    Omega_C(ii) = fooSim.Omega_C(ii);
                end
                
                Omegas = [Omega_C;Omega_T];

                syms w1 w2 OD1 OD2 p q3 q7 d1 d5 d9 d8
                w = 66; OD = 0.5; tankflow = 0;
                
                
% Omegas = subs(Omegas,[w1 w2],[w w]);
% Omegas = subs(Omegas,[OD1 OD2],[OD OD]);
% Omegas = subs(Omegas,d6,tankflow);

ResistancePart = fooGraph.Phi*Omegas;
HeightPart = (fooGraph.Psi*(fooSim.NodeHeights))*(fooSim.rho*fooSim.g)/(10^5);
PressurePart = (fooGraph.I*(p-0));

dqdt = fooSim.P*(-ResistancePart+HeightPart+PressurePart);

jacobq = jacobian(dqdt,[q3 q7 d1 d5 d9 d8]);
jacobw = jacobian(dqdt,[w1 w2]);

dqdt = subs(dqdt,[w1 w2],[w w]);
dqdt = subs(dqdt,[OD1 OD2],[OD OD]);
dqdt = subs(dqdt,d8,tankflow);

eqpoint = solve(dqdt == 0);
q0 = struct2array(eqpoint);
pt = q0(4);
q0 = [q0(5:6) q0(1:3) 0]';

w = 66; w = 66;
OD = 0.5; OD = 0.5;

A = subs(jacobq,[OD1 OD2 q3 q7 d1 d5 d9 d8],[OD OD q0']);
B = subs(jacobw,symvar(jacobw),[w w]);

eig(A)

% Just to check that f(x0) is actually 0
[dqdt,pbar,pt_new] = fooSim.Model_TimeStep([w w],[OD OD],[q0(1) q0(2)],[q0(3) q0(4) q0(5) 0],q0(6),pt);
double(dqdt)

%% simulering
A = double(A);
B = double(B);
C = [0 0 -1 -1 -1 -1];

ts = 0.025;

LinSys = ss(A,B,C,[]);
LinSys = c2d(LinSys,ts,'Tustin');



f_p = 2/(24*10); % Assumed pump frequency

t = 0:ts:(2*60); % Time vector corresponding to 24 hours
% w = 50*(1+square(12*pi*f_p*t,50))/2+50*(1+square(12*pi*f_p*t,100))/2; % Pump waveforms
% OD =(1+sin(pi*f_p*t+t(end)/2))/2; % Valve waveforms


clear flow tankpres df d_t pt w1 w2 OD1 OD2 pressures q_lin pt_lin flowchange LinPump2
q0 = double(q0);
w1 = 66*ones(length(t),1); w2 = 66*ones(length(t),1);
% w1 = w; w2 = w; 
% OD1 = OD; OD2 = OD;
OD1 = 0.5*ones(length(t),1); OD2 = 0.5*ones(length(t),1);
qc = [0;0];
% qc = q0(1:2);
% df = [q0(3:end-1);0];
df = [0;0;0;0]; 
pt =0*4*rho*g/10^5;
pt_lin = 0;
% d_t = 0;
d_t = 0;

w0 = 66;

q_lin(:,1) = zeros(6,1);

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
    
    flowchange(:,ii) = dqdt;

    pt = pt_new;
    flow(:,ii) = [qc;df;d_t];
    flowL = [qc;df(1:3);d_t];
    pressures(:,ii) = double(pbar);
    tankpres(ii) = pt;
    
    q_lin(:,ii+1) = LinSys.A*(flowL-q0)+LinSys.B*([w1(ii);w2(ii)]-w0);
    pt_lin(ii+1) = pt_lin(ii)-0.000096*q_lin(end,ii);
end

% Plots

q_lin = q_lin+q0;

close all

figure()
subplot(2,2,1)
plot(t,flow(1:2,:))
hold on
plot(t,q_lin(1:2,1:end-1),'--')
hold off
legend('Chord 1','Chord 2','LinChord1','LinChord2')
subplot(2,2,2)
plot(t,flow(3,:))
hold on
plot(t,flow(6,:))
plot(t,q_lin(3,1:end-1),'--')
hold off
legend('Pump 1','Pump 2','LinPump1')
subplot(2,2,3)
plot(t,flow(4:5,:))
hold on
plot(t,q_lin(4:5,1:end-1),'--')
hold off
legend('Consumer 1','Consumer 2','LinConsumer1','LinConsumer2')
subplot(2,2,4)
plot(t,flow(end,:))
hold on
plot(t,q_lin(end,1:end-1),'--r')
hold off
legend('Tank','Linearized Tank')
figure()
plot(t,pressures)
figure()
subplot(3,1,1)
plot(t,tankpres)
hold on
plot(t,pt_lin(1:end-1),'--')
legend('Tank Pressure','Linearized Tank Pressure')
subplot(3,1,2)
plot(t,OD1)
hold on
plot(t,OD2)
hold off
legend('Valve 1','Valve 2')
subplot(3,1,3)
plot(t,w1)
hold on
plot(t,w2)
hold off
legend('Pump 1','Pump 2')


%%

% syms d1 d5 d8 d9 q3 q7 omega
% qn = [q3 q7 d1 d5 d9 0 d8]';




% omega_tilde = omega - w;



% qn_tilde = qn - q0;
% qn_lin_dot = (a1(1)*w+(abs(q0)+sign(q0)*q0')*(ResistancePart+a2(1)+1/(Kv*OD)^2)*qn_ntilde)-((a1(1)*q0+2*a0(1))*omega_tilde)
% %%

% ts = 0.00001;
% 
% 
% t = 0:ts:(6*60)/10; % Time vector corresponding to 6 min
% 
% clear flow tankpres df d_t pt w1 w2 OD1 OD2 pressures
% 
% w1 = 66*ones(length(t),1); w2 = 66*ones(length(t),1);
% OD1 = 0.5*ones(length(t),1); OD2 = 0.5*ones(length(t),1);
% 
% qc = [0;0];
% df = [0;0;0;0]; 
% pt = 4*rho*g/10^5; % 4 bar normalt i boliger
% d_t = 0;
% 
% for ii = 1:length(t)
%     [dqdt,pbar,pt_new] = fooSim.Model_TimeStep([w1(ii) w2(ii)],[OD1(ii) OD2(ii)],[qc(1) qc(2)],[df(1) df(2) df(3) df(4)],d_t,pt);
%     qc(1) = qc(1)+dqdt(1)*ts;
%     qc(2) = qc(2)+dqdt(2)*ts;
%     df(1) = df(1)+dqdt(3)*ts; 
%     df(2) = df(2)+dqdt(4)*ts; 
%     df(3) = df(3)+dqdt(5)*ts; 
%     d_t = d_t+dqdt(end)*ts;
% 
%     % Must obey basic mass conservation if no leaks are assumed
%     masscon = -cumsum([df(1:3);d_t]);
%     df(4) = masscon(end);
% 
%     pt = pt_new;
%     flow(:,ii) = [qc;df;d_t];
%     pressures(:,ii) = double(pbar);
%     tankpres(ii) = pt;
% 
% 	dqdt_saved(:,ii) = dqdt;
% 
% end
% 
% %%
% % 
% % W0 = [66;0;0;0;0;0;0;0;0;0;0;0;0;66];
% % a0_array = [a0(1);0;0;0;0;0;0;0;0;0;0;0;0;a0(2)];
% % a1_array = [a1(1);0;0;0;0;0;0;0;0;0;0;0;0;a1(2)];
% % a2_array = [a2(1);0;0;0;0;0;0;0;0;0;0;0;0;a2(2)];
% % 
% % syms d1 d5 d8 d9 q3 q7
% % % NB!!! tilføj cords!
% % q_0_no_chords = subs(fooSim.q_T,[d1 d5 d8 d9 q3 q7],[eqpoint.d1 eqpoint.d5 eqpoint.d8 eqpoint.d9 eqpoint.q3 eqpoint.q7]);
% % q_0 = [q_0_no_chords(1:2);eqpoint.q3;q_0_no_chords(3:6);eqpoint.q7;q_0_no_chords(7:12)];
% % OD0 =[0;0;0;0.5;0;0;0;0;0;0;0.5;0;0;0];
% % 
% % for i = 1:length(Pipes)
% %     K_lampda_only_pipes(i) = ((Pipes(i).hm+Pipes(i).hf)*Pipes(i).rho)/(10^5*3600^2);
% % end
% % 
% % K_lampda = [0,K_lampda_only_pipes(1:2),0,K_lampda_only_pipes(3:8),0,K_lampda_only_pipes(9:10),0];
% % Kv_array = [zeros(1,3),1,zeros(1,6),1,zeros(1,3)];
% % 
% % 
% % %dd_pipes = a1_array.*W0+(abs(q_0)+sign(q_0).*q_0)*(K_lampda'+a2_array+1./(Kv_array'.*OD0).^2)
% % %dd_pipes = (a1_array(2)*W0(1)+(abs(q_0(2))+sign(q_0(2)).*q_0(2))*(K_lampda(2)+a2_array(2)+1/(Kv_array(2)*OD0(2)).^2))
% % % dd_pump = a1*q0+2*a0*W0
% 
% 
% %%
% close all; clear ax
% 
% x      = 4;   % Screen position
% y      = 3;   % Screen position
% width  = 30; % Width of figure
% height = 25; % Height of figure (by default in pixels)
% 
% figures = []
% 
% legends = {	["Chord 1","Chord 2"];
% 			["Pump 1","Pump 2"];
% 			["Consumer 1","Consumer 2"];
% 			["Tank Flow"]}
% 
% tits = {	["Chord Flows (q_C)"];
% 			["Pump Flows (d_p)"];
% 			["Consumer Flows (d_c)"];
% 			["Tank flow (d_t)"]}
% 
% ylab = "Flow [m^3/h]"
% xlab = "Time [sec]"
% 
% flowindices = {	[1,2];	% chords
% 				[3,6];	% pumps
% 				[4,5];	% consumers
% 				[7]}	% tank
% 
% 
% 
% figures = [figure( 'Color', 'white', 'Units','centimeters','Position', [x y width height])]
% for i=1:4
% 	ax(i) = subplot(4,1,i)
% 	plot(t,flow(flowindices{i},:))
% 	legend(legends{i})
% 	xlabel(xlab)
% 	ylabel(ylab)
% 	title(tits{i})
% 	grid on
% end
% linkaxes(ax,'x')
% xlim([0 2.05])
% % sgtitle(['Flows in Water network; initial tank pressure = ', num2str(pt_init)], 'Fontsize',20)
% 
% 
% %% tank, flow, OD/omega
% clear ax
% output = {	tankpres; 
% 			flow(7,:);
% 			[w1, w2]';
% 			[OD1, OD2]';}
% 
% legends = {	["Tank Pressure"];
% 			["Tank Flow"];
% 			["Pump1 Speed","Pump1 Speed"];
% 			["Consumer 1","Consumer 2"]}
% 
% tits = {	["Tank Pressure"];
% 			["Tank Flow"];
% 			["Inputs"];
% 			["Disturbances"]}
% 
% ylab = {	["Tank Pressure [bar]"];
% 			["Flow [m^3/h]"];
% 			["Pump speed [%]"];
% 			["Opening degree"]}
% xlab = "Time [sec]"
% ylims = {	[-10 10];
% 			[-10 10];
% 			[-0.5 100.5];
% 			[-0.005 1.005]}
% 
% 
% figures = [figures; figure( 'Color', 'white', 'Units','centimeters','Position', [x y width height])]
% for i=1:4
% 	ax(i) = subplot(4,1,i)
% 	plot(t,output{i})
% 	legend(legends{i})
% 	xlabel(xlab)
% 	ylabel(ylab{i})
% 	title(tits{i})
% 	ylim(ylims{i})
% 	grid on
% end
% linkaxes(ax,'x')
% xlim([0 2.05])
% % sgtitle(['Tank measures affected by inputs and disturbances; initial tank pressure = ', num2str(pt_init)], 'Fontsize',20)
% 
% %% 
% % idea afterwards: identify which pipe has bigger inertia ie. pressure
% % drop in each pipe and flow in each pipe
% clear ax
% %				H_bar' 11x9 * 9x1 = 
% delta_p = fooSim.Graph.H_bar'*pressures;
% % manual way - should be doable automatically
% flowsall = zeros([length(fooSim.Graph.edges),length(t)]);
% flowsall(fooSim.Graph.chords,:) = flow(1:2,:); % chords
% 
% output = {	[delta_p];
% 			[flowsall]} % pressures and 
% legends = {}
% 
% for i = 1:length(fooSim.Graph.edges)
% 	legends{i,1} = [sprintf("Pressure drop: edge %d", fooSim.Graph.edges(i)), ...
% 		sprintf("Flow: edge %d", fooSim.Graph.edges(i))]
% end
% 
% 
% % tits = {	["Tank Pressure"];
% % 			["Tank Flow"];
% % 			["Inputs"];
% % 			["Disturbances"]}
% % 
% % 
% ylab = "Pressure [bar]/ flow[m^3/h]"
% xlab = "Time [sec]"
% % ylims = {	[-inf inf];
% % 			[-inf inf];
% % 			[-0.5 100.5];
% % 			[-0.005 1.005]}
% 
% 
% figures = [figures; figure( 'Color', 'white', 'Units','centimeters','Position', [x y width height])]
% 
% for i=1:11
% 	ax(i) = subplot(4,3,i)
% 	plot(t,output{1,1}(i,:))
% 	hold on
% 	plot(t,output{2,1}(i,:))
% 	legend(legends{i})
% 	xlabel(xlab)
% 	ylabel(ylab)
% % 	title(tits{i})
% % 	ylim(ylims{i})
% 	grid on
% end
% % linkaxes(ax,'x')
% xlim([-2 t(end)+1])
% sgtitle(['Pressure and flow in edges', ], 'Fontsize',20)
% 
% 
% 
% %%
% % close all
% % 
% % figure()
% % subplot(2,2,1)
% % plot(t,flow(1:2,:))
% % legend('Chord 1','Chord 2')
% % subplot(2,2,2)
% % plot(t,flow(3,:))
% % hold on
% % plot(t,flow(6,:))
% % hold off
% % legend('Pump 1','Pump 2')
% % subplot(2,2,3)
% % plot(t,flow(4:5,:))
% % legend('Consumer 1','Consumer 2')
% % subplot(2,2,4)
% % plot(t,flow(end,:))
% % legend('Tank')
% % 
% % figure()
% % plot(t,pressures)
% % 
% % figure()
% % subplot(3,1,1)
% % plot(t,tankpres)
% % legend('Tank Pressure')
% % subplot(3,1,2)
% % plot(t,OD1)
% % hold on
% % plot(t,OD2)
% % hold off
% % legend('Valve 1','Valve 2')
% % subplot(3,1,3)
% % plot(t,w1)
% % hold on
% % plot(t,w2)
% % hold off
% % legend('Pump 1','Pump 2')
% % 
