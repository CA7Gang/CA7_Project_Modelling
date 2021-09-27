Length = 10; rho = 1; Area = 0.5; Diameter = 1.5; g = 9.82; eta = 0.03; Reynolds = 125000; kf = 1.8; dz = -3;

dP = 1; q = 1; % Units are bar and m^3/h

dPSI = dP*10e5; qSI = q*1/3600;  % Units are pascal and m^3/s 

fooPipe = PipeComponent(Length,rho,Area,Diameter,g,eta,Reynolds,kf,dz)

dQ = calcdQ(fooPipe,q,dP)

dQSI = calcdQSI(fooPipe,qSI,dPSI)

dQSI/dQ % Should correspond exactly to 1/3600
