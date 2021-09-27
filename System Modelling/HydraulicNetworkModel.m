% Incidence matrix

H = [1 0 1 0 0 0 0; -1 1 0 0 0 0 0; 0 0 -1 1 1 0 0; 0 -1 0 -1 0 1 0;...
    0 0 0 0 -1 0 -1; 0 0 0 0 0 -1 1];

% Define reference point and form reduced incidence matrix

vref = 4; 

Hbar = H;
Hbar(4,:) = []

n = size(H,1);

% Chords are e1 and e4. Spanning tree consists of remaining edges.
% Define partitioned matrices

Hc = [Hbar(:,1),Hbar(:,4)]

nc = size(Hc,2); % Number of chords

Ht = [Hbar(:,2) Hbar(:,3) Hbar(:,5) Hbar(:,6) Hbar(:,7)]

% Construct loop matrix from partitioned matrices

B = [eye(nc,nc) -(Hc'*inv(Ht)')]