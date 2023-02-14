clear all, close all

% Calcul equilibre du modèle simplifié (sans la limitation du consommateur)

umax1=1.9872;   % maximum growth rates of P1 (d^{-1})
umax2=2.7648;   % maximum growth rates of P2 (d^{-1})
gmax1=1.4226;   % maximum grazing rates of Z on P1 (d^{-1})
gmax2=1.4226;   % maximum grazing rates of Z on P2 (d^{-1})
mZ=0.007;		% Z2 quadratic mortality rate (mmolC^{-1} m^{3} d^{-1})
eZ=0.1;	        % zoo excretion rate (Z1 and Z2) (d^{-1})
gamma1=0.7;     % conversion factor from P1 to Z
gamma2=0.7;     % conversion factor from P2 to Z
epsilon=0.75 ;  % fraction of Z excretion that is available as regenerated PO4
Psupply = 0.01;


%% Pour P1 = 0

%Valeurs simulées :
PO4 = 0.4428;
P2 = 0.009;
Z = 1.1886;

Z_barre = P2*gmax2*(gamma2-eZ+eZ*gamma2)/mZ;
PO4_barre = (Psupply+epsilon*eZ*P2*gmax2*(1-gamma2))/(umax2*P2);
P2_barre = Psupply/(PO4*umax2+epsilon*eZ*gmax2*(gamma2-1));

%% Pour P2 = 0

%Valeurs simulées :
% PO4 = 0.4281;
% P1 = ?;
% Z = 1.1886;
% 
% Z_barre = (PO4/gmax1)*umax1;
% PO4_barre = (Psupply+epsilon*eZ*(P1*gmax1-gamma1*P1*gmax1))/(umax1*P1);
% P1_barre = Psupply/(PO4*umax1-epsilon*eZ*gmax1+epsilon*eZ*gamma1*gmax1);









