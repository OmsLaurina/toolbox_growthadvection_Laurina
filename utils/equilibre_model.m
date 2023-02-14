clear all, close all

% Calcul equilibre du modèle simplifié (sans la limitation du consommateur)

%Parametres
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

%Valeurs calculées : 
Z_barre = P2*gmax2*(gamma2-eZ+eZ*gamma2)/mZ; 
PO4_barre = (Psupply+epsilon*eZ*P2*gmax2*(1-gamma2))/(umax2*P2);
P2_barre = Psupply/(PO4*umax2+epsilon*eZ*gmax2*(gamma2-1));

% Matrice jacobienne

j11 = PO4_barre*umax1-gmax1*Z_barre;
j12 = 0;
j13 = 0;
j14 = 0;

j21 = 0;
j22 = PO4_barre*umax2-gmax2*Z_barre;
j23 = PO4_barre*umax2*P2_barre-gmax2*P2_barre;
j24 = umax2*P2_barre-P2_barre*gmax2*Z_barre;

j31 = gamma1*gmax1*gamma2*P2_barre*gmax2*Z_barre-eZ*((1-gamma1)*gmax1*Z_barre+(1-gamma2)*P2_barre*gmax2*Z_barre)-mZ*Z_barre^2;
j32 = gamma2*gmax2-eZ*(1-gamma2)*gmax2-mZ*Z_barre^2;
j33 = P2_barre*gmax2*gamma2-eZ*gamma2*gmax2-2*mZ*Z_barre;
j34 = 0;

j41 = Psupply+epsilon*eZ*((1-gamma1)*gmax1*Z_barre+(1-gamma2)*P2_barre*gmax2*Z_barre)-PO4_barre*umax1+PO4_barre*P2_barre*umax2;
j42 = Psupply+epsilon*eZ*((1-gamma2)*gmax2*Z_barre)-PO4_barre*umax2;
j43 = Psupply+epsilon*eZ*((1-gamma2)*P2_barre*gmax2)-PO4_barre*umax2*P2_barre;
j44 = Psupply+epsilon*eZ*((1-gamma2)*P2_barre*gmax2*Z_barre)-umax2*P2_barre;

J_barre1 = [j11 j12 j13 j14; j21 j22 j23 j24; j31 j32 j33 j34; j41 j42 j43 j44];

%Determinant matrice non carré 


%% Pour P2 = 0

%Valeurs simulées :
% PO4 = 0.4281;
% P1 = ?;
% Z = 1.1886;
% 
% Z_barre = (PO4/gmax1)*umax1;
% PO4_barre = (Psupply+epsilon*eZ*(P1*gmax1-gamma1*P1*gmax1))/(umax1*P1);
% P1_barre = Psupply/(PO4*umax1-epsilon*eZ*gmax1+epsilon*eZ*gamma1*gmax1);









