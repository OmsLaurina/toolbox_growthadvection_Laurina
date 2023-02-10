%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       growth model -1D TEST (temporal evolution only)      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Starting script 1D (time only) of different versions of    %%%
%%% adapted model                                              %%%
%%% Use function : growth_model_*.m                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Author : Laurina Oms                                       %%%
%%% Creation : 14/11/2022                                      %%%     
%%% Region : Mediteranean Sea                                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

clear all, close all

n_days = 2000; % nb de jours de simulation
dt = 0.2; % pas de temps
time=(0:dt:n_days); % temps

% Psupply (mmolC m^{-3} d^{-1}) = apport externe de nutriment
Psupply_moy = 0.01; %0.06 = seuil a partir duquel Pbig prend le dessus et gagne

%Psupply constant
Psupply = ones(length(time),1)*Psupply_moy;

% %Psupply variable
% w = 0.5; %periode (day)
% Psupply_sin = time*NaN;
% Psupply = time*NaN;
% 
% %Sinusoïdalement
% b = 0.01; %amplitude (mmolC/m3)
% for i=1:length(time)
%     Psupply_sin(i)=Psupply_moy+b*sin(w*i);
%     Psupply(i) = Psupply_sin(i);
% end

% %Pulsé
% T = 4; %(pulses en jours)
% b = (Psupply_moy*0.25*T)/(1-exp(-0.25*T))-Psupply_moy;
% for i=1:length(time)
%     Psupply_sin(i)=Psupply_moy+b*sin(w*i);
%     Psupply(i) = Psupply_sin(i);
% end

% Concentration initiale (mmolC m^{-3})
C_nut = 0.5;

%% --------------- Adapted model version 7
% output_1 = growth_model_2P2Z_v7(Psupply,C_nut,Pini1,Pini2,'kP_small',1,'kP_big',2,'gmax_small',1.4226,'gmax_big',1.1226,'time',time,'plot');

%% --------------- Adapted model version 8
%  output_2 = growth_model_2P1Z_v8(Psupply,C_nut,'time',time,'plot');

%% --------------- Adapted model version 9 
output_3 = growth_model_2P1Z_v9(Psupply,C_nut,'time',time,'plot');

toc
