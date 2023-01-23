%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       growth model -1D TEST (temporal evolution only)      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Use function : ga_model_*.m                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Author : Laurina Oms                                       %%%
%%% Creation : 14/11/2022                                      %%%     
%%% Region : Mediteranean Sea                                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
clear all, close all

n_days = 600; %nb de jours de simulation
dt = 0.2; %pas de temps
time=(0:dt:n_days); %temps

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

% Concentration initiale de phosphate (mmolC m^{-3})
C_nut = 0.5;
Pini1 = 0.6;
Pini2 = 0.1;

%% --------------- Initial model (without changing)
% use script ga_model_2P2Z_v1.m 
%  output_1 = ga_model_2P2Z_v1(Psupply,C_nut,'kNreg',13.3,'kNnew',13.3,'time',time,'plot');

%% --------------- Improve model to the Med sea = Temperature dependent model (without changing )
% % use script ga_model_2P2Z_v2.m
% output_2 = ga_model_2P2Z_v2(Psupply,SSTsupply,Nini,psmall_ini,pbig_ini,'kNreg',0.5/16*106,'kNnew',0.75/16*106,'time',time,'plot');

%% --------------- Initial model (with changing : usmall = f(Nnew,Nreg))
% % use script ga_model_2P2Z_v3.m  
% output_3 = ga_model_2P2Z_v3(Psupply,C_nut,'kNreg',13.3,'kNnew',13.3,'time',time,'plot');
% 
%% --------------- Initial model (with changing : usmall= f(Nnew+Nreg))
% use script ga_model_2P2Z_v4.m  
% output_4 = ga_model_2P2Z_v4(Psupply,C_nut,'kNreg',13.3,'kNnew',13.3,'time',time,'plot');

%% --------------- TEST : Initial model (with changing PPsmall = f(Nnew) NO Nreg)
% % use script ga_model_2P2Z_v5.m  
% output_5 = ga_model_2P2Z_v5(Psupply,C_nut,'kNreg',13.3,'kNnew',13.3,'time',time,'plot');
% 
%% --------------- Initial model (with changing : usmall,u_big= f(Nnew+Nreg))
% % use script ga_model_2P2Z_v6.m  
% output_6 = ga_model_2P2Z_v6(Psupply,C_nut,'kNreg',13.3,'kNnew',13.3,'time',time,'plot');

%% --------------- Initial model (with changing : usmall,u_big= f(PO4)
% use script ga_model_2P2Z_v7.m  
output_7 = ga_model_2P2Z_v7(Psupply,C_nut,Pini1,Pini2,'kP_small',1,'kP_big',2,'gmax_small',1.4226,'gmax_big',1.1226,'time',time,'plot');

%% --------------- Initial model (with changing : 1Z)
% use script ga_model_2P1Z_v7.m  
%  output_8 = ga_model_2P1Z_v8(Psupply,C_nut,'time',time,'plot');

%% --------------- pcolor ratio  
% output_9 = ga_model_pcolor_param_u_KP_v8(Psupply,C_nut,'time',time,'plot');

%% --------------- pcolor ratio
% output_10 = ga_model_pcolor_param_gmax_v8(Psupply,C_nut,'time',time,'plot');

toc
