%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       growth model -1D TEST (temporal evolution only)      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Use function : ga_model_2P2Z_fromPsupply*.m                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Author : Laurina Oms                                       %%%
%%% Creation : 14/11/2022                                      %%%     
%%% Region : Mediteranean Sea                                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all

n_days = 5000;
time = datenum(0,1,1):datenum(0,1,1)+n_days;

% parameters from Nathan's master's thesis
Psupply = ones(n_days+1,1)*10E-3;
%SSTsupply = ones(n_days+1,1)+14;
SSTsupply= 13 + (15-13) .* rand(n_days+1,1); %T° variable
C_nut = 0.1;
Nini = 1;%0.6; 
psmall_ini = 0.6;
pbig_ini = 0.1;

%% --------------- TEST1 = initial model without changing 
% use script ga_model_2P2Z_fromPsupply_Med_1D.m   
output_1 = ga_model_2P2Z_fromPsupply_Med_1D(Psupply,C_nut,'time',time,'plot');

%% --------------- TEST2 = initial model with changing (PPsmall = f(Nnew))
% use script ga_model_2P2Z_fromPsupply_Med_PPsmall_fct_Nnew.m  
output_2 = ga_model_2P2Z_fromPsupply_Med_PPsmall_fct_Nnew(Psupply,C_nut,'time',time,'plot');

%% --------------- TEST3 = Temperature dependent model without changing 
% use script ga_model_2P2Z_fromPsupply_SST_Biss_Med.m
output_3 = ga_model_2P2Z_fromPsupply_SST_Biss_Med(Psupply,SSTsupply,Nini,psmall_ini,pbig_ini,'time',time,'plot');

%% --------------- TEST4 = Temperature dependent model with  changing (PPsmall = f(Nnew,Nreg))
% use script ga_model_2P2Z_fromPsupply_SST_Biss_Med_PPsmall_fct_Nnew.m
output_4 = ga_model_2P2Z_fromPsupply_SST_Biss_Med_PPsmall_fct_Nnew(Psupply,SSTsupply,Nini,psmall_ini,pbig_ini,'time',time,'plot');

%% --------------- TEST5 = initial model with changing (new usmall equation)
% use script ga_model_2P2Z_fromPsupply_Med_PPsmall_fct_Nnew_umodify.m  
output_5 = ga_model_2P2Z_fromPsupply_Med_PPsmall_fct_Nnew_umodify(Psupply,C_nut,'time',time,'plot');

%% --------------- TEST6 = Initial model with changing (PPsmall = f(Nnew) (without Nreg))
% use script ga_model_2P2Z_fromPsupply_Med_PPsmall_fct_Nnew_umodify.m  
output_5 = ga_model_2P2Z_fromPsupply_Med_PPsmall_fct_Nnew_only_umodify(Psupply,C_nut,'time',time,'plot');



