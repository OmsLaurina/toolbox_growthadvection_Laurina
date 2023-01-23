function varargout=Berhenfeld_model_v1(n_days,dt,varargin)
%% Berhenfeld model adapted 2P1Z
% Default parameters ... ?)
% output is a structure containing 
%	.time, .P_small, .P_big, .Z, .u_small, .u_big, .units, .attributs
%
% Optional inputs:
% 'nbdays_advec'	number of days during which the model is run
% 'dt'				time step
% 'time'			can replace nbdays_advec and dt (required if Nsupply is a vector
% 'plot'			displays the plankton model outputs as a function of time

%% -------------- Default parameters and initial conditions

default_parameters={...
'c1_small',0.1,... % herbivore grazing rate for P_small
'c1_big',0.2,... % herbivore grazing rate for P_big
'c2',0.5,... % ingestion efficiency
'c3',0.005,... % predatory loss rate
'P_small_ini',0.6,... % initial biomass (mmolC m^{-3})
'P_big_ini',0.1,... 
'Z_ini',0.3,...
'u_small_ini', 1.9872,... % initial growth rate (d^-1)
'u_big_ini', 2.7648};

[arg,flag]=ga_read_varargin(varargin,[{'nbdays_advec',n_days,'dt',dt,'time',[]},default_parameters],{'plot'});

%% -------------- Time

if isempty(arg.time), time=(0:arg.dt:arg.nbdays_advec)'; 
else, time=arg.time(:); arg.dt=time(2)-time(1); 
end
nb_time=length(time); 

%% -------------- Initial conditions

P_small=time*NaN; 			% small phyto biomass (mmolC m^{-3})
P_big=time*NaN; 			% large phyto biomass (mmolC m^{-3})
Z=time*NaN; 			    % zoo biomass (mmolC m^{-3})
u_big=time*NaN; 	        % P_big growth rate (Nnew-limited) (d^{-1})
u_small=time*NaN; 			% P_small growth rate (Nreg-limited) (d^{-1})

% Initial values
P_small(1)=arg.P_small_ini; 
P_big(1)=arg.P_big_ini; 
Z(1)=arg.Z_ini; 	
u_small(1)=arg.u_small_ini;
u_big(1) = arg.u_big_ini;

%lag-time between division and loss rates
j = 1;

%% -------------- Loop on time

for t=2:nb_time

    %Phytoplankton biomass
    P_small(t) = P_small(t-1) + u_small(t-1)*arg.dt - (arg.c1_small*P_small(t-1)*Z(t-1))*arg.dt;
    P_big(t) = P_big(t-1) + u_big(t-1)*arg.dt - (arg.c1_big*P_big(t-1)*Z(t-1))*arg.dt;

    %Zooplankton biomass
    Z(t) = Z(t-1) + (arg.c1_small*arg.c2*P_small(t-1))*arg.dt + (arg.c1_big*arg.c2*P_big(t-1))*arg.dt - (arg.c3*Z(t-1)^2)*arg.dt;
        
    %Growth rate calculated with j (ie: lag-time between division and loss rates)
    if mod(t,j)==0
        u_small(t) = u_small(t-j) + (1/P_small(t)) * (P_small(t)-P_small(t-j))/arg.dt;
        u_big(t) = u_big(t-j) + (1/P_big(t)) * (P_big(t)-P_big(t-j))/arg.dt;
    end
end

%% -------------- Ouputs

units=struct('time','days',...
	'P_small','mmolC m^{-3}','P_big','mmolC m^{-3}','Z','mmolC m^{-3}',...
	'u_small','d^{-1}','u_big','d^{-1}');
output=struct('units',units,'time',time,...
	'P_small',P_small,'P_big',P_big,'Z',Z,...
	'u_small',u_small,'u_big',u_big, 'attributs',struct('arg',arg));
varargout={output}; varargout=varargout(1:nargout);

%% -------------- Figures

if flag.plot

    %Temporal evolution of biomass
	figure, hold on
	plot(output.time,output.P_small,'LineWidth',2)
	plot(output.time,output.P_big,'LineWidth',2)
	plot(output.time,output.Z,'LineWidth',2)
	legend({'P\_small','P\_big','Z'})
	title('Model output (plankon concentration over time)')
    xlim([min(output.time) max(output.time)]); 

    %Temporal evolution of growth rate
	figure, hold on
	plot(output.time,output.u_small,'LineWidth',2)
	plot(output.time,output.u_big,'LineWidth',2)
	legend({'u\_small','u\_big'})
	title('Model output (plankon concentration over time)')
    xlim([min(output.time) max(output.time)]);
end

return