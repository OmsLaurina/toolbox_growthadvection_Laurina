function varargout=ga_model_2P1Z_v8(Nsupply,C_nut,Pini1, Pini2,varargin)
%% GA_MODEL_2P1Z_FROMPSUPPLY
% Default parameters are based on Messié & Chavez (2017).
% 
% Use:
% [output=]ga_model_2P1Z_*(Nsupply,varargin)
%
% output is a structure containing 
%	.time, .Nsupply, .P_1, .P_2, .Z_1, .Z_2, 
%	.Nnew, .Nreg, .Chl, .PP, .u_1, .u_2, .g_1, .g_2, 
%	.units, .attributs
%
% Required input:
% 	Nsupply expressed in mmolC/m3/d, all units are carbon-based. Note - Nnew and Nreg represent new and regenerated PO4rients, respectively
%		(NO3 and NH4 as a simplification) but are termed "Nnew" and "Nreg" to limit confusion with the unit, since they are expressed in carbon.
% 	Nsupply can be either a number, corresponding to the rate observed during upw_duration (default 1 day following Messié & Chavez 2017)
%					or a vector (then time needs to be given), that provides Nsupply as a function of time along a current trajectory,
%					for instance (useful to take Ekman pumping into account)
%
% Optional inputs:
% 'nbdays_advec'	number of days during which the model is run
% 'dt'				time step
% 'time'			can replace nbdays_advec and dt (required if Nsupply is a vector)
% 'upw_duration'	number of days during which Nsupply happens (default 1 day, not used if Nsupply is a vector)
% 'plot'			displays the plankton model outputs as a function of time
%
% Monique Messié, 2021 for public version
% Reference: Messié, M., & Chavez, F. P. (2017). PO4rient supply, surface currents, and plankton dynamics predict zooplankton hotspots 
%					in coastal upwelling systems. Geophysical Research Letters, 44(17), 8979-8986, https://doi.org/10.1002/2017GL074322
% Differences with Messié and Chavez (2017):
%		all Z1 excretion is now availabe as regenerated nutrients (ie no export on Z1 excretion)
%		Z2 grazing formulation is different with the half-saturation constant applying to Z_1+P_2 in both cases


%% -------------- Default parameters (see Messié & Chavez, 2017 suppl inf)

% 1 = Psmall
% 2 = Pbig

default_parameters={...
'umax_1',1.9872,'umax_2',2.7648,'gmax1',1.4226,'gmax2',1.1126,...% maximum growth and grazing rates (d^{-1}) (Baklouti et al., 2021)
'cChl_1',200,'cChl_2',50,...								% C:Chl ratios for P1 and P2 (only used to calculate Chl) (Lazzari et al., 2012)
'kP_1',1,...											% half-saturation constant for P1 on PO4 (mmolC m^{-3}) (Pulido-Villena et al., 2021) )
'kP_2',2,...												% half-saturation constant for P2 on PO4 (mmolC m^{-3}) (Pulido-Villena et al., 2021)
'kG',5,...										% half-saturation constant for Z on P (mmolC m^{-3}) (Auger et al., 2011)
'mP',0,...														% P2 mortality rate (default 0 ie no P2 sinking) (d^{-1}) (Auger et al., 2011)
'mZ',0.005,...											% Z2 quadratic mortality rate (mmolC^{-1} m^{3} d^{-1})
'eZ',0.1,...													% zoo excretion fraction (Z1 and Z2) (d^{-1}) (Baklouti et al., 2021)
'epsilon',0.75 ,...												% fraction of Z excretion that is available as regenerated PO4rients
'P_1_ini',Pini1,'P_2_ini',Pini2,'Z_ini',0.6};	% initial biomass (mmolC m^{-3}) (In situ for P_1_ini with method from Marrec et al., 2018 and Tzortzis et al., 2021 & further)

[arg,flag]=ga_read_varargin(varargin,[{'nbdays_advec',20,'dt',0.2,'time',[]},default_parameters],{'plot'});
if length(Nsupply)>1 && isempty(arg.time), error('Give time if Nsupply is a vector'), end


%% -------------- Time

if isempty(arg.time), time=(0:arg.dt:arg.nbdays_advec)'; 
else, time=arg.time(:); arg.dt=time(2)-time(1); 
end
nb_time=length(time); 



%% -------------- Nsupply

if length(Nsupply)==1
	Nsupply=zeros(nb_time,1);
end

%% -------------- Initial conditions

%Biomasses
PO4=time*NaN; 				% PO4rients (PO4 expressed in carbon) (mmolC m^{-3})
P_1=time*NaN; 			% 1 phyto biomass (mmolC m^{-3})
P_2=time*NaN; 			% large phyto biomass (mmolC m^{-3})
Z=time*NaN; 			% 1 zoo biomass (mmolC m^{-3})
u_1=time*NaN; 			% P_1 growth rate (Nreg-limited) (d^{-1})
u_2=time*NaN; 	        % P_2 growth rate (Nnew-limited) (d^{-1})
g=time*NaN; 			% Z_2 grazing rate (P_2 and Z_1 limited) (d^{-1})
PP_1=time*NaN;         % primary production for P_1 (mmolC m^{-3} d^{-1})  
PP_2=time*NaN; 			% primary production for P_2 (mmolC m^{-3} d^{-1})

%Flux
G1=time*NaN; 			% grazing from Z onto P2 (mmolC m^{-3} d^{-1})
G2=time*NaN; 			% grazing from Z onto P1 (mmolC m^{-3} d^{-1})
excretion_Z=time*NaN; 	% Z excretion (mmolC m^{-3} d^{-1})
death_Z=time*NaN; 		% Z2 mortality term (mmolC m^{-3} d^{-1})
regeneration_PO4=time*NaN;	% PO4rient regeneration (mmolC m^{-3} d^{-1})
p_fec=time*NaN;
Export=time*NaN;
Budget=time*NaN;

%Valeurs initiales
P_1(1)=arg.P_1_ini; 
P_2(1)=arg.P_2_ini; 
Z(1)=arg.Z_ini;
PO4(1)=C_nut;

%% -------------- Loop on time

for t=2:nb_time

	% growth and grazing rates
    u_1=PO4(t-1)/(arg.kP_1+PO4(t-1))*arg.umax_1;
    u_2=PO4(t-1)/(arg.kP_2+PO4(t-1))*arg.umax_2;
  
    g1=P_1(t-1)/(arg.kG+P_1(t-1)+P_2(t-1))*arg.gmax1; 	
	g2=P_2(t-1)/(arg.kG+P_1(t-1)+P_2(t-1))*arg.gmax2;
	g(t-1)=g1+g2;

	% fluxes
    PP_1(t)=u_1*P_1(t-1); 
    PP_2(t)=u_2*P_2(t-1);
    
	G1(t)=g1*Z(t-1); 
	G2(t)=g2*Z(t-1); 
	death_Z(t)=arg.mZ*Z(t-1)^2;
	excretion_Z(t)=arg.eZ*Z(t-1);
	regeneration_PO4(t)=arg.epsilon*excretion_Z(t);
    p_fec(t) = (1-arg.epsilon)*excretion_Z(t);%-regeneration_PO4(t);
    
    max_available_PO4 = PO4(t-1)+Nsupply(t)*arg.dt+regeneration_PO4(t)*arg.dt;
	if PP_1(t)>max_available_PO4/arg.dt, PP_1(t)=max_available_PO4/arg.dt; end
    max_available_PO4=max_available_PO4-PP_1(t)*arg.dt;
	if PP_2(t)>max_available_PO4/arg.dt, PP_2(t)=max_available_PO4/arg.dt; end
    
    % export 
    Export(t) =p_fec(t)+death_Z(t); 
    
	% carbon-based nutrients and biomass
	PO4(t)=PO4(t-1)+Nsupply(t)*arg.dt+regeneration_PO4(t)*arg.dt-PP_1(t)*arg.dt-PP_2(t)*arg.dt; PO4(PO4<=0)=0;
    
	P_1(t)=P_1(t-1)+PP_1(t)*arg.dt-G1(t)*arg.dt; P_1(P_1<=0)=0;
    P_2(t)=P_2(t-1)+PP_2(t)*arg.dt-G2(t)*arg.dt ;P_2(P_2<=0)=0;
    Z(t)=Z(t-1)+G1(t)*arg.dt+G2(t)*arg.dt-excretion_Z(t)*arg.dt-death_Z(t)*arg.dt; Z(Z<=0)=0;
    
end

% save('outputs/maxPP_21', 'maxPP_21');

%% -------------- Ouputs

units=struct('time','days','Nsupply','mmolC m^{-3} d^{-1}',...
	'P_1','mmolC m^{-3}','P_2','mmolC m^{-3}','Z','mmolC m^{-3}',...
	'PO4','mmolC m^{-3}','Chl','mg m^{-3}','PP','gC m^{-3}/yr',...
	'u_1','d^{-1}','u_2','d^{-1}','g_1','d^{-1}','g_2','d^{-1}', 'PP_1', 'gC m^{-3}/yr', 'PP_2','gC m^{-3}/yr', 'Export', 'mmolC m^{-3} d^{-1}','Budget','mmolC m^{-3} d^{-1}','death_Z', 'mmolC m^{-3} d^{-1}', 'p_fec','mmolC m^{-3} d^{-1}');

output=struct('units',units,'time',time,'Nsupply',Nsupply,...
	'P_1',P_1,'P_2',P_2,'Z', Z,'PO4',PO4,...
	'Chl',P_1*12/arg.cChl_1+P_2*12./arg.cChl_2,'PP',(PP_2+PP_1)*1E-3*12*365.25,...
	'u_1',u_1,'u_2',u_2,'g',g, 'PP_1', PP_1,'PP_2', PP_2, 'Export', Export,'Budget', Budget,'death_Z',death_Z, 'p_fec', p_fec, 'attributs',struct('arg',arg));
varargout={output}; varargout=varargout(1:nargout);

%% -------------- Figures

if flag.plot

    %Temporal evolution
	figure, hold on
	plot(output.time,output.P_1,'LineWidth',2, 'Color','g')
	plot(output.time,output.P_2,'LineWidth',2,'Color', '#77AC30')
	plot(output.time,output.Z,'LineWidth',2,'Color','#00FFFF')
    plot(output.time,output.PO4,'LineWidth',2, 'Color','m')
	ylabel(output.units.P_1)
	if min(output.time)>datenum(1900,1,1), datetick('x','keeplimits'), xlabel('Time')
	else, xlabel('Time (days)')
	end
	legend({'P\_small','P\_big','Z','PO4'})
	title('Model output (plankon concentration over time)')
    xlim([min(output.time) max(output.time)]);
    ylim([0 1.2]);
    
    %Temporal evolution budget
	figure, hold on
	plot(output.time,output.Export,'LineWidth',2, 'Color','k')
    plot(output.time,output.Nsupply,'LineWidth',2, 'Color','r')
	if min(output.time)>datenum(1900,1,1), datetick('x','keeplimits'), xlabel('Time')
	else, xlabel('Time (days)')
	end
	legend({'Export', 'Psupply'})
	title('Model output (Budget over time)')
    xlim([min(output.time) max(output.time)]);

    %Monod equation
    figure, hold on
    N_theo_PO4 = linspace(0,10,length(arg.time));
    f_monod2(N_theo_PO4,arg.umax_1,arg.kP_1);
    f_monod2(N_theo_PO4,arg.umax_2,arg.kP_2);
    f_monod2(output.PO4,arg.umax_1,arg.kP_1);
    f_monod2(output.PO4,arg.umax_2,arg.kP_2);
    xlabel('[PO4] mmolC.m³')
    ylabel('\mu d^{-1}');
    legend({'P\_1\_theo','P\_2\_theo','P\_1\_mod','P\_2\_mod' })
    
    %Portrait de phase
    figure, hold on
    % Calcul des dérivées de P_1 et P_2 par rapport au temps
    dp1_dt = gradient(output.P_1, arg.dt);
    dp2_dt = gradient(output.P_2, arg.dt);
    
    % Tracé du portrait de phase
    plot(output.P_2, output.P_1)
    xlabel('P\_big')
    ylabel('P\_small')
    
    % Tracé du champ de vecteur sur le portrait de phase
    quiver(output.P_2, output.P_1, dp2_dt, dp1_dt)
     
end

return