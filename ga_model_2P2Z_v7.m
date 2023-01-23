function varargout=ga_model_2P2Z_v7(Nsupply,C_nut,Pini1, Pini2,varargin)
%% GA_MODEL_2P2Z_FROMNSUPPLY: plankton model used to model zooplankton hotspots
% Default parameters are based on Messié & Chavez (2017).
% 
% Use:
% [output=]ga_model_2P2Z_fromNsupply(Nsupply,varargin)
%
% output is a structure containing 
%	.time, .Nsupply, .P_small, .P_big, .Z_small, .Z_big, 
%	.Nnew, .Nreg, .Chl, .PP, .u_small, .u_big, .g_small, .g_big, 
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
%		all Zsmall excretion is now availabe as regenerated nutrients (ie no export on Zsmall excretion)
%		Zbig grazing formulation is different with the half-saturation constant applying to Z_small+P_big in both cases


%% -------------- Default parameters (see Messié & Chavez, 2017 suppl inf)


default_parameters={...
'umax_small',1.9872,'umax_big',2.7648,'gmax_small',1.4226,'gmax_big',1.1120,...% maximum growth and grazing rates (d^{-1}) (Baklouti et al., 2021)
'cChl_small',200,'cChl_big',50,...								% C:Chl ratios for Psmall and Pbig (only used to calculate Chl) (Lazzari et al., 2012)
'kP_small',13.3,...											% half-saturation constant for Psmall on PO4 (mmolC m^{-3}) (Pulido-Villena et al., 2021) )
'kP_big',15.3,...												% half-saturation constant for Pbig on PO4 (mmolC m^{-3}) (Pulido-Villena et al., 2021)
'kG_small',5,...											% half-saturation constant for Zsmall on Psmall (mmolC m^{-3}) (Auger et al., 2011)
'kG_big',5,...										% half-saturation constant for Zbig on Pbig and Zsmall (mmolC m^{-3}) (Auger et al., 2011)
'mP',0,...														% Pbig mortality rate (default 0 ie no Pbig sinking) (d^{-1}) (Auger et al., 2011)
'mZ',0.005,...											% Zbig quadratic mortality rate (mmolC^{-1} m^{3} d^{-1})
'eZ',0.1,...													% zoo excretion fraction (Zsmall and Zbig) (d^{-1}) (Baklouti et al., 2021)
'epsilon',0.75 ,...												% fraction of Zbig excretion that is available as regenerated PO4rients
'P_small_ini',Pini1,'P_big_ini',Pini2,'Z_small_ini',0.3,'Z_big_ini',0.3};	% initial biomass (mmolC m^{-3}) (In situ for P_small_ini with method from Marrec et al., 2018 and Tzortzis et al., 2021 & further)

[arg,flag]=ga_read_varargin(varargin,[{'nbdays_advec',20,'dt',0.2,'time',[],'upw_duration',1},default_parameters],{'plot'});
if length(Nsupply)>1 && isempty(arg.time), error('Give time if Nsupply is a vector'), end


%% -------------- Time

if isempty(arg.time), time=(0:arg.dt:arg.nbdays_advec)'; 
else, time=arg.time(:); arg.dt=time(2)-time(1); 
end
nb_time=length(time); 



%% -------------- Nsupply

if length(Nsupply)==1
	Nsupply_max=Nsupply; 
	Nsupply=zeros(nb_time,1);
	Nsupply(time<arg.upw_duration)=Nsupply_max; 
end



%% -------------- Initial conditions

PO4=time*NaN; 				% PO4rients (PO4 expressed in carbon) (mmolC m^{-3})
P_small=time*NaN; 			% small phyto biomass (mmolC m^{-3})
P_big=time*NaN; 			% large phyto biomass (mmolC m^{-3})
Z_small=time*NaN; 			% small zoo biomass (mmolC m^{-3})
Z_big=time*NaN;				% large zoo biomass (mmolC m^{-3})
u_big=time*NaN; 	        % P_big growth rate (Nnew-limited) (d^{-1})
u_small=time*NaN; 			% P_small growth rate (Nreg-limited) (d^{-1})
g_big=time*NaN; 			% Z_big grazing rate (P_big and Z_small limited) (d^{-1})
g_small=time*NaN;			% Z_small grazing rate (P_small limited) (d^{-1})


%%%% MODIFS
PP_big=time*NaN; 			% primary production for P_big (mmolC m^{-3} d^{-1})
PP_small=time*NaN;         % primary production for P_small (mmolC m^{-3} d^{-1})        
%%%% MODIFS

G_big1=time*NaN; 			% grazing from Z_big onto P_big (mmolC m^{-3} d^{-1})
G_big2=time*NaN; 			% grazing from Z_big onto Z_small (mmolC m^{-3} d^{-1})
G_small=time*NaN; 			% grazing from Z_small onto P_small (mmolC m^{-3} d^{-1})
excretion_Zbig=time*NaN; 	% Zbig excretion (mmolC m^{-3} d^{-1})
excretion_Zsmall=time*NaN; 	% Zsmall excretion (mmolC m^{-3} d^{-1})
death_Zbig=time*NaN; 		% Zbig mortality term (mmolC m^{-3} d^{-1})
death_Pbig=time*NaN; 		% Pbig mortality term (mmolC m^{-3} d^{-1})
regeneration_PO4=time*NaN;	% PO4rient regeneration (mmolC m^{-3} d^{-1})

P_small(1)=arg.P_small_ini; 
P_big(1)=arg.P_big_ini; 
Z_small(1)=arg.Z_small_ini; 
Z_big(1)=arg.Z_big_ini; 	
PO4(1)=C_nut;



%% -------------- Loop on time

for t=2:nb_time

	% growth and grazing rates

    %%%% MODIFS
    u_big=PO4(t-1)/(arg.kP_big+PO4(t-1))*arg.umax_big;
    u_small=PO4(t-1)/(arg.kP_small+PO4(t-1))*arg.umax_small;
    %%%% MODIFS
  
	g_big1=P_big(t-1)/(arg.kG_big+Z_small(t-1)+P_big(t-1))*arg.gmax_big; 	
	g_big2=Z_small(t-1)/(arg.kG_big+Z_small(t-1)+P_big(t-1))*arg.gmax_big;
	g_big(t-1)=g_big1+g_big2;
	g_small(t-1)=P_small(t-1)/(arg.kG_small+P_small(t-1))*arg.gmax_small;

	% fluxes
	
    %%%% MODIFS
    PP_big(t)=u_big*P_big(t-1);
    PP_small(t)=u_small*P_small(t-1); 
    %%%% MODIFS
    
	G_big1(t)=g_big1*Z_big(t-1); 
	G_big2(t)=g_big2*Z_big(t-1); 
	G_small(t)=g_small(t-1)*Z_small(t-1);
	death_Zbig(t)=arg.mZ*Z_big(t-1)^2;
	death_Pbig(t)=arg.mP*P_big(t-1);
	excretion_Zbig(t)=arg.eZ*Z_big(t-1);
	excretion_Zsmall(t)=arg.eZ*Z_small(t-1);
	regeneration_PO4(t)=arg.epsilon*excretion_Zbig(t)+excretion_Zsmall(t);
    
    
    %%%% MODIFS
    max_available_PO4 = PO4(t-1)+Nsupply(t)*arg.dt+regeneration_PO4(t)*arg.dt;
	if PP_small(t)>max_available_PO4/arg.dt, PP_small(t)=max_available_PO4/arg.dt; end
    max_available_PO4=max_available_PO4-PP_small(t)*arg.dt;
	if PP_big(t)>max_available_PO4/arg.dt, PP_big(t)=max_available_PO4/arg.dt; end
    %%%% MODIF
	
    
	% carbon-based nutrients and biomass
    
    %%%% MODIFS
	PO4(t)=PO4(t-1)+Nsupply(t)*arg.dt+regeneration_PO4(t)*arg.dt-PP_small(t)*arg.dt-PP_big(t)*arg.dt; PO4(PO4<=0)=0;
    %%%% MODIFS
    
    %%%% MODIFS
    P_big(t)=P_big(t-1)+PP_big(t)*arg.dt-G_big1(t)*arg.dt-death_Pbig(t)*arg.dt ;P_big(P_big<=0)=0;
	P_small(t)=P_small(t-1)+PP_small(t)*arg.dt-G_small(t)*arg.dt; P_small(P_small<=0)=0;
    %%%% MODIFS
    
	Z_small(t)=Z_small(t-1)+G_small(t)*arg.dt-G_big2(t)*arg.dt-excretion_Zsmall(t)*arg.dt; Z_small(Z_small<=0)=0;
	Z_big(t)=Z_big(t-1)+G_big1(t)*arg.dt+G_big2(t)*arg.dt-excretion_Zbig(t)*arg.dt-death_Zbig(t)*arg.dt; Z_big(Z_big<=0)=0;

end

% save('outputs/maxPP_big1', 'maxPP_big1');
% save('outputs/maxPP_big2', 'maxPP_big2');
% save('outputs/maxPP_small1', 'maxPP_small1');
% save('outputs/maxPP_small2', 'maxPP_small2');

%% -------------- Ouputs

units=struct('time','days','Nsupply','mmolC m^{-3} d^{-1}',...
	'P_small','mmolC m^{-3}','P_big','mmolC m^{-3}','Z_small','mmolC m^{-3}','Z_big','mmolC m^{-3}',...
	'PO4','mmolC m^{-3}','Chl','mg m^{-3}','PP','gC m^{-3}/yr',...
	'u_small','d^{-1}','u_big','d^{-1}','g_small','d^{-1}','g_big','d^{-1}', 'PP_small', 'gC m^{-3}/yr', 'PP_big','gC m^{-3}/yr');
output=struct('units',units,'time',time,'Nsupply',Nsupply,...
	'P_small',P_small,'P_big',P_big,'Z_small',Z_small,'Z_big',Z_big,'PO4',PO4,...
	'Chl',P_small*12/arg.cChl_small+P_big*12./arg.cChl_big,'PP',(PP_big+PP_small)*1E-3*12*365.25,...
	'u_small',u_small,'u_big',u_big,'g_small',g_small,'g_big',g_big, 'PP_small', PP_small,'PP_big', PP_big, 'attributs',struct('arg',arg));
varargout={output}; varargout=varargout(1:nargout);


% save('outputs/u_big','u_big')
save('outputs/P_big','P_big')
% save('outputs/PP_big','PP_big')
%% -------------- Figures

if flag.plot

    %Temporal evolution
	figure, hold on
	plot(output.time,output.P_small,'LineWidth',2, 'Color','g')
	plot(output.time,output.P_big,'LineWidth',2,'Color', '#77AC30')
	plot(output.time,output.Z_small,'LineWidth',2,'Color','#4DBEEE')
	plot(output.time,output.Z_big,'LineWidth',2, 'Color','b')
    plot(output.time,output.PO4,'LineWidth',2, 'Color','m')
	ylabel(output.units.P_small)
	if min(output.time)>datenum(1900,1,1), datetick('x','keeplimits'), xlabel('Time')
	else, xlabel('Time (days)')
	end
	legend({'P\_small','P\_big','Z\_small','Z\_big', 'PO4'})
	title('Model output (plankon concentration over time)')
    xlim([min(output.time) max(output.time)]);
    ylim([0 1]);

    %Monod equation
    figure, hold on
    N_theo_PO4 = linspace(0,10,length(arg.time));
    f_monod2(N_theo_PO4,arg.umax_small,arg.kP_small);
    f_monod2(N_theo_PO4,arg.umax_big,arg.kP_big);
    f_monod2(output.PO4,arg.umax_small,arg.kP_small);
    f_monod2(output.PO4,arg.umax_big,arg.kP_big);
    xlabel('[PO4] mmolC.m³')
    ylabel('\mu d^{-1}');
    legend({'P\_small\_theo','P\_big\_theo','P\_small\_mod','P\_big\_mod' })
     
end

return