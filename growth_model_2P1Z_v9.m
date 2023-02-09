function varargout=growth_model_2P1Z_v9(Nsupply,C_nut,varargin)
%% GA_MODEL_2P1Z_FROMPSUPPLY

% Growth-advection model by Messié & Chavez (2017), adapted to the Mediterranean. 
% Required input:
% 	Nsupply expressed in mmolC/m3/d, all units are carbon-based.
%
% Optional inputs:
% 'nbdays_advec'	number of days during which the model is run
% 'dt'				time step
% 'time'			can replace nbdays_advec and dt (required if Nsupply is a vector)
% 'plot'			displays the plankton model outputs as a function of time
%
% Monique Messié, 2021 for public version
% Reference: Messié, M., & Chavez, F. P. (2017). PO4rient supply, surface currents, and plankton dynamics predict zooplankton hotspots 
%					in coastal upwelling systems. Geophysical Research Letters, 44(17), 8979-8986, https://doi.org/10.1002/2017GL074322
% Differences with Messié and Chavez (2021):
%		Nutrient supply is Phosphate instead Nitrate
%		One nutrient compartiment
%       One zooplankton compartiment
%       Add conversion factor : ressources to consummator

%% -------------- Default parameters

% For now we define :
% 1 = Small phytoplankton (SYNECO, PICO) < 2um
% 2 = Big phytoplankton (MICRO) >20um

default_parameters={...
'umax1',1.9872,...  % maximum growth rates of P1 (d^{-1})
'umax2',2.7648,...  % maximum growth rates of P2 (d^{-1})
'gmax1',1.4226,...  % maximum grazing rates of Z on P1 (d^{-1})
'gmax2',1.4226,...  % maximum grazing rates of Z on P2 (d^{-1})
'kP1',1,...	        % half-saturation constant for P1 on PO4 (mmolC m^{-3})
'kP2',2,...	        % half-saturation constant for P2 on PO4 (mmolC m^{-3})
'kG',5,...			% half-saturation constant for Z on P (mmolC m^{-3})
'mP',0,...			% P2 mortality rate (default 0 ie no P2 sinking) (d^{-1})
'mZ',0.007,...		% Z2 quadratic mortality rate (mmolC^{-1} m^{3} d^{-1})
'eZ',0.1,...	    % zoo excretion rate (Z1 and Z2) (d^{-1})
'alpha1',0.7,...    % conversion factor from PO4 to P1
'alpha2',0.7,...    % conversion factor from PO4 to P2
'gamma1',0.7,...    % conversion factor from P1 to Z
'gamma2',0.5,...    % conversion factor from P2 to Z
'epsilon',0.75 ,... % fraction of Z excretion that is available as regenerated PO4
'P1_ini',0.6,...  % initial biomass of P1 (mmolC m^{-3})
'P2_ini',0.1,...  % initial biomass of P2 (mmolC m^{-3})
'Z_ini',0.6;	    % initial biomass of Z (mmolC m^{-3})
};	

% Call parameters as arg.parameters
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
PO4=time*NaN; 			% PO4 (PO4 expressed in carbon) (mmolC m^{-3})
P1=time*NaN; 			% P1 biomass (mmolC m^{-3})
P2=time*NaN; 			% P2 biomass (mmolC m^{-3})
Z=time*NaN; 			% Z biomass (mmolC m^{-3})
u1=time*NaN; 			% P1 growth rate (Nreg-limited) (d^{-1})
u2=time*NaN; 	        % P2 growth rate (Nnew-limited) (d^{-1})
g1=time*NaN; 			% Functional response on P1 to Z: Holling type II (P2 and P1 limited=) (d^{-1})
g2=time*NaN; 			% Functional response on P2 to Z: Holling type II (P2 and P1 limited=) (d^{-1})
PP1=time*NaN;           % primary production for P1 (mmolC m^{-3} d^{-1})  
PP2=time*NaN; 			% primary production for P2 (mmolC m^{-3} d^{-1})

%Flux
G1=time*NaN; 			% grazing from Z onto P1 (mmolC m^{-3} d^{-1})
G2=time*NaN; 			% grazing from Z onto P2 (mmolC m^{-3} d^{-1})
exc=time*NaN; 	        % Z excretion (mmolC m^{-3} d^{-1})
d_Z=time*NaN; 		    % Z mortality term (mmolC m^{-3} d^{-1})
reg=time*NaN;	        % PO4 regeneration (mmolC m^{-3} d^{-1})
p_fec=time*NaN;         % Fecal pelets (mmolC m^{-3} d^{-1})
Export=time*NaN;        % Export (biomass out of system) mmolC m^{-3} d^{-1})

%Initial values
P1(1)=arg.P1_ini; 
P2(1)=arg.P2_ini; 
Z(1)=arg.Z_ini;
PO4(1)=C_nut;

%% -------------- Loop on time

for t=2:nb_time

	% growth rates
    u1=PO4(t-1)/(arg.kP1+PO4(t-1))*arg.umax1;
    u2=PO4(t-1)/(arg.kP2+PO4(t-1))*arg.umax2;
  
    % functional response
    g1=P1(t-1)/(arg.kG+P1(t-1)+P2(t-1))*arg.gmax1; 	
	g2=P2(t-1)/(arg.kG+P1(t-1)+P2(t-1))*arg.gmax2;

	%%% FLUXES
    
    % primary production
    PP1(t)=arg.alpha1*u1*P1(t-1); 
    PP2(t)=arg.alpha2*u2*P2(t-1);
    
    % predation
	G1(t)=arg.gamma1*g1*Z(t-1); 
	G2(t)=arg.gamma2*g2*Z(t-1); 
    
    % mortality
	d_Z(t)=arg.mZ*Z(t-1)^2;
    
    % excretion
	exc(t)=arg.eZ*((1-arg.gamma1)*g1*Z(t-1)+(1-arg.gamma2)*g2*Z(t-1)); 
    
    % regeneration
	reg(t)=arg.epsilon*exc(t);
    
    % export
    p_fec(t) = (1-arg.epsilon)*exc(t);
    Export(t) =p_fec(t)+d_Z(t); 
    
    max_available_PO4 = PO4(t-1)+Nsupply(t)*arg.dt+reg(t)*arg.dt;
	if PP1(t)>max_available_PO4/arg.dt, PP1(t)=max_available_PO4/arg.dt; end
    max_available_PO4=max_available_PO4-PP1(t)*arg.dt;
	if PP2(t)>max_available_PO4/arg.dt, PP2(t)=max_available_PO4/arg.dt; end
    
	%%% BIOMASSES
    
    % phosphate
	PO4(t)=PO4(t-1)+Nsupply(t)*arg.dt+reg(t)*arg.dt-PP1(t)*arg.dt-PP2(t)*arg.dt; PO4(PO4<=0)=0;
    
    % phytoplanktons
	P1(t)=P1(t-1)+PP1(t)*arg.dt-G1(t)*arg.dt; P1(P1<=0)=0;
    P2(t)=P2(t-1)+PP2(t)*arg.dt-G2(t)*arg.dt ;P2(P2<=0)=0;
    
    % zooplankton
    Z(t)=Z(t-1)+G1(t)*arg.dt+G2(t)*arg.dt-exc(t)*arg.dt-d_Z(t)*arg.dt; Z(Z<=0)=0;
    
end

% save('outputs/maxPP21', 'maxPP21');

%% -------------- Ouputs

units=struct('time','days','Nsupply','mmolC m^{-3} d^{-1}',...
	'P1','mmolC m^{-3}','P2','mmolC m^{-3}','Z','mmolC m^{-3}',...
	'PO4','mmolC m^{-3}','PP','gC m^{-3}/yr',...
	'u1','d^{-1}','u2','d^{-1}','g1','d^{-1}','g2','d^{-1}', 'PP1', 'gC m^{-3}/yr', 'PP2','gC m^{-3}/yr', 'Excretion','mmolC m^{-3} d^{-1}', 'Regeneration','mmolC m^{-3} d^{-1}', 'Export', 'mmolC m^{-3} d^{-1}','d_Z', 'mmolC m^{-3} d^{-1}', 'p_fec','mmolC m^{-3} d^{-1}');

output=struct('units',units,'time',time,'Nsupply',Nsupply,...
	'P1',P1,'P2',P2,'Z', Z,'PO4',PO4,...
    'PP',(PP2+PP1)*1E-3*12*365.25,...
	'u1',u1,'u2',u2,'g1',g1,'g2',g2, 'PP1', PP1,'PP2', PP2, 'exc', exc, 'reg', reg, 'Export', Export,'d_Z',d_Z, 'p_fec', p_fec, 'attributs',struct('arg',arg));

varargout={output}; varargout=varargout(1:nargout);

%% -------------- Figures

if flag.plot

    %Temporal evolution
	figure, hold on
	plot(output.time,output.P1,'LineWidth',2, 'Color','g')
	plot(output.time,output.P2,'LineWidth',2,'Color', '#77AC30')
	plot(output.time,output.Z,'LineWidth',2,'Color','#00FFFF')
    plot(output.time,output.PO4,'LineWidth',2, 'Color','m')
	ylabel(output.units.P1)
	if min(output.time)>datenum(1900,1,1), datetick('x','keeplimits'), xlabel('Time')
	else, xlabel('Time (days)')
	end
	legend({'P\_small','P\_big','Z','PO4'})
	title('Model output (plankon concentration over time)')
    xlim([min(output.time) max(output.time)]);
    %ylim([0 1.2]);
    
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
    f_monod2(N_theo_PO4,arg.umax1,arg.kP1);
    f_monod2(N_theo_PO4,arg.umax2,arg.kP2);
    f_monod2(output.PO4,arg.umax1,arg.kP1);
    f_monod2(output.PO4,arg.umax2,arg.kP2);
    xlabel('[PO4] mmolC.m³')
    ylabel('\mu d^{-1}');
    legend({'P\_1\_theo','P\_2\_theo','P\_1\_mod','P\_2\_mod' })
    
%     %Portrait de phase
%     figure, hold on
%     % Calcul des dérivées de P1 et P2 par rapport au temps
%     dp1_dt = gradient(output.P1, arg.dt);
%     dp2_dt = gradient(output.P2, arg.dt);
%     
%     % Tracé du portrait de phase
%     plot(output.P2, output.P1)
%     xlabel('P\_big')
%     ylabel('P\_small')
%     
%     % Tracé du champ de vecteur sur le portrait de phase
%     quiver(output.P2, output.P1, dp2_dt, dp1_dt)
     
end

return