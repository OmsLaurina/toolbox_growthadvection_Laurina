function varargout=ga_model_pcolor_param_gmax_v8(Nsupply,C_nut,varargin)

%% -------------- Default parameters (see Messié & Chavez, 2017 suppl inf)

% 1 = Psmall
% 2 = Pbig

default_parameters={...
'umax_1',1.9872,'umax_2',2.7648,'gmax1',1.4226,'gmax2',1.4226,...% maximum growth and grazing rates (d^{-1}) (Baklouti et al., 2021)
'cChl_1',200,'cChl_2',50,...								% C:Chl ratios for P1 and P2 (only used to calculate Chl) (Lazzari et al., 2012)
'kP_1',1,...											% half-saturation constant for P1 on PO4 (mmolC m^{-3}) (Pulido-Villena et al., 2021) )
'kP_2',2,...												% half-saturation constant for P2 on PO4 (mmolC m^{-3}) (Pulido-Villena et al., 2021)
'kG',5,...										% half-saturation constant for Z on P (mmolC m^{-3}) (Auger et al., 2011)
'mP',0,...														% P2 mortality rate (default 0 ie no P2 sinking) (d^{-1}) (Auger et al., 2011)
'mZ',0.005,...											% Z2 quadratic mortality rate (mmolC^{-1} m^{3} d^{-1})
'eZ',0.1,...													% zoo excretion fraction (Z1 and Z2) (d^{-1}) (Baklouti et al., 2021)
'epsilon',0.75 ,...												% fraction of Z excretion that is available as regenerated PO4rients
'P_1_ini',0.6,'P_2_ini',0.1,'Z_ini',0.6};	% initial biomass (mmolC m^{-3}) (In situ for P_1_ini with method from Marrec et al., 2018 and Tzortzis et al., 2021 & further)

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
u_1=time*NaN; 			% P_1 growth rate (d^{-1})
u_2=time*NaN; 	        % P_2 growth rate (d^{-1})
PP_1=time*NaN;         % primary production for P_1 (mmolC m^{-3} d^{-1}) 
PP_2=time*NaN; 			% primary production for P_2 (mmolC m^{-3} d^{-1})

%Flux
G1=time*NaN; 			% grazing from Z onto P2 (mmolC m^{-3} d^{-1})
G2=time*NaN; 			% grazing from Z onto P1 (mmolC m^{-3} d^{-1})
excretion_Z=time*NaN; 	% Z excretion (mmolC m^{-3} d^{-1})
death_Z=time*NaN; 		% Z2 mortality term (mmolC m^{-3} d^{-1})
regeneration_PO4=time*NaN;	% PO4 regeneration (mmolC m^{-3} d^{-1})
p_fec=time*NaN;
Export=time*NaN;

%Valeurs initiales
P_1(1)=arg.P_1_ini; 
P_2(1)=arg.P_2_ini; 
Z(1)=arg.Z_ini;
PO4(1)=C_nut;

%% -------------- Loop on time

n = 10;

min_gmax2 = 1;
max_gmax2 = 5;
l_gmax2 = linspace(min_gmax2,max_gmax2,n);

min_gmax1 = 1;
max_gmax1 = 5;
l_gmax1 = linspace(min_gmax1,max_gmax1,n);

ratio=zeros(length(l_gmax2),length(l_gmax1));
ratio2=zeros(length(l_gmax2),length(l_gmax1));

i = 1;
for gmax2=l_gmax2
    
    j=1;
    for gmax1=l_gmax1
        
        for t=2:nb_time

        % growth and grazing rates
        u_1=PO4(t-1)/(arg.kP_1+PO4(t-1))*arg.umax_1;
        u_2=PO4(t-1)/(arg.kP_2+PO4(t-1))*arg.umax_2;

        g1=P_1(t-1)/(arg.kG+P_1(t-1)+P_2(t-1))*gmax1;
        g2=P_2(t-1)/(arg.kG+P_1(t-1)+P_2(t-1))*gmax2;

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
        
        ratio(i,j) = P_1(length(P_1)-1)/(P_1(length(P_1)-1)+P_2(length(P_2)-1));
        ratio2(i,j) = P_2(length(P_2)-1)/(P_1(length(P_1)-1)+P_2(length(P_2)-1));
    j = j+1;
    end
    
    i = i+1;
    i
end

save('outputs/P_1','P_1');

%% -------------- Ouputs

output=struct('ratio', ratio, 'ratio2', ratio2,'attributs',struct('arg',arg));
varargout={output}; varargout=varargout(1:nargout);

% save('outputs/output','output')
%% -------------- Figures

if flag.plot
    
    figure,hold on
    
    x = linspace(min_gmax2,max_gmax2,n);
    y = linspace(min_gmax1,max_gmax1,n);
    [X,Y] = meshgrid(x,y);
    pcolor(X,Y,output.ratio)
    xlabel('gmax\_big')
    ylabel('gmax\_small')
    shading interp
    colorbar
    caxis([0 1])
    title('P\_small/(P\_small+P\_big)')
    
    figure,hold on
    
    x = linspace(min_gmax2,max_gmax2,n);
    y = linspace(min_gmax1,max_gmax1,n);
    [X,Y] = meshgrid(x,y);
    pcolor(X,Y,output.ratio2)
    xlabel('gmax\_big')
    ylabel('gmax\_small')
    shading interp
    colorbar
    caxis([0 1])
    title('P\_big/(P\_big+P\_small)')
    
end

return