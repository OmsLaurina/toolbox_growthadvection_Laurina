function varargout=pcolor_umax1_Psupp_v9(Nsupply,C_nut,varargin)

%% -------------- Default parameters

% For the time we define :
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
'gamma1',0.7,...    % conversion factor from P1 to Z
'gamma2',0.7,...    % conversion factor from P2 to Z
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

n = 50;

min_umax1 = 0;
max_umax1 = 5;
l_umax1 = linspace(min_umax1,max_umax1,n);

min_Psupp = 0;
max_Psupp = 1;
l_Psupp = linspace(min_Psupp,max_Psupp,n);

P_1=zeros(length(l_umax1),length(l_Psupp));
P_2=zeros(length(l_umax1),length(l_Psupp));

i = 1;
for umax1=l_umax1
    
    j=1;
    for Psupp=l_Psupp
        
        for t=2:nb_time
            
            % growth rates
            u1=PO4(t-1)/(arg.kP1+PO4(t-1))*umax1;
            u2=PO4(t-1)/(arg.kP2+PO4(t-1))*arg.umax2;


            % functional response
            g1=P1(t-1)/(arg.kG+P1(t-1)+P2(t-1))*arg.gmax1; 	
            g2=P2(t-1)/(arg.kG+P1(t-1)+P2(t-1))*arg.gmax2;


            %%% FLUXES

            % primary production
            PP1(t)=u1*P1(t-1); 
            PP2(t)=u2*P2(t-1);

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

            max_available_PO4 = PO4(t-1)+Nsupply(t)*Psupp*arg.dt+reg(t)*arg.dt;
            if PP1(t)>max_available_PO4/arg.dt, PP1(t)=max_available_PO4/arg.dt; end
            max_available_PO4=max_available_PO4-PP1(t)*arg.dt;
            if PP2(t)>max_available_PO4/arg.dt, PP2(t)=max_available_PO4/arg.dt; end

            %%% BIOMASSES

            % phosphate
            % Define Nsupply = 1
            PO4(t)=PO4(t-1)+Nsupply(t)*Psupp*arg.dt+reg(t)*arg.dt-PP1(t)*arg.dt-PP2(t)*arg.dt; PO4(PO4<=0)=0;

            % phytoplanktons
            P1(t)=P1(t-1)+PP1(t)*arg.dt-G1(t)*arg.dt; P1(P1<=0)=0;
            P2(t)=P2(t-1)+PP2(t)*arg.dt-G2(t)*arg.dt ;P2(P2<=0)=0;

            % zooplankton
            Z(t)=Z(t-1)+G1(t)*arg.dt+G2(t)*arg.dt-exc(t)*arg.dt-d_Z(t)*arg.dt; Z(Z<=0)=0;
        
        end
        
        P_1(i,j) = P1(length(P1)-1);
        P_2(i,j) = P2(length(P2)-1);
        j=j+1;
    end
    i = i+1;
    i
end

save('outputs/P_2', 'P_2');
save('outputs/P_1', 'P_1');

%% -------------- Ouputs


output=struct('P_1', P_1, 'P_2', P_2, 'time',time,'PO4', PO4, 'attributs',struct('arg',arg));
varargout={output}; varargout=varargout(1:nargout);


%% -------------- Figures

if flag.plot
    
    figure,hold on
    box on
    x = linspace(min_umax1,max_umax1,n);
    y = linspace(min_Psupp,max_Psupp,n);
    [X,Y] = meshgrid(x,y);
    pcolor(X,Y,output.P_1)
    xlabel('umax1')
    ylabel('Psupp')
    shading interp
    colorbar
    %caxis([0 0.01])
    title('P1')
    
    figure,hold on
    box on
    x = linspace(min_umax1,max_umax1,n);
    y = linspace(min_Psupp,max_Psupp,n);
    [X,Y] = meshgrid(x,y);
    pcolor(X,Y,output.P_2)
    xlabel('umax1')
    ylabel('Psupp')
    shading interp
    colorbar
    %caxis([0 0.01])
    title('P2')
    
end

return