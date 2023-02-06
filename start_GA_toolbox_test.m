%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 growth-advection model - TEST              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% from https://github.com/messiem/toolbox_GrowthAdvection    %%%
%%% See start_GA_toolbox.m                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Use function : ga_*.m                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Author : Laurina Oms                                       %%%
%%% Creation : 18/10/2022                                      %%%     
%%% Region : California coastal current                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

clear all, close all

%%Directories and paths%%
%addpath('Documents/toolbox_GrowthAdvection_Laurina/')

disp('start_GA_toolbox_test : Start simulation');

%directory of the Ariane workplace
global dir_ariane_global
dir_ariane_global='Ariane_workplace/';

%directory where output will be save
global dir_output_global
dir_output_global='outputs/';

%choose de variable you want to observe (Z_big, Z_small, P_big, P_small,
%Chl...)
var_model = 'P_small';

%loop to define caxis limits in function of the variable choose
if strcmp(var_model, 'Z_big') 
    caxis_lim = [0 20];
elseif strcmp(var_model, 'Z_small') 
    caxis_lim = [1 2];
elseif strcmp(var_model, 'P_big') 
    caxis_lim = [0 30];
elseif strcmp(var_model, 'P_small') 
    caxis_lim = [0 1];
elseif strcmp(var_model, 'Chl') 
    caxis_lim = [0 20];
end
    

%% -----Run the plankton model

n_days = 60; %nb de jours de simulation
dt = 0.2; %pas de temps
time=(0:dt:n_days); %temps
Psupply = ones(length(time),1)*0.01;
output=ga_model_2P2Z_v7(Psupply,0.5,0.6,0.1,'kP_small',1,'kP_big',2,'gmax_small',1.4226,'gmax_big',1.1226,'time',time,'plot');
%print('-djpeg','-r300','Documents/toolbox_GrowthAdvection_Monique/outputs/plankton_model.jpg')

%% -----Compute position over time using Ariane
lon_min = 0;
lat_min = 36;
lon_max = 7;
lat_max = 40;
dt = 0.2;
nbdays_advec = 13;
time0 = datenum(2018,4,30);
name_curr='PROTEVS_bio';	% using currents into Ariane_workplace/currents_data

%Coastlines
load('inputs/gumby.mat')
coast_x = ncst(:,1);
coast_y = ncst(:,2);

%Initialization of particule's position
xvalue = 2:0.1:5; %lat
yvalue = 37:0.1:40;     %lon
[X,Y] = meshgrid(xvalue, yvalue);
point_1 = [X(:) Y(:)];

%use ga_advection_ariane
positions=ga_advection_ariane(point_1,name_curr,'dt',dt,'time0',time0,'nbdays_advec',nbdays_advec);
time2D=repmat(positions.time',size(point_1,1),1);

%For just 3 particules
% mat_positions_ini=[[1 37];[1 37.5];[2 39]];
% positions=ga_advection_ariane(mat_positions_ini,'PROTEVS_bio2','dt',0.2,'time0',datenum(2018,4,30),'nbdays_advec',13);
% time2D=repmat(positions.time',size(mat_positions_ini,1),1);

figure, hold on
        scatter(positions.lon2D(:),positions.lat2D(:),10,time2D(:),'filled')
        plot(coast_x,coast_y,'k')
        xlim([lon_min lon_max]), ylim([lat_min lat_max])
        hbar=colorbar; datetick(hbar,'y','keeplimits')
        set(get(hbar,'title'),'string','time');
        xlabel('Longitude'), ylabel('Latitude')
        title('Current trajectories initialized on April 4, 2018')

%% -----Compute one daily run of the GA model

% % Load inputs and set run date
% load('inputs/Nsupply_PROTEVS.mat')% load the Nsupply forcing
% 
% load('inputs/gumby.mat')
% coast_x = ncst(:,1);
% coast_y = ncst(:,2);	
% 
% options_plankton_model={'kP_small',1,'kP_big',2,'gmax_small',1.4226,'gmax_big',1.1226};										% using currents toolbox_* into Ariane_workplace/currents_data	
% 
% % Construct init structure
% init=struct();
% %init.Nsupply=nan(length(Nsupply_PROTEVS.lat),1);
% [LON,LAT] = meshgrid(Nsupply_PROTEVS.lon, Nsupply_PROTEVS.lat);
% Nsupply_q = interp2(LON,LAT,Nsupply_PROTEVS.Nsupply_ini',xvalue,yvalue);
% 
% init.lon=xvalue';
% init.lat=yvalue';
% init.Nsupply=Nsupply_q';
% 
% %for ilat=1:length(init.lat), init.Nsupply(ilat)=interp1(Nsupply.time,Nsupply.Nsupply(ilat,:),time0); end
% 
% % for ilat=1:length(Nsupply.lat), init.Nsupply(ilat)=interp1(Nsupply.time,Nsupply.Nsupply(ilat,:),time0); end
% 
% % Run the growth-advection program
% zoo=ga_growthadvection(init,name_curr,time0,'options_plankton_model',options_plankton_model);
% 
% %to get the variable field wanted
% values = getfield(zoo,var_model);
% units = getfield(zoo.units,var_model);
% 
% % Figure
% % Note: pixels over land are due to current interpolation to the coast (so pixels overlay land). 
% % The functions that concatenes daily runs to generate maps moves them over to the coastline again (see ga_concatene_runs)
% figure, hold on
% 	scatter(flipud(zoo.lon2D(:)),flipud(zoo.lat2D(:)),5,flipud(values(:)),'filled')
% 	plot(coast_x,coast_y,'k')
% 	xlim([lon_min lon_max]), ylim([lat_min lat_max])
% 	hbar=colorbar; caxis(caxis_lim)
% 	set(get(hbar,'title'),'string',units);
% 	xlabel('Longitude'), ylabel('Latitude')
% 	title([var_model, ' from GA run initialized'], 'Interpreter', 'none')

%% -----Compute full GA model

% Load inputs and set run options
 load('inputs/Nsupply_PROTEVS.mat')                               % load the Nsupply forcing for the 2008 season (Nsupply.name = CCMP3km)

options_plankton_model={'kP_small',1,'kP_big',2,'gmax_small',1.4226,'gmax_big',1.1226}; % krill parameterization

name_curr='PROTEVS_bio';                % using currents toolbox_* into Ariane_workplace/currents_data

% Run the full GA model
ga_full_GArun(Nsupply_PROTEVS,name_curr,'options_plankton_model',options_plankton_model) 

% Look at outputs: example May 2008 (by default days are all at 15) - Reproduces Fig. 1c
% Note - positions are shifted for pcolor so that pixels are centered on the position (pcolor considers positions to be the bottom left corner)
load('outputs/zoo_Lagrangian.mat')
load('inputs/gumby.mat')
coast_x = ncst(:,1);
coast_y = ncst(:,2);	
itime=zoo_all.time==datenum(2018,4,30);

%to get the variable field wanted
values2D = getfield(zoo_all,var_model);
units2D = getfield(zoo_all.units,var_model);

figure, hold on
pcolor(zoo_all.lon2D-0.125/2,zoo_all.lat2D-0.125/2,values2D(:,:,itime)), shading flat
plot(coast_x,coast_y,'k')
xlim([lon_min lon_max]), ylim([lat_min lat_max])
hbar=colorbar; caxis(caxis_lim)
set(get(hbar,'title'),'string',units2D);
xlabel('Longitude'), ylabel('Latitude'), title([var_model, ' mapped output'], 'Interpreter', 'none')
%print('-djpeg','-r300','Documents/toolbox_GrowthAdvection/outputs/GAmap_monthly_200805.jpg')

disp('start_GA_toolbox_test : done');
toc

