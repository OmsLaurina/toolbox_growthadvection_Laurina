%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                     growth-advection model - toolbox Laurina                                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% https://bitbucket.org/kientzn/toolbox_stage/src/master/                                                            %%%
%%% Use function : ga_*.m                                                                                              %%%
%%% Author : adapted by Laurina Oms                                                                                    %%%
%%% Creation : 14/11/2022                                                                                              %%%     
%%% Region : South Balearic Island                                                                                     %%%
                            
%%% First part : Phytoplancton biomasse (in situ and extrapolated) advected along lagrangian trajectories from Ariane  %%% 
%%% Second part : Phytoplankton biomasse (modelized from growth model & insitu) advected along lagrangian trajectories %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all

tic

disp('start_GA_toolbox_test : ...');

%run scripts to make variables needed
run('transfo_currents.m');
% run('chl_map.m');
% run('sst_map.m');
% run('abundance_map.m');
% run('Abundance_relations.m');
% run('relation_SST_T_tsg.m');
% run('test_ini_map.m');

%% ------------------ Set up

%%% Set directory 

% where Ariane is installed (used in most functions)
% There must be a directory "currents_data" inside dir_ariane_global where currents netcdf files are saved (see ga_write_ariane_currents)
global dir_ariane_global
dir_ariane_global='Ariane_workplace/';

% where outputs will be saved (used in ga_full_GArun)
global dir_output_global
dir_output_global='outputs/';

%%% Define variables

%Simulation time 
time_start = '2018,1,1';
time_end = '2018,12,31';

%coordonates of the simulate region
lon_min = 0;
lat_min = 36;
lon_max = 7;
lat_max = 40;

%Choose the size of your mesh grid (in °)
m_g = 0.1;

%Index of cytometry data during hippodrome
index_hipp = 407:511;

%% GENERATE 1 PARTICULE FOR EACH GRID MESH BEFORE ADVECTION

% %mesh grid 
% xvalue = lon_min:m_g:lon_max; %1 particule tout les m_g°
% yvalue = lat_min:m_g:lat_max;
% [X,Y] = meshgrid(xvalue, yvalue);
% point_1 = [X(:) Y(:)]; %initial positions
% 
% SST_MUR=struct();
% SST_MUR.sst = sst_MUR;
% SST_MUR.lon_sst = lon_sst_MUR;
% SST_MUR.lat_sst = lat_sst_MUR;
% SST_MUR.time_sst = datenum(time_start):datenum(time_end);
% save('inputs/SST_MUR','SST_MUR')
% save('inputs/time','time')
% 
% %mask
% iok = ~isnan(mat_micro);
% mat_positions_ini=[X(iok) Y(iok)];
% %save('inputs/mat_positions_ini','mat_positions_ini')
% positions=ga_advection_ariane(mat_positions_ini,'PROTEVS_','dt',0.2,'time0',time(120),'nbdays_advec',1);
% %save('inputs/positions','positions');
% time2D=repmat(positions.time',size(point_1,1),1);
% save('inputs/iok', 'iok')
% save('inputs/time2D', 'time2D')

% %% ------------------ FIRST PART 
% %%% COMPUTE ABUNDANCE OVER TIME USING ARIANE (ONLY ADVECTION = no growth model )
% 
% %Define phytoplankton settings for this scripts
% name_group = 'pico';
% parameters = phytoP_settings(name_group);
% 
% %position of phytoP group in dataCYTONEW
% pos_group = parameters.pos_group;
% 
% %define caxis limits in function of the variable choose
% caxis_lim = parameters.caxis_lim;
% 
% %Define the mean size magnitude of phytoplankton cells (um)
% size_magnitude= parameters.size_magnitude;
% %---------------------- SST ----------------------%
% 
% tok = ~isnan(mat_sst);
% %mat_sst = mat_sst(tok);
% sst_matrix = reshape(mat_sst(:,:),[size(mat_sst,1)*size(mat_sst,2),1]);
% for it = 1:length(time2D(1,:))-1
%     sst_matrix= [sst_matrix reshape(mat_sst(:,:),[size(mat_sst,1)*size(mat_sst,2),1])];
% end
%     
% %---------------------- Particles advected ----------------------%
% 
% pok = ~isnan(mat_sst);
% % pok = pok.*(positions.lon2D(:,1));
% particles_matrix = reshape(pok(:,:),[size(pok,1)*size(pok,2),1]);
% for it = 1:length(time2D(1,:))-1
%     particles_matrix = [particles_matrix reshape(pok(:,:),[size(pok,1)*size(pok,2),1])];
% end  
% 
% %Averaged abundance
% %Compute the mean on rectangular patch from scattered data
% Datasize = length(lat_MUR)*length(lon_MUR);
% Lat = positions.lat2D(:,end-7);
% Lon = positions.lon2D(:,end-7);
% 
% %extrapolated biomasses
% Data = parameters.biomasse_extp.biomasse;
% n_array = 6;
% lat_edge = lat_min+(m_g*n_array)/2:m_g*n_array:lat_max+(m_g*n_array)/2;
% lon_edge = lon_min+(m_g*n_array)/2:m_g*n_array:lon_max+(m_g*n_array)/2;
% meanData = histcn([Lat(:) Lon(:)], lat_edge, lon_edge, 'AccumData', Data, 'Fun', @mean);
% meanData(meanData==0)=NaN;
% %in situ biomasses
% ab_cyto = parameters.biomasse_insitu.biomasse;
% 
% adv_biom=struct();
% adv_biom.lon=lon_edge;
% adv_biom.lat=lat_edge;
% adv_biom.time=positions.time;
% adv_biom.biom = Data;
% adv_biom.lon2D = positions.lon2D;
% adv_biom.lat2D = positions.lat2D;
% save('outputs/adv_biom','adv_biom')
% 
% %%% Figures
% figure('DefaultAxesFontSize',22)
% m_proj('mercator','lon',[lon_min  lon_max],'lat',[lat_min lat_max]);
% m_pcolor(lon_edge(1:end-1), lat_edge(1:end-1), meanData)
% shading flat
% hold on
% m_scatter(lon_abundance(index_hipp),lat_abundance(index_hipp),30,ab_cyto,'filled');
% m_quiver(X_curr,Y_curr,mean(u0(:,:,1:13),3,'omitnan'),mean(v0(:,:,1:13),3,'omitnan'),3,'r');
% m_usercoast('gumby','patch','w');
% m_grid('box','fancy','linestyle','-','gridcolor','none','backcolor','none');
% hold on
% cbar = colorbar;
% caxis(caxis_lim);
% colorTitleHandle = get(cbar,'Title');
% titleString = ({'Biomasse (mmol C.m⁻³)'});
% set(colorTitleHandle ,'String',titleString);
% title(name_group);

% %% ------------------ SECOND PART
%%% RUN THE FULL GA MODEL
dataCYTONEW = dlmread('inputs/data_CYTO_NEW.txt');
%%% Set run options
load('inputs/Nsupply_PROTEVS.mat')
options_plankton_model={'kP_small',1,'kP_big',2,'gmax_small',1.4226,'gmax_big',1.1226};
name_curr='PROTEVS_bio';									% using currents toolbox_* into Ariane_workplace/currents_data

%%% Run the full GA model
ga_full_GArun(Nsupply_PROTEVS,name_curr,'options_plankton_model',options_plankton_model)

load('outputs/zoo_Lagrangian')
load('outputs/zoo')
zoo_all.P_small = zoo_all.P_1;
zoo_all.P_big = zoo_all.P_2;

%%% Define mesh grid
load('inputs/curr_struc')
reso = curr.lon(2)-curr.lon(1);
xvalue_curr = min(curr.lon):reso:max(curr.lon);
yvalue_curr = min(curr.lat):reso:max(curr.lat);
[X_curr,Y_curr] = meshgrid(xvalue_curr, yvalue_curr);
Lat = zoo_all.lat2D(:,end-7); %equivalent to the middle of hippodrom period
Lon = zoo_all.lon2D(:,end-7);

%%%---------------------- Pico ----------------------%%%

%convert modeled biomasse mmol C/m3 into abundance (cell/cm3)
QC_pico = 0.26*exp(-5.8702)*(0.9*10000)^(0.9228)*0.86*1000;
QC_pico = QC_pico*1E-12/12.106;

zoo_ab = struct();
zoo_ab.P_small = (zoo_all.P_small./QC_pico)./1000000; %convert biomasse mmol C/m3 into abundance (cell/cm3)

%convert in situ abundances cell/cm³ into biomasse (mmol C/m3)
pico_ab = dataCYTONEW(index_hipp,10);
pico_ab = pico_ab.*1000000; %cell/m3
pico_biom = pico_ab.*QC_pico;

Data = zoo_all.P_small(:,end-15:end);

n_array = 3;
lat_edge = min(Lat)+(0.02*n_array)/2:0.02*n_array:max(Lat)+(0.02*n_array)/2;
lon_edge = min(Lon)+(0.02*n_array)/2:0.02*n_array:max(Lon)+(0.02*n_array)/2;
%meanData_marged = NaN(length(lat_edge),length(lon_edge));
meanData = histcn([Lat(:) Lon(:)], lat_edge, lon_edge, 'AccumData', Data, 'Fun', @mean);
meanData(meanData==0)=NaN;
% for i = 2:size(meanData,1)+1
%     for j = 1:size(meanData,2)
%         meanData_marged(i,j) = meanData(i-1,j);
%     end
% end      
%meanData_marged(meanData_marged==0)=NaN;

pico = struct();
pico.biom_mod = meanData;
pico.lon_mod = lon_edge;
pico.lat_mod = lat_edge;
pico.biom_situ = pico_biom;

for i = 1:length(dataCYTONEW)
lat_abundance(i) = dataCYTONEW(i,end);
lon_abundance(i) = dataCYTONEW(i,end-1);
end
pico.lon_situ = lon_abundance(index_hipp);
pico.lat_situ = lat_abundance(index_hipp);

pico.lon_interp = min(pico.lon_mod):0.02:max(pico.lon_mod);
pico.lat_interp = min(pico.lat_mod):0.02:max(pico.lat_mod);
[Xq, Yq] = meshgrid(pico.lon_interp,pico.lat_interp);
pico.biom_interp = interp2(pico.lon_mod(1:end-1),pico.lat_mod(1:end-1),pico.biom_mod,Xq, Yq);

save('outputs/pico_struc','pico');

%%%---------------------- FIGURES ----------------------%%%
figure('DefaultAxesFontSize',22)
m_proj('mercator','lon',[lon_min lon_max],'lat',[lat_min lat_max]);
%m_pcolor(lon_edge(1:size(meanData,2)),lat_edge(1:size(meanData,1)),meanData)
m_pcolor(Xq, Yq, pico.biom_interp)
shading flat
hold on
m_scatter(lon_abundance(index_hipp),lat_abundance(index_hipp),30,pico_biom,'filled');
m_quiver(X_curr,Y_curr,mean(u0(:,:,1:13),3),mean(v0(:,:,1:13),3),3,'r');
m_usercoast('gumby','patch','w');
m_grid('box','fancy','linestyle','-','gridcolor','none','backcolor','none');
hold on
cbar = colorbar;
caxis([0.7 0.8]);
colorTitleHandle = get(cbar,'Title');
titleString = ({'Biomasse (mmol C.m⁻³)'});
set(colorTitleHandle ,'String',titleString);

%%%---------------------- Micro ----------------------%%%

%convert modeled biomasse mmol C/m3 into abundance (cell/cm3)
QC_micro = 0.26*exp(-5.8702)*(90*10000)^(0.9228)*0.86*1000;
QC_micro = QC_micro*1E-12/12.106;

zoo_ab = struct();
zoo_ab.P_big = (zoo_all.P_big./QC_micro)./1000000; %convert biomasse mmol C/m3 into abundance (cell/cm3)

%convert in situ abundances cell/cm³ into biomasse (mmol C/m3)
micro_ab = dataCYTONEW(index_hipp,14);
micro_ab = micro_ab.*1000000; %cell/m3
micro_biom = micro_ab.*QC_micro;

Data = zoo_all.P_big(:,end);

n_array = 3;
lat_edge = min(Lat)+(0.02*n_array)/2:0.02*n_array:max(Lat)+(0.02*n_array)/2;
lon_edge = min(Lon)+(0.02*n_array)/2:0.02*n_array:max(Lon)+(0.02*n_array)/2;
%meanData_marged = NaN(length(lat_edge),length(lon_edge));
meanData = histcn([Lat(:) Lon(:)], lat_edge, lon_edge, 'AccumData', Data, 'Fun', @mean);
meanData(meanData==0)=NaN;
% for i = 2:size(meanData,1)+1
%     for j = 1:size(meanData,2)
%         meanData_marged(i,j) = meanData(i-1,j);
%     end
% end      
%meanData_marged(meanData_marged==0)=NaN;

micro = struct();
micro.biom_mod = meanData;
micro.lon_mod = lon_edge;
micro.lat_mod = lat_edge;
micro.biom_situ = micro_biom;
micro.lon_situ = lon_abundance(index_hipp);
micro.lat_situ = lat_abundance(index_hipp);

micro.lon_interp = min(micro.lon_mod):0.01:max(micro.lon_mod);
micro.lat_interp = min(micro.lat_mod):0.01:max(micro.lat_mod);
[Xq, Yq] = meshgrid(micro.lon_interp,micro.lat_interp);
micro_interp = interp2(micro.lon_mod(1:end-1),micro.lat_mod(1:end-1),micro.biom_mod,Xq, Yq);

micro.biom_interp = micro_interp;
save('outputs/micro_struc','micro');


%%%---------------------- FIGURES ----------------------%%%
figure('DefaultAxesFontSize',22)
m_proj('mercator','lon',[lon_min lon_max],'lat',[lat_min lat_max]);
%m_pcolor(lon_edge(1:size(meanData,2)),lat_edge(1:size(meanData,1)),meanData)
m_pcolor(Xq, Yq, micro.biom_interp)
shading flat
hold on
m_scatter(lon_abundance(index_hipp),lat_abundance(index_hipp),30,micro_biom,'filled');
m_quiver(X_curr,Y_curr,mean(u0(:,:,1:13),3,'omitnan'),mean(v0(:,:,1:13),3,'omitnan'),3,'r');
m_usercoast('gumby','patch','w');
m_grid('box','fancy','linestyle','-','gridcolor','none','backcolor','none');
hold on
cbar = colorbar;
caxis([0.3 0.4]);
colorTitleHandle = get(cbar,'Title');
titleString = ({'Biomasse (mmol C.m⁻³)'});
set(colorTitleHandle ,'String',titleString);

disp('start_GA_toolbox : done');

toc