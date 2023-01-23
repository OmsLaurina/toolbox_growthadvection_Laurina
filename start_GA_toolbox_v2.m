%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                     growth-advection model - toolbox Laurina                                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% https://bitbucket.org/kientzn/toolbox_stage/src/master/                                                            %%%
%%% Use function : ga_*.m                                                                                              %%%
%%% Author : adapted by Laurina Oms                                                                                    %%%
%%% Creation : 23/11/2022                                                                                              %%%     
%%% Region : South Balearic Island                                                                                     %%%
%%% Make the map of micro and pico biomasses from the full GA model                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all

tic

disp('start_GA_toolbox_test : ...');

% Run scripts to make all variables needed

run('transfo_currents.m');
run('chl_map.m');
%run('sst_map.m');
run('abundance_map.m');
run('Abundance_relations.m');
%run('relation_SST_T_tsg.m');
run('test_ini_map.m');

%% ------------------ Set up ------------------ %%

%%% Set directory 

% where Ariane is installed (used in most functions)

% There must be a directory "currents_data" inside dir_ariane_global where currents netcdf files are saved (see ga_write_ariane_currents)
global dir_ariane_global
dir_ariane_global='Ariane_workplace/';

% where outputs will be saved (used in ga_full_GArun)

global dir_output_global
dir_output_global='outputs/';

%%% Set the simulation's context

time_0 = datenum(2018,4,30);        % Time start
m_g = 0.02;                         % space step
dt = 0.2;                           % time step
nbdays_advec = 13;                  % number of days advected
nbdays_wanted = 3;                  % number of days I want to visualize (3 correspond to the days of the hippodrome)
nb_dt = nbdays_wanted/dt;           % number of time step for nbdays_wanted
dt_1 = ((nbdays_advec/dt)+1)-nb_dt; % start/end position in the matrice of the simulation (here correspond to the last 3 days of the cruise (=hippodrome))
dt_end = (nbdays_advec/dt)+1;       % but for ex if I want the first x days: dt_1 = 1; dt_end = nb_dt (where nbdays_wanted = x)

%coordonates of the simulate region

lon_min = 0;
lat_min = 36;
lon_max = 7;
lat_max = 40;

%Index of cytometry data during hippodrome

index_hipp = 407:511;

%% ------------------ RUN THE FULL GA MODEL ------------------ %%

%%% Set run options

load('inputs/Nsupply_PROTEVS.mat')
options_plankton_model={'dt', dt, 'nbdays_advec', nbdays_advec};
name_curr='PROTEVS_bio2';          % using currents toolbox_* into Ariane_workplace/currents_data

%%% Run the full GA model

ga_full_GArun(Nsupply_PROTEVS,name_curr,time_0,m_g,dt,nbdays_advec,'options_plankton_model',options_plankton_model);

load('outputs/zoo_Lagrangian')
load('outputs/zoo')

%%% Define mesh grid 

load('inputs/curr_struc')
reso = curr.lon(2)-curr.lon(1);
xvalue_curr = min(curr.lon):reso:max(curr.lon);
yvalue_curr = min(curr.lat):reso:max(curr.lat);
[X_curr,Y_curr] = meshgrid(xvalue_curr, yvalue_curr);
Lat = zoo_all.lat2D(:,end-nb_dt);
Lon = zoo_all.lon2D(:,end-nb_dt);

%%%%%%%%%%%%%% Pico %%%%%%%%%%%%%% 

%convert modeled biomasse mmol C/m3 into abundance (cell/cm3)

QC_pico = 0.26*exp(-5.8702)*(0.9*10000)^(0.9228)*0.86*1000;
QC_pico = QC_pico*1E-12/12.106;

zoo_ab = struct();
zoo_ab.P_small = (zoo_all.P_small./QC_pico)./1000000; %convert biomasse mmol C/m3 into abundance (cell/cm3)

%convert in situ abundances cell/cm³ into biomasse (mmol C/m3)

pico_ab = dataCYTONEW(index_hipp,10);
pico_ab = pico_ab.*1000000; %cell/m3
pico_biom = pico_ab.*QC_pico;

Data = zoo_all.P_small(:,dt_1:dt_end);

n_array = 3;
lat_edge = min(Lat)+(m_g*n_array)/2:m_g*n_array:max(Lat)+(m_g*n_array)/2;
lon_edge = min(Lon)+(m_g*n_array)/2:m_g*n_array:max(Lon)+(m_g*n_array)/2;
meanData = histcn([Lat(:) Lon(:)], lat_edge, lon_edge, 'AccumData', Data, 'Fun', @mean);
meanData(meanData==0)=NaN;

pico = struct();

%Modelized data

pico.biom_mod = meanData;
pico.lon_mod = lon_edge;
pico.lat_mod = lat_edge;

%In situ data

pico.biom_situ = pico_biom;
pico.lon_situ = lon_abundance(index_hipp);
pico.lat_situ = lat_abundance(index_hipp);

%Interp modelized data on the mesh grid

pico.lon_interp = min(pico.lon_mod):m_g:max(pico.lon_mod);
pico.lat_interp = min(pico.lat_mod):m_g:max(pico.lat_mod);
[Xq, Yq] = meshgrid(pico.lon_interp,pico.lat_interp);
pico.biom_interp = interp2(pico.lon_mod(1:end-1),pico.lat_mod(1:end-1),pico.biom_mod,Xq, Yq);

save('outputs/pico_struc','pico');

%%%%%%%%%%%%%% FIGURES %%%%%%%%%%%%%%

figure('DefaultAxesFontSize',22)
m_proj('mercator','lon',[lon_min lon_max],'lat',[lat_min lat_max]);
m_pcolor(Xq, Yq, pico.biom_interp)
shading flat
hold on
m_scatter(lon_abundance(index_hipp),lat_abundance(index_hipp),30,pico_biom,'filled');
m_quiver(X_curr,Y_curr,mean(u0(:,:,1:13),3),mean(v0(:,:,1:13),3),3,'r');
m_usercoast('gumby','patch','w');
m_grid('box','fancy','linestyle','-','gridcolor','none','backcolor','none');
hold on
cbar = colorbar;
%caxis([min(pico_biom) max(pico_biom)]);
caxis([0.435 0.455])
colorTitleHandle = get(cbar,'Title');
titleString = ({'Biomasse (mmol C.m⁻³)'});
set(colorTitleHandle ,'String',titleString);

%%%%%%%%%%%%%% Micro %%%%%%%%%%%%%%

%convert modeled biomasse mmol C/m3 into abundance (cell/cm3)

QC_micro = 0.26*exp(-5.8702)*(90*10000)^(0.9228)*0.86*1000;
QC_micro = QC_micro*1E-12/12.106;

zoo_ab = struct();
zoo_ab.P_big = (zoo_all.P_big./QC_micro)./1000000; %convert biomasse mmol C/m3 into abundance (cell/cm3)

%convert in situ abundances cell/cm³ into biomasse (mmol C/m3)

micro_ab = dataCYTONEW(index_hipp,14);
micro_ab = micro_ab.*1000000; %cell/m3
micro_biom = micro_ab.*QC_micro;

Data = zoo_all.P_big(:,dt_1:dt_end);

n_array = 3;
lat_edge = min(Lat)+(m_g*n_array)/2:m_g*n_array:max(Lat)+(m_g*n_array)/2;
lon_edge = min(Lon)+(m_g*n_array)/2:m_g*n_array:max(Lon)+(m_g*n_array)/2;
meanData = histcn([Lat(:) Lon(:)], lat_edge, lon_edge, 'AccumData', Data, 'Fun', @mean);
meanData(meanData==0)=NaN;

micro = struct();

%Modelized data

micro.biom_mod = meanData;
micro.lon_mod = lon_edge;
micro.lat_mod = lat_edge;

%In situ data

micro.biom_situ = micro_biom;
micro.lon_situ = lon_abundance(index_hipp);
micro.lat_situ = lat_abundance(index_hipp);

%Interp modelized data on the mesh grid

micro.lon_interp = min(micro.lon_mod):m_g:max(micro.lon_mod);
micro.lat_interp = min(micro.lat_mod):m_g:max(micro.lat_mod);
[Xq, Yq] = meshgrid(micro.lon_interp,micro.lat_interp);
micro_interp = interp2(micro.lon_mod(1:end-1),micro.lat_mod(1:end-1),micro.biom_mod,Xq, Yq);
micro.biom_interp = micro_interp;

save('outputs/micro_struc','micro');

%%%%%%%%%%%%%% FIGURES %%%%%%%%%%%%%%

figure('DefaultAxesFontSize',22)
m_proj('mercator','lon',[lon_min lon_max],'lat',[lat_min lat_max]);
m_pcolor(Xq, Yq, micro.biom_interp)
shading flat
hold on
m_scatter(lon_abundance(index_hipp),lat_abundance(index_hipp),30,micro_biom,'filled');
m_quiver(X_curr,Y_curr,mean(u0(:,:,1:13),3,'omitnan'),mean(v0(:,:,1:13),3,'omitnan'),3,'r');
m_usercoast('gumby','patch','w');
m_grid('box','fancy','linestyle','-','gridcolor','none','backcolor','none');
hold on
cbar = colorbar;
%caxis([min(micro_biom) max(micro_biom)]);
caxis([0.125 0.14]);
colorTitleHandle = get(cbar,'Title');
titleString = ({'Biomasse (mmol C.m⁻³)'});
set(colorTitleHandle ,'String',titleString);

disp('start_GA_toolbox : done');

toc