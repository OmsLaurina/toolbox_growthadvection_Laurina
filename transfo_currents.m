%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% CREATE FRIENDLY CURRENT MATRIX FOR ARIANE2D %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% N. KIENTZ 2022 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% a 3D matrix is needed with longitude, latitude and time for dimensions
%%% as 'lat', 'lon', 'time' and a 'curr' variable as the complexe of uo+ivo
%%% YOU NEED TO RUN IT BEFORE start_GA_toolbox.m
clear all;
close all;
addpath(genpath('SST'))

disp('transfo current ...')
%%%%%% IMPORT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%file = 'cmems_obs-sl_eur_phy-ssh_my_allsat-l4-duacs-0.125deg_P1D_1642067031887.nc';
%file = 'currents_data/cmems_obs-sl_eur_phy-ssh_my_allsat-l4-duacs-0.125deg_P1D_1646929918181.nc';
file_curr = 'cmems_obs-sl_eur_phy-ssh_my_allsat-l4-duacs-0.125deg_P1D_1647597462795.nc';
%ncinfo(file) %display info of the structure of the file (e.g., total dimension, # of variables, ...)
%ncdisp(file) % %display info of the file (e.g., variables, dimensions, ...)

%%%%%% EXTRACTING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lon = double(ncread(file_curr,'lon_bnds')) ;
lon = lon(1:end);
%lon(end+1) = 6;
lat = double(ncread(file_curr,'lat_bnds')) ;
lat = lat(1:end);

time = double(ncread(file_curr,'time')) ;
time = time(1:end);
%to avoid repetition of lon
lon = lon(lon>=0);
lon_ok = lon(1:2:length(lon));
%lon_ok = lon_ok(lon_ok <=7);
lat_ok = lat(1:2:length(lat));
%lat_ok = lat_ok(lat_ok >= 36 & lat_ok <=40);
iter_lat=[];
for i = 1:length(lat_ok)
    if lat_ok(i) >= 36 && lat_ok(i) <=40
        iter_lat(length(iter_lat)+1)=i;
    end
end
iter_lon=[];
for i = 1:length(lon_ok)
    if lon_ok(i) >= 0 && lon_ok(i) <= 7
        iter_lon(length(iter_lon)+1)=i;
    end
end
% lat = lat_ok(lat_ok>=36.02 & lat_ok<=39.98);
% lon = lon_ok(lon_ok>=0 & lon_ok<=6.98);

u0 = ncread(file_curr,'ugos') ;
u0 = u0(iter_lon,iter_lat,:);
u0 = permute(u0,[2 1 3]);
v0= ncread(file_curr,'vgos') ;
v0 = v0(iter_lon,iter_lat,:);
v0 = permute(v0,[2 1 3]);


%%%%%% FORMATING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curr=struct('lat',lat_ok,'lon',lon_ok,'time',time,'curr',complex(u0,v0));
curr = permute(curr,[2 1]);

disp('transfo current : done')