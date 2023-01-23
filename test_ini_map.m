%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% INITIALIZE ABUNDANCE DATA AND CREATE ABUNDANCE MAP %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% N. KIENTZ, 2022 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
disp('test_ini_map : ...');
%clear all;
%close all;

% Add needed utility functions. Remember to add the path where histcn is installed.
% addpath(genpath('utils'))
% addpath(genpath('outputs'))
% addpath(genpath('Ariane_workplace'))

%run('transfo_currents.m')
%run('relation_SST_T_tsg.m');
%%

% %initial sst matrix
% lon_shape_iter = (lon_MUR <= 7 & lon_MUR >= 0);
% lat_shape_iter = (lat_MUR <= 40 & lat_MUR >= 36);
% lon_shaped=[];
% lat_shaped=[];
% for i = 1:length(lon_shape_iter)
%     if lon_shape_iter(i)>=0
%         lon_shaped(length(lon_shaped)+1)=i;
%     end
% end
% for i = 1:length(lat_shape_iter)
%     if lat_shape_iter(i)>=0
%         lat_shaped(length(lat_shaped)+1)=i;
%     end
% end
% sst_shaped = sst_MUR(lon_shaped,lat_shaped,120);
% mat_sst = sst_shaped;
% 
% xvalue = min(lon):0.01:max(lon);
% yvalue = min(lat):0.01:max(lat);
% [X,Y] = meshgrid(xvalue, yvalue);
% point_1 = [X(:) Y(:)];


%% Import Abundance data
dataCYTONEW = dlmread('data_CYTO_NEW.txt');
file = importdata('inputs/data_CYTO.txt');
file = file.textdata;

%% Abundance at t0
lon_ini = [];
iter_lon = [];
lat_ini = [];
iter_lat = [];
ab_micro = [];
ab_pico = [];
ab_nano = [];

mat_micro = nan(length(yvalue),length(xvalue));
mat_pico = nan(length(yvalue),length(xvalue));
mat_nano = nan(length(yvalue),length(xvalue));

%initial abundance matrix
n=2;
for i = 1:n-1
    lon_ini(i)=dataCYTONEW(i,end-1);
    iter_lon(i) = floor((lon_ini(i)-min(lon_MUR))/0.01);
    lat_ini(i)=dataCYTONEW(i,end);
    iter_lat(i) = floor((lat_ini(i)-min(lat_MUR))/0.01);
    ab_micro(i) = dataCYTONEW(i,14);
    ab_pico(i) = dataCYTONEW(i,16);
    ab_nano(i) = dataCYTONEW(i,17);
    mat_micro(iter_lat(i),iter_lon(i)) = ab_micro(i);
    mat_pico(iter_lat(i),iter_lon(i)) = ab_pico(i);
    mat_nano(iter_lat(i),iter_lon(i)) = ab_nano(i);
end

%initial sst matrix
sst_files = dir('SST/MUR_SST.nc');
sst_MUR = double(ncread(sst_files.name,'analysed_sst'));
% NOAA_2_n = double(ncread(NOAA_files_n(2).name,'n_an'));
lat_MUR= double(ncread(sst_files.name,'latitude'));
lon_MUR= double(ncread(sst_files.name,'longitude'));
% lon_shape_iter = (lon_MUR <= max(lon) & lon_MUR >= min(lon));
% lat_shape_iter = (lat_MUR <= max(lat) & lat_MUR >= min(lat));
% lon_shaped=[];
% lat_shaped=[];
% for i = 1:length(lon_shape_iter)
%     if lon_shape_iter(i)>0
%         lon_shaped(length(lon_shaped)+1)=i;
%     end
% end
% for i = 1:length(lat_shape_iter)
%     if lon_shape_iter(i)>0
%         lat_shaped(length(lat_shaped)+1)=i;
%     end
% end

%sst_shaped = sst_MUR(lon_shaped,lat_shaped,120);

pente_t_tsg = 0.983139 ;

% for ix = 1:length(lon_shaped)
%     for iy = 1:length(lat_shaped)
% %micro
%             pente_micro = -7.08748*pente_t_tsg ;
%             delta_sst = sst_shaped(ix,iy)-sst_shaped(iter_lon,iter_lat);
%             extp_micro = mat_micro(iter_lat,iter_lon)+delta_sst*pente_micro;
%             mat_micro(iy,ix) = extp_micro;
% %pico
%             pente_pico = 390.01*pente_t_tsg;
%             delta_sst = sst_shaped(ix,iy)-sst_shaped(iter_lon,iter_lat);
%             extp_pico = mat_pico(iter_lat,iter_lon)+delta_sst*pente_pico;
%             mat_pico(iy,ix) = extp_pico;
%      end
% end

for ix = 1:length(lon_shaped)
    for iy = 1:length(lat_shaped)
%micro
            delta_sst = sst_shaped(ix,iy)-sst_shaped(iter_lon,iter_lat);
            pente_micro = (0.6467*(delta_sst)^4 - 2.682*(delta_sst)^3 + 2.28*(delta_sst)^2 - 2.1*delta_sst + 13.4)*pente_t_tsg ;
            extp_micro = pente_micro;
            if extp_micro >= 0
            mat_micro(iy,ix) = extp_micro;
            end

%pico
            delta_sst = sst_shaped(ix,iy)-sst_shaped(iter_lon,iter_lat);
            pente_pico = (24.51*(delta_sst)^4 + 53.51*(delta_sst)^3 + 24.56*(delta_sst)^2 + 260.3*delta_sst + 1022)*pente_t_tsg;
            extp_pico = pente_pico;
            if extp_pico >= 0
            mat_pico(iy,ix) = extp_pico;
            end
    end
end
%convert modeled biomasse mmol C/m3 into abundance (cell/cm3)
QC_pico = 0.26*exp(-5.8702)*(0.9*10000)^(0.9228)*0.86*1000;
QC_pico = QC_pico*1E-12/12.106;

%convert in situ abundances cell/cm³ into biomasse (mmol C/m3)
pico_ab = mat_pico(:,:);
pico_ab = pico_ab.*1000000; %cell/m3
pico_biom = pico_ab.*QC_pico;

%convert modeled biomasse mmol C/m3 into abundance (cell/cm3)
QC_micro = 0.26*exp(-5.8702)*(90*10000)^(0.9228)*0.86*1000;
QC_micro = QC_micro*1E-12/12.106;

%convert in situ abundances cell/cm³ into biomasse (mmol C/m3)
micro_ab = mat_micro(:,:);
micro_ab = micro_ab.*1000000; %cell/m3
micro_biom = micro_ab.*QC_micro;

micro_biom(micro_biom<=0)=NaN;
pico_biom(pico_biom<=0)=NaN;

biom_ini = struct();
biom_ini.micro = micro_biom;
biom_ini.pico = pico_biom;
save('inputs/biom_ini','biom_ini')

% figure('DefaultAxesFontSize',22);
% m_proj('mercator','lon',[0 7],'lat',[36 40]);
% hold on
% m_pcolor(min(lon_MUR):0.01:max(lon_MUR),min(lat_MUR):0.01:max(lat_MUR),pico_biom)
% shading flat
% hold on
% %m_gshhs_h('save','gumby');
% m_usercoast('gumby','patch','w'); 
% m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor','none');
% hold on
% cbar = colorbar;
% caxis([0.15 0.35]);
% colorTitleHandle = get(cbar,'Title');
% titleString = ({'Biomasse (mmol C.m⁻³)'});
% set(colorTitleHandle ,'String',titleString);
% 
% 
% figure('DefaultAxesFontSize',22);
% m_proj('mercator','lon',[0 7],'lat',[36 40]);
% hold on
% m_pcolor(lon_MUR,lat_MUR,micro_biom)
% shading flat
% hold on
% m_usercoast('gumby','patch','w'); 
% m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor','none');
% hold on
% cbar = colorbar;
% caxis([0.15 0.35]);
% colorTitleHandle = get(cbar,'Title');
% titleString = ({'Biomasse (mmol C.m⁻³)'});
% set(colorTitleHandle ,'String',titleString);

% figure;
% m_proj('mercator','lon',[0 7],'lat',[36 40]);
% hold on
% m_pcolor(lon_MUR,lat_MUR,sst_shaped')
% shading flat
% hold on
% m_gshhs_h('save','gumby');
% m_usercoast('gumby','patch','w'); 
% m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor','none');
% hold on
% colorbar

%%% initial chl_map
chl_files = dir('CHL/dataset*.nc');
chl_CMEMS = double(ncread(chl_files.name,'CHL'));
lat_CMEMS= double(ncread(sst_files.name,'latitude'));
lon_CMEMS= double(ncread(sst_files.name,'longitude'));

disp('test_ini_map : done')