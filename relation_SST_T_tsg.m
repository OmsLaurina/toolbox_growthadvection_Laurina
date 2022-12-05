%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% INITIALIZE CHL DATA AND CREATE CHL MAP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% N. KIENTZ, 2022 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%clear all;
close all;

% Add needed utility functions. Remember to add the path where histcn is installed.
addpath(genpath('utils'))
addpath(genpath('outputs'))
addpath(genpath('TEMP'))

dataCYTONEW = dlmread('data_CYTO_NEW.txt');
file = importdata('inputs/data_CYTO.txt');
file = file.textdata;

sst_MUR_file = dir('SST/MUR_SST.nc');
sst_MUR = double(ncread(sst_MUR_file.name,'analysed_sst'));
% NOAA_2_n = double(ncread(NOAA_files_n(2).name,'n_an'));
lat_MUR= double(ncread(sst_MUR_file.name,'latitude'));
lon_MUR= double(ncread(sst_MUR_file.name,'longitude'));
time_MUR = double(ncread(sst_MUR_file.name,'time'));
STRUC_MUR = struct();
STRUC_MUR.lat = lat_MUR;
STRUC_MUR.lon = lon_MUR;
STRUC_MUR.sst = sst_MUR;
STRUC_MUR.time = time_MUR;
save('STRUC_MUR','STRUC_MUR')

% sst_CMEMS_file = dir('SST/SST*.nc');
% sst_CMEMS = double(ncread(sst_CMEMS_file.name,'analysed_sst'));
% lat_CMEMS= double(ncread(sst_CMEMS_file.name,'lat'));
% lon_CMEMS= double(ncread(sst_CMEMS_file.name,'lon'));


% list_sst = dir('SST/nrt*.mat');
% file_start=1; %sst file start
% file_end =21; %sst file end
% sst_file = [];
% for i = file_start:file_end
%     files_SST = list_sst(i).name; %get group names
%     D=load(files_SST);
%     [lon_sst,lat_sst,sst]=deal(D.lon,D.lat,D.sst); 
%     sst_file(:,:,i)=sst; %matrix with data from each sst file
% end
T = dataCYTONEW(:,end-3);


%% SST in function of abundance
lat_sst=[];
lon_sst=[];
sst_ab=[];
sst_from_averaged=[];
start = 1;
stop = 548;

% for i = 1:length(dataCYTONEW)
% date = dataCYTONEW(i,3);
%     if date == 30 
%         sst_data = sst_file(:,:,2);
%     else
%         sst_data = sst_file(:,:,date+2);
%     end   
%         [d, iy ] = min( abs( D.lat-dataCYTONEW(i,end)) );         
%         lat_sst(i)=iy;%indice des latitudes
%         [d, ix ] = min( abs( D.lon-dataCYTONEW(i,end-1)));
%         lon_sst(i)=ix;%indice des longitudes
%         sst_ab(i)=sst_data(lat_sst(i),lon_sst(i));
%         if isnan(sst_data(lat_sst(i),lon_sst(i)))
%              sst_ab(i) = NaN;
%         else
%             sst_ab(i) = sst_data(lat_sst(i),lon_sst(i));
%         end
% end



%% SST_MUR in function of T_tsg
lat_sst_MUR=[];
lon_sst_MUR=[];
sst_t_MUR=[];

for i = 1:length(dataCYTONEW)
    date = dataCYTONEW(i,3);
    [d, iy ] = min( abs( lat_MUR-dataCYTONEW(i,end) )); 
    lat_sst_MUR(i)=iy;%indice des latitudes
    [d, ix ] = min( abs( lon_MUR-dataCYTONEW(i,end-1)));
    lon_sst_MUR(i)=ix;%indice des longitudes
%     sst_t(i)=sst_MUR(lon_sst_clim(i),lat_sst_clim(i),1);
    if date == 30 
            sst_t_MUR(i)=sst_MUR(lon_sst_MUR(i),lat_sst_MUR(i),120);
        else 
            sst_t_MUR(i)=sst_MUR(lon_sst_MUR(i),lat_sst_MUR(i),120+date);
        end
       
end

mdl = fitlm(T,sst_t_MUR );

% figure;
% plotAdded(mdl)
% ylabel('SST MUR')
% xlabel('T\_tsg')
% title(mdl.Rsquared.Ordinary)

% m_proj('mercator','lon',[min(X(1,:)) max(X(1,:))],'lat',[min(Y(:,1)) max(Y(:,1))]);
% hold on
% m_pcolor(X(1,:),Y(:,1),mat_nano(:,:,2))
% shading flat
% hold on
% m_gshhs_h('save','gumby');
% m_usercoast('gumby','patch','w'); 
% m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor','none');
% hold on

%% SST_CMEMS in function of T_tsg
% lat_sst_CMEMS=[];
% lon_sst_CMEMS=[];
% sst_t_CMEMS=[];
% 
% for i = 1:length(dataCYTONEW)
%     date = dataCYTONEW(i,3);
%     [d, iy ] = min( abs( lat_CMEMS-dataCYTONEW(i,end) )); 
%     lat_sst_CMEMS(i)=iy;%indice des latitudes
%     [d, ix ] = min( abs( lon_CMEMS-dataCYTONEW(i,end-1)));
%     lon_sst_CMEMS(i)=ix;%indice des longitudes
% %     sst_t(i)=sst_MUR(lon_sst_clim(i),lat_sst_clim(i),1);
%     if date == 30 
%             sst_t_CMEMS(i)=sst_CMEMS(lon_sst_CMEMS(i),lat_sst_CMEMS(i),120)-273.15;
%         else 
%             sst_t_CMEMS(i)=sst_CMEMS(lon_sst_CMEMS(i),lat_sst_CMEMS(i),120+date)-273.15;
%         end
%        
% end
% 
% mdl = fitlm(T,sst_t_CMEMS );
% 
% figure;
% plotAdded(mdl)
% ylabel('SST CMEMS')
% xlabel('T\_tsg')
% title(mdl.Rsquared.Ordinary)
% 
% mdl = fitlm(sst_t_MUR,sst_t_CMEMS);
% figure;
% plotAdded(mdl)
% ylabel('SST CMEMS')
% xlabel('SST MUR')
% title(mdl.Rsquared.Ordinary)