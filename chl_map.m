%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% INITIALIZE CHL DATA AND CREATE CHL MAP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% N. KIENTZ, 2022 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%clear all;
close all;
disp('chl_map : ...')
% Add needed utility functions. Remember to add the path where histcn is installed.
addpath(genpath('utils'))
addpath(genpath('outputs'))
addpath(genpath('CHL'))

dataCYTONEW = dlmread('data_CYTO_NEW.txt');
file = importdata('inputs/data_CYTO.txt');
file = file.textdata;

%% IMPORT CHL DATA 
list = dir('CHL/nrt*.mat');

%%% PLOT PARAMETERS
chl_level = linspace(0.1,0.5,8); %set the chl min, max and number of indent for contourf
l = 0.1; %to enlarge axis
file_start=1; %chl file start
file_end =21; %chl file end

%figure;
for i = file_start:file_end
    files_CHL = list(i).name; %get group names
    D=load(files_CHL);
    [lon,lat,chl]=deal(D.lon,D.lat,D.chl); %get lat and lon indices of chl files
    chl_file(:,:,i)=chl; %matrix with data from each chl file
    
%     subplot(1,3,i-12)
%     m_pcolor(lon,lat,chl)
%     hold on
%      m_proj('mercator','lon',[min(lon) max(lon)],'lat',[min(lat) max(lat)]);
%     m_gshhs_h('save','gumby');
%     m_usercoast('gumby','patch','w'); 
%     m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor','none');
%     hold on
%     m_plot(dataCYTONEW(407:511,end-1),dataCYTONEW(407:511,end),'.r');
%     hold on
%     cbar = colorbar
%     caxis([0.1 0.6]);
%     title('11 may CHL map');
%     colorTitleHandle = get(cbar,'Title');
%     titleString = ({'Chl','[mg m-3]'});
%     set(colorTitleHandle ,'String',titleString);
end

averaged_chl = nanmean(chl_file,3); %mean of chl accross every selected files

%% CHL MAP PLOTS
%%% PCOLOR PLOT OF AVERAGED CHL
% figure;
% m_proj('mercator','lon',[min(lon) max(lon)],'lat',[min(lat) max(lat)]);
% m_coast('patch',[.9 .9 .9],'edgecolor','none'); 
% m_grid('tickdir','out','yaxislocation','left',...
%             'xaxislocation','bottom','xlabeldir','end','ticklen',.02);
% hold on
% m_pcolor(lon,lat,averaged_chl);
% shading flat
% caxis([min(chl_level) max(chl_level)]);
% cbar = colorbar;
% title('Averaged Chl between 2018-04-28 and 2018-05-20');
% colorTitleHandle = get(cbar,'Title');
% titleString = ({'Chl','[ \it{µg.L^{-1}} ]',''});
% set(colorTitleHandle ,'String',titleString);

%%% CONTOURF PLOT OF AVERAGED CHL
% figure;
% m_proj('mercator','lon',[min(lon)-l max(lon)+l],'lat',[min(lat)-l max(lat)+l]);
% m_gshhs_h('save','gumby');
% m_usercoast('gumby','patch','w'); 
% m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor','none');
% hold on
% m_contourf(lon,lat,averaged_chl,chl_level,'ShowText','on');
% hold on
% scatter(dataCYTONEW(:,end-1),dataCYTONEW(:,end),dataCYTONEW(:,10))
% caxis([min(chl_level) max(chl_level)]);
% cbar = colorbar;
% title('Averaged Chl between 2018-04-28 and 2018-05-20');
% %subtitle(file(1,4));
% colorTitleHandle = get(cbar,'Title');
% titleString = ({'Chl','[ \it{µg.L^{-1}} ]'});
% set(colorTitleHandle ,'String',titleString);

%% Standard deviation for CHL
% square_dist=[];
% N = ones(199,350)*(file_end-file_start);
% SD = [];
% summed_N = [];
% for it = file_start:file_end
%             N = N - isnan(chl_file(:,:,it));
% end
% for iy = 1:length(lat)
%     for ix = 1:length(lon)
%         for it = file_start:file_end
%             square_dist(iy,ix,it) = abs((chl_file(iy,ix,it)-averaged_chl(iy,ix))).^2;
%         end
%         
%         summed_N(iy,ix) = sum(N(iy,ix,:));
%         SD(iy,ix) = sqrt(nansum(square_dist(iy,ix,:))/summed_N(iy,ix));      
%     end
% end
% SD(SD==0)= NaN; % for better looking map
% N(N==0)= NaN; % same

%% SD MAP PLOT

% figure, hold on
%     subplot(1,2,1)
%     m_proj('mercator','lon',[min(lon)-l max(lon)+l],'lat',[min(lat)-l max(lat)+l]);
%     m_gshhs_h('save','gumby');
%     m_usercoast('gumby','patch','w'); 
%     m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor','none');
%     hold on
% m_pcolor(lon,lat,SD);
% colorbar;
% caxis([0 1]);
% title('Standard deviation for chl');
    
% subplot(1,2,2)
%     m_proj('mercator','lon',[min(lon)-l max(lon)+l],'lat',[min(lat)-l max(lat)+l]);
%     m_gshhs_h('save','gumby');
%     m_usercoast('gumby','patch','w'); 
%     m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor','none');
%     hold on
% m_pcolor(lon,lat,N);
% colorbar;
% caxis([0 21]);
% title('Number of measurement for chl along the period');

disp('chl_map : done')