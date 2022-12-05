%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% INITIALIZE SST DATA AND CREATE SST MAP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% N. KIENTZ, 2022 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%clear all;
close all;
disp('sst_map : ...')
% Add needed utility functions. Remember to add the path where histcn is installed.
addpath(genpath('utils'))
addpath(genpath('outputs'))
addpath(genpath('SST'))


%% IMPORT SST DATA
list = dir('SST/nrt*.mat');
%%% PARAMETERS
l = 0.1; %to enlarge axis
file_start=1; %sst file start
file_end =21; %sst file end

%figure;
for i = file_start:file_end
    files_SST = list(i).name; %get group names
    D=load(files_SST);
    [lon,lat,sst]=deal(D.lon,D.lat,D.sst); %get lat and lon indices of chl files
    sst_file(:,:,i)=sst; %matrix with data from each chl file
    
%     subplot(3,7,i)
%     m_proj('mercator','lon',[min(lon) max(lon)],'lat',[min(lat) max(lat)]);
%     m_pcolor(lon,lat,sst)
%     cbar = colorbar('eastoutside');
%     cbar.Title.String = {'SST','[°C]'};
%     m_coast('patch',[.9 .9 .9],'edgecolor','none'); 
%     m_grid('tickdir','out','yaxislocation','left',...
%             'xaxislocation','bottom','xlabeldir','end','ticklen',.02);
%     title(i);
end

%sst_focus = sst_file(:,:,13);

averaged_sst = nanmean(sst_file(:,:,13:14),3); %mean of chl accross every selected files

% figure;
% for i = 13
%     files_SST = list(i).name; %get group names
%     D=load(files_SST);
%     [lon,lat,chl]=deal(D.lon,D.lat,D.sst); %get lat and lon indices of chl files
%     sst_file(:,:,i)=sst; %matrix with data from each chl file
%     
% %     subplot(1,3,i-12)
%     m_pcolor(lon,lat,averaged_sst)
%     hold on
%      m_proj('mercator','lon',[min(lon) max(lon)],'lat',[min(lat) max(lat)]);
%     m_gshhs_h('save','gumby');
%     m_usercoast('gumby','patch','w'); 
%     m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor','none');
%     hold on
%     m_plot(dataCYTONEW(407:511,end-1),dataCYTONEW(407:511,end),'.r');
%     hold on
%     cbar = colorbar
%     %caxis([0.1 0.6]);
%     title('Averaged sst between 12 and 13 may');
%     colorTitleHandle = get(cbar,'Title');
%     titleString = ({'SST','[°C]'});
%     set(colorTitleHandle ,'String',titleString);
% end

% averaged_sst = nanmean(sst_file,3); %mean of chl accross every selected files
% 
% %% Standard deviation for SST
% square_dist=[];
% N = ones(199,350)*(file_end-file_start); %199 = length(lat) and 350 = length(lon)
% SD = [];
% summed_N = [];
% for it = file_start:file_end
%             N = N - isnan(sst_file(:,:,it));
% end
% for iy = 1:length(lat)
%     for ix = 1:length(lon)
%         for it = file_start:file_end
%             square_dist(iy,ix,it) = abs((sst_file(iy,ix,it)-averaged_sst(iy,ix))).^2;
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
% title('Standard deviation for SST');
%     
% subplot(1,2,2)
%     m_proj('mercator','lon',[min(lon)-l max(lon)+l],'lat',[min(lat)-l max(lat)+l]);
%     m_gshhs_h('save','gumby');
%     m_usercoast('gumby','patch','w'); 
%     m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor','none');
%     hold on
% m_pcolor(lon,lat,N);
% colorbar;
% caxis([0 21]);
% title('Number of measurement for SST along the period');

disp('sst_map : done')