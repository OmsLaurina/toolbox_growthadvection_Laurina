%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% INITIALIZE ABUNDANCE DATA AND CREATE ABUNDANCE MAP %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% N. KIENTZ, 2022 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%clear all;
close all;
disp('abundance_map : ...')
% Add needed utility functions. Remember to add the path where histcn is installed.
addpath(genpath('utils'))
addpath(genpath('outputs'))

%% Import Abundance data
dataCYTONEW = dlmread('data_CYTO_NEW.txt');
file = importdata('inputs/data_CYTO.txt');
file = file.textdata;

%% SET ABUNDANCE DISTRIBUTION IN FUNCTION OF CHL DISTRIBUTION
%%% HAVE TO RUN chl_map FIRST
%%% find position of chl > 2 or < 2 (water type 1 and 2 refered to the
%%% front from Tzortzis et al. (2021))

% set names, abundances and empty matrix 
Z_name = {'syne','pico1', 'pico2'}; 
Z1_val = [22000 1100 1700];
Z2_val = [10000 500 900];
abundance_matrix=ones(length(lat),length(lon),length(Z_name));

% figure;
% hold on
for i = 1:length(Z_name)
    
    averaged_chl(averaged_chl==0)=NaN;
    z1 = averaged_chl >= 0.21429;
    z1 = z1*Z1_val(i);

    z2 = averaged_chl < 0.21429;
    z2 = z2*Z2_val(i);
        
    mat = z1+z2; 
    mat(mat == 0) = NaN;
    
    abundance_matrix(:,:,i)=mat;
    
    
% %     % Abundance patch distribution map
%     figure;
%     %subplot(abs(4-length(Z_name)),3,i)
%   
%     m_proj('mercator','lon',[min(lon)-l max(lon)+l],'lat',[min(lat)-l max(lat)+l]);
%      m_gshhs_h('save','gumby');
%     m_usercoast('gumby','patch','w'); 
%     %m_gshhs_i('color','k');
%     %m_grid('linestyle','none','tickdir','out','linewidth',3);
%     m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor','none');
%     hold on
%     m_contourf(lon,lat,mat);
%     title(Z_name(i));
%     %subtitle(Z_name(i));
%     cbar = colorbar;
%     caxis([Z2_val(i) Z1_val(i)]);
%     colorTitleHandle = get(cbar,'Title');
%     titleString = ({'Abundance','[\it{cells.cm^{-3}}]',' '});
%     set(colorTitleHandle ,'String',titleString);
end

%abundance_matrix(end+1:200,:,:)=missing; %add an empty line to match dimension of other matrix (multiple of 20x35) i.e. 20 lat and 35 lon 

%% CONTOUR CHL X ABUNDANCE PLOT
%it=1;
% figure ;
% hold on
% for ip = 7:19
% it=it+1;
% subplot(2,7,it-1)
% m_proj('mercator','lon',[min(lon)-l max(lon)+l],'lat',[min(lat)-l max(lat)+l]);
% m_gshhs_h('save','gumby');
% m_usercoast('gumby','patch','w'); 
% m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor','none');
% hold on
% m_contour(lon,lat,averaged_chl,chl_level,'ShowText','on');
% hold on
% m_scatter(dataCYTONEW(:,end-1),dataCYTONEW(:,end),40,dataCYTONEW(:,ip),'filled');
% cbar = colorbar;
% colormap jet
% %caxis([10 10000]);
% title(file(1,it));
% colorTitleHandle = get(cbar,'Title');
% titleString = ({'Abundance','[ \it {cells.cm^{-3}} ]',''});
% set(colorTitleHandle ,'String',titleString);
% end

disp('abundance_map : done')