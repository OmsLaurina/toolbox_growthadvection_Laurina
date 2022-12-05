%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% PLOTS RELATION BETWEEN ABUNDANCE AND PARAMETERS %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% (S, T, SST, T/S, CHL, PAR) %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% N. KIENTZ, 2022 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
disp('Abundance_relations : ...')
%%
%%% RUN chl_map AND sst_map AND abundance_map FIRST TO INITIALIZE DATA %%%
%
%% CHL in function of abundance
lat_abundance = []; %latitude of each measurement of abundance
lon_abundance = []; %longitude of each measurement of abundance
group_ab = [];
lat_chl=[];
lon_chl=[];
chl_ab=[];
delta_lat=[];
delta_lon=[];
chl_from_averaged=[];
lat_chl_from_averaged = [];
lon_chl_from_averaged = [];
averaged_chl_iter = [];
start=1;
stop = 548;
% figure;
% hold on
for n = 11
for i = 1:length(dataCYTONEW)
lat_abundance(i) = dataCYTONEW(i,end);
lon_abundance(i) = dataCYTONEW(i,end-1);
group_ab(n,i) = dataCYTONEW(i,6+n); 
date = dataCYTONEW(i,3);
    if date == 30 
        chl_data = chl_file(:,:,2);
    else
        chl_data = chl_file(:,:,date+2);
    end 
    
        [d, iy ] = min( abs( D.lat-lat_abundance(i) ) );        
        lat_chl(i)=iy;%indice des latitudes
        delta_lat(i) = D.lat(lat_chl(i))-lat_abundance(i);
        [d, ix ] = min( abs( D.lon-lon_abundance(i) ) );
        lon_chl(i)=ix;%indice des longitudes
        delta_lon(i) = D.lon(lon_chl(i))-lon_abundance(i);
        chl_ab(i)=chl_data(lat_chl(i),lon_chl(i));
        if isnan(chl_data(lat_chl(i),lon_chl(i)));
            %averaged_chl_iter(length(averaged_chl_iter)+1) = i;
            chl_ab(i) = NaN; %replace meaned data by nan when initially no
            %chl_ab(i) = averaged_chl(lat_chl(i),lon_chl(i));
            %disp('o')
            %chl_from_averaged(length(chl_from_averaged)+1) = averaged_chl(lat_chl(i),lon_chl(i));
            %lat_chl_from_averaged(length(lat_chl_from_averaged)+1) = lat_chl(i);
            %lon_chl_from_averaged(length(lon_chl_from_averaged)+1)= lon_chl(i);
            %disp('averaged data used at iteration : ')
            %disp(i)
        else
            chl_ab(i) = chl_data(lat_chl(i),lon_chl(i));
        end
end
nano = group_ab(n,:)';
%mdl = fitlm(chl_ab(:),group_ab(n,:));

%     figure;
%     plotAdded(mdl);
%     hold on
%     %scatter(chl_ab(start:stop),group_ab(n,start:stop),30,datenum(dataCYTONEW(start:stop,1),dataCYTONEW(start:stop,2), dataCYTONEW(start:stop,3), dataCYTONEW(start:stop,4), dataCYTONEW(start:stop,5), dataCYTONEW(start:stop,6)),'filled');
%     colorbar;
%     cbdate;
%     grid on
%     grid minor
%     %plot(chl_ab(1:548),group_ab(n,1:501),'or') %%% hippodrome from 407 to 511
%     %hold on
%     %plot(chl_from_averaged(:),group_ab(n,averaged_chl_iter(:)),'.g');
%     ylabel('Abundance (cell.cm-3)')
%     xlabel('Chl (µg/L)')
%     title([file(1,n+1) mdl.Rsquared.Ordinary]);
%     subtitle('R² for period, replacement with 5 days mean')
 end


%% SST in function of abundance
% lat_sst=[];
% lon_sst=[];
% sst_ab=[];
% sst_from_averaged=[];
% start = 319;
% stop = 548;
% %figure;
% %hold on
% m = 0;
% for n = [8]
%     m=m+1;
% for i = 1:length(dataCYTONEW)
% group_ab(n,i) = dataCYTONEW(i,6+n); 
% date = dataCYTONEW(i,3);
%     if date == 30 
%         sst_data = sst_file(:,:,2);
%     else
%         sst_data = sst_file(:,:,date+2);
%     end   
%     
%         [d, iy ] = min( abs( D.lat-dataCYTONEW(i,end)) );         
%         lat_sst(i)=iy;%indice des latitudes
%         %delta_lat_sst(i) = D.lat(lat_chl(i))-lat_abundance(i);
%         [d, ix ] = min( abs( D.lon-dataCYTONEW(i,end-1)));
%         lon_sst(i)=ix;%indice des longitudes
%         %delta_lon_sst(i) = D.lon(lon_sst(i))-lon_abundance(i);
%         sst_ab(i)=sst_data(lat_sst(i),lon_sst(i));
%         if isnan(sst_data(lat_sst(i),lon_sst(i)))
%             %averaged_sst_iter(length(averaged_sst_iter)+1) = i;
%             %sst_ab(i) = averaged_sst(lat_sst(i),lon_sst(i));
%              sst_ab(i) = NaN;
% %             sst_from_averaged(length(sst_from_averaged)+1) = averaged_sst(lat_sst(i),lon_sst(i));
% %             lat_sst_from_averaged(length(lat_sst_from_averaged)+1) = lat_sst(i);
% %             lon_sst_from_averaged(length(lon_sst_from_averaged)+1)= lon_sst(i);
%         else
%             sst_ab(i) = sst_data(lat_sst(i),lon_sst(i));
%         end
%  end
%  mdl = fitlm(sst_ab(start:stop),group_ab(n,start:stop));
% % R2_matrix(n,1)=mdl.Rsquared.Ordinary;
% %     %subplot(3,1,m)
%     figure;
%     plot(mdl);
%     hold on
%     grid on
% grid minor
%     scatter(sst_ab(start:stop),group_ab(n,start:stop),30,datenum(dataCYTONEW(start:stop,1),dataCYTONEW(start:stop,2), dataCYTONEW(start:stop,3), dataCYTONEW(start:stop,4), dataCYTONEW(start:stop,5), dataCYTONEW(start:stop,6)),'filled');
%     colorbar;
%     cbdate;
%     %hold on
%     %plot(x1,p1,'r');
%     %hold on
%     %plot(sst_ab(:),group_ab(n,:),'or')
%     %hold on
%     %plot(sst_from_averaged(:),group_ab(n,averaged_sst_iter(:)),'.g');
%     ylabel('Abundance (cell.cm-3)')
%     xlabel('SST (°C)')
%     title([file(1,n+1) mdl.Rsquared.Ordinary]);
%     subtitle('R² for period, no replacement')
 %end
% 
%% T_tsg in function of abundance
T = dataCYTONEW(:,end-3);
S = dataCYTONEW(:,end-4);
ts=[];

for n = [4,8, 10]
for i = 1:length(dataCYTONEW)
group_ab(n,i) = dataCYTONEW(i,6+n); 
ts(i) = T(i)/S(i);
end
%mdl = fitlm(T(:),group_ab(n,:));

%     figure;
%     plotAdded(mdl)
%     hold on
%     grid('on')
%     grid('minor')
%     ylabel('Abundance (cell.cm3)')
%     xlabel('T')
%     title([file(1,n+1) mdl.Rsquared.Ordinary]);
end

pico2 = group_ab(4,:)';
micro = group_ab(8,:)';
pico = group_ab(10,:)';
% x = T;
% y = micro;
% % Régression linéaire: y = a*x +b
% nb_coefs=2;
% X=[x , ones(size(x))];
% coef_poly=X\y;  % tel que a=coef_poly(1); b=coef_poly(2);
% mat_yfit=repmat(coef_poly',length(x),1);
% yfit=sum(mat_yfit.*X,2);
% % r obtained by correlating yfit and y.
% 
% % Régression exp: y = a*exp(x) + b
% % idem en remplaçant la ligne X=[x , ones(size(x))]; par 
% X=[exp(x) , ones(size(x))];
% coef_poly_exp=X\y;  % tel que a=coef_poly(1); b=coef_poly(2);
% mat_yfit_exp=repmat(coef_poly_exp',length(x),1);
% yfit_exp=sum(mat_yfit_exp.*X,2);
% 
% x = mldivide(X,y);

% subplot(3,4,12)
% plot(S(407:511),T(407:511),'or')
% hold on
% plot(S(start:stop),T(start:stop),'.b')
% xlabel('Salinity')
% ylabel('Temperature')
% title('TS diagram') 
% legend('hippodrome','total')
% % 

 %% Focus on hippodrome in an area of 1°x1° around the hippodrome
%between 2.5 and 3.5 for lon and 38 and 39 for lat
% index_focus=[];
% for i = 1:length(dataCYTONEW) %detect index of values in the focus area
%     if (38 < dataCYTONEW(i,end) && 39 > dataCYTONEW(i,end) && 2.5 < dataCYTONEW(i,end-1) && dataCYTONEW(i,end-1) <3.5)
%         index_focus(length(index_focus)+1)=i;  
%     end
% end
% lat_focus = dataCYTONEW(:,end);
% lon_focus = dataCYTONEW(:,end-1);
% 
% %% CHL
% group_ab = [];
% lat_chl=[];
% lon_chl=[];
% chl_ab=[];
% delta_lat=[];
% delta_lon=[];
% 
% %
% chl_from_averaged=[];
% lat_chl_from_averaged = [];
% lon_chl_from_averaged = [];
% 
% %
% chl_from_pre_averaged=[];
% lat_chl_from_pre_averaged = [];
% lon_chl_from_pre_averaged = [];
% averaged_pre_chl_iter = [];
% date_pre_ab = [];
% start=1;
% stop = 548;
% % figure;
% % hold on
%  for n = 1:11
% %     averaged_chl_iter = [];
% %     chl_iter = [];
% % %for i = 1:length(index_focus) %%% to only have data of focus area
%  for i = 1:548 %%% to have all the data (have to be changed for the section "CHL relation with focus on the time evolution of the data")
% % %group_ab(n,i) = dataCYTONEW(index_focus(i),6+n);
% % %chl_mat_ab_time(i,1,n) = datenum(dataCYTONEW(i,1),dataCYTONEW(i,2), dataCYTONEW(i,3), dataCYTONEW(i,4), dataCYTONEW(i,5), dataCYTONEW(i,6));
% % group_ab(n,i) = dataCYTONEW(i,6+n); %replace i by index_focus(i) to only focus area
%  date = dataCYTONEW(i,3);
% % date_pre_ab(i) = datenum(dataCYTONEW(1,1),dataCYTONEW(1,2),dataCYTONEW(1,3),dataCYTONEW(1,4),dataCYTONEW(1,5),dataCYTONEW(1,6));
%     if date == 30 
%         chl_data = chl_file(:,:,2);
%         averaged_pre_chl = nanmean(chl_file(:,:,1:3),3); %mean of chl for the 30th april (-1/+1 from the date)
%     else
%         chl_data = chl_file(:,:,date+2);
%         averaged_pre_chl = nanmean(chl_file(:,:,date+1:date+3),3); %mean of chl (-1 / +1 from the date)
%     end   
%     
%         [d, iy ] = min( abs( D.lat-dataCYTONEW(i,end)) );        
%         lat_chl(i)=iy;%indice des latitudes
%         %delta_lat(i) = D.lat(lat_chl(i))-lat_focus(i);
%         [d, ix ] = min( abs( D.lon-dataCYTONEW(i,end-1)));
%         lon_chl(i)=ix;%indice des longitudes
%         %delta_lon(i) = D.lon(lon_chl(i))-lon_focus(i);
%         %chl_ab(i)=chl_data(lat_chl(i),lon_chl(i));
%         if isnan(chl_data(lat_chl(i),lon_chl(i)))
% %             
% %             % averaged data for all the period
% % %             averaged_chl_iter(length(averaged_chl_iter)+1) = i;
% % %             chl_ab(i) = averaged_chl(lat_chl(i),lon_chl(i));
% % %             chl_from_averaged(length(chl_from_averaged)+1) = averaged_chl(lat_chl(i),lon_chl(i));
% % %             lat_chl_from_averaged(length(lat_chl_from_averaged)+1) = lat_chl(i);
% % %             lon_chl_from_averaged(length(lon_chl_from_averaged)+1)= lon_chl(i);
% % 
% %             % more precise averaged data (-1/+1 from the day)
% %             averaged_pre_chl_iter(length(averaged_pre_chl_iter)+1) = i;
% %             %chl_ab(i) = averaged_pre_chl(lat_chl(i),lon_chl(i));
%              chl_ab(i) = NaN;
% %             chl_from_pre_averaged(length(chl_from_pre_averaged)+1) = averaged_pre_chl(lat_chl(i),lon_chl(i));
% %             lat_chl_from_pre_averaged(length(lat_chl_from_pre_averaged)+1) = lat_chl(i);
% %             lon_chl_from_pre_averaged(length(lon_chl_from_pre_averaged)+1)= lon_chl(i);
%          else
%              chl_ab(i) = chl_data(lat_chl(i),lon_chl(i));
% %             chl_iter(length(chl_iter)+1) = i;
%          end
%  end
% % % % mdl = fitlm(chl_ab(:),group_ab(n,:));
% % % % R2_matrix(n,1)=mdl.Rsquared.Ordinary;
% % % %     subplot(3,4,n)
% % % %     %plot(chl_ab,group_ab(n,index_focus),'.')
% % % %     %hold on
% % % %     plot(mdl);
% % % % %     hold on
% % % % %     plot(chl_from_pre_averaged(:),group_ab(n,averaged_pre_chl_iter(:)),'.c');
% % % %     ylabel('Abundance (cell.cm-3)')
% % % %     xlabel('Chl (µg/L)')
% % % %     title([file(1,n+1) mdl.Rsquared.Ordinary]);
% % % %     subtitle('focus area, precise replacement');
%  end
% 
% % %% SST
% index_focus=[];
% for i = 1:length(dataCYTONEW) %detect index of values in the focus area
%     if (38 < dataCYTONEW(i,end) && 39 > dataCYTONEW(i,end) && 2.5 < dataCYTONEW(i,end-1) && dataCYTONEW(i,end-1) <3.5)
%         index_focus(length(index_focus)+1)=i;  
%     end
% end
% % lat_focus = dataCYTONEW(index_focus,end);
% % lon_focus = dataCYTONEW(index_focus,end-1);
% % lat_sst=[];
% % lon_sst=[];
% % sst_ab=[];
% % delta_lat_sst=[];
% % delta_lon_sst=[];
% % sst_from_averaged=[];
% % lat_sst_from_averaged = [];
% % lon_sst_from_averaged = [];
% % averaged_sst_iter = [];
% % 
% % sst_from_pre_averaged=[];
% % lat_sst_from_pre_averaged = [];
% % lon_sst_from_pre_averaged = [];
% % averaged_pre_sst_iter = [];
% % figure;
% % hold on
% % for n = 1:11
% %     averaged_sst_iter = [];
% %     sst_iter = [];
% % for i = 1:548
% % group_ab(n,i) = dataCYTONEW(i,6+n); 
% % date = dataCYTONEW(i,3);
% %     if date == 30 
% %         sst_data = sst_file(:,:,2);
% %         averaged_pre_sst = nanmean(sst_file(:,:,1:3),3);
% %     else
% %         sst_data = sst_file(:,:,date+2);
% %     end   
% %     
% %         [d, iy ] = min( abs( D.lat-lat_abundance(i) ) );        
% %         lat_sst(i)=iy;%indice des latitudes
% %         delta_lat_sst(i) = D.lat(lat_sst(i))-lat_abundance(i);
% %         [d, ix ] = min( abs( D.lon-lon_abundance(i) ) );
% %         lon_sst(i)=ix;%indice des longitudes
% %         delta_lon_sst(i) = D.lon(lon_sst(i))-lon_abundance(i);
% %         %sst_ab(i)=sst_data(lat_sst(i),lon_sst(i));
% %         if isnan(sst_data(lat_sst(i),lon_sst(i)))
% %            % standard mean 
% % %             averaged_sst_iter(length(averaged_sst_iter)+1) = i;
% % %             sst_ab(i) = averaged_sst(lat_sst(i),lon_sst(i));
% % %             sst_from_averaged(length(sst_from_averaged)+1) = averaged_sst(lat_sst(i),lon_sst(i));
% % %             lat_sst_from_averaged(length(lat_sst_from_averaged)+1) = lat_sst(i);
% % %             lon_sst_from_averaged(length(lon_sst_from_averaged)+1)= lon_sst(i);
% % 
% %             % more precise averaging
% %             averaged_pre_sst_iter(length(averaged_pre_sst_iter)+1) = i;
% %             %sst_ab(i) = averaged_pre_sst(lat_sst(i),lon_sst(i));
% %             sst_ab(i) = NaN;
% %             sst_from_pre_averaged(length(sst_from_pre_averaged)+1) = averaged_pre_sst(lat_sst(i),lon_sst(i));
% %             lat_sst_from_pre_averaged(length(lat_sst_from_pre_averaged)+1) = lat_sst(i);
% %             lon_sst_from_pre_averaged(length(lon_sst_from_pre_averaged)+1)= lon_sst(i);
% %         else
% %              sst_ab(i) = sst_data(lat_sst(i),lon_sst(i));
% %              sst_iter(length(sst_iter)+1) = i;
% %         end
% % end
% % % mdl = fitlm(sst_ab(:),group_ab(n,:));
% % % R2_matrix(n,1)=mdl.Rsquared.Ordinary;
% % %     subplot(3,4,n)
% % %     plot(mdl);
% % %     %hold on
% % %     %plot(sst_from_averaged(:),group_ab(n,averaged_sst_iter(:)),'.g');
% % %     ylabel('Abundance (cell.cm-3)')
% % %     xlabel('SST (°C)')
% % %     title([file(1,n+1) mdl.Rsquared.Ordinary]);
% % %     subtitle('focus area, no precise replacement');
% % end
% 
% 
% %% CHL relation with focus on the time evolution of the data
% chl_mat_ab_time = [];
% for n = 1:11
%     for i = 1:548
%         chl_mat_ab_time(i,1,n) = datenum(dataCYTONEW(i,1),dataCYTONEW(i,2), dataCYTONEW(i,3), dataCYTONEW(i,4), dataCYTONEW(i,5), dataCYTONEW(i,6));   %time
%         chl_mat_ab_time(i,2,n) = dataCYTONEW(i,6+n); %abundance
%         chl_mat_ab_time(i,3,n) = chl_ab(i); %chl
%         chl_mat_ab_time(i,4,n) = dataCYTONEW(i,end); %lat
%         chl_mat_ab_time(i,5,n) = dataCYTONEW(i,end-1); %lon
%     end
% end
% 
% % 
% figure;
% m_proj('mercator','lon',[min(lon) max(lon)],'lat',[min(lat)-0.1 max(lat)+0.2]);
% m_gshhs_h('save','gumby');
% m_usercoast('gumby','patch','w'); 
% m_grid('box','fancy','linestyle','-','gridcolor','k','backcolor','none');
% hold on
% m_plot(chl_mat_ab_time(:,5,1), chl_mat_ab_time(:,4,1),'k');
% hold on
% m_scatter(chl_mat_ab_time(:,5,1), chl_mat_ab_time(:,4,1), 50, chl_mat_ab_time(:,1,1), 'filled');
% m_northarrow(max(lon)-0.5,min(lat)+0.4,0.5,'type',2);
% ax = gca;
% ax.FontSize = 17; 
% xlabel('Longitude','FontSize',17);
% ylabel('Latitude','Fontsize',17);
% hcb=colorbar;
% hcb.Title.String = 'Date des mesures';
% caxis([chl_mat_ab_time(1,1,1) chl_mat_ab_time(end,1,1)]);
% cbdate('dd-mm-yyyy');
% 
% 
% % figure;
% % hold on
% % 
% % for i = 1:11
% %     mdl = fitlm(chl_mat_ab_time(:,3,i),chl_mat_ab_time(:,2,i)); %replace ':' by 'index_focus' for only focus area
% %     subplot(3,4,i)
% %     plot(mdl)
% %     hold on
% %     scatter(chl_mat_ab_time(:,3,i),chl_mat_ab_time(:,2,i),15,chl_mat_ab_time(:,1,i),'filled');
% %     ylabel('Abundance (cell.cm-3)')
% %     xlabel('Chl (µg/L)')
% %     colorbar;
% %     cbdate('yy-mm-dd');
% %     title([file(1,i+1) mdl.Rsquared.Ordinary])
% %     subtitle('total area, no replacement')
% % end   
% % 
% % %% SST relation with focus on the time evolution of the data
% % sst_mat_ab_time = [];
% % for n = 1:11
% %     for i = 1:548
% %         sst_mat_ab_time(i,1,n) = datenum(dataCYTONEW(i,1),dataCYTONEW(i,2), dataCYTONEW(i,3), dataCYTONEW(i,4), dataCYTONEW(i,5), dataCYTONEW(i,6));   %time
% %         sst_mat_ab_time(i,2,n) = datenum(dataCYTONEW(i,6+n)); %abundance
% %         sst_mat_ab_time(i,3,n) = sst_ab(i); %sst
% %         sst_mat_ab_time(i,4,n) = dataCYTONEW(i,end); %lat
% %         sst_mat_ab_time(i,5,n) = dataCYTONEW(i,end-1); %lon
% %     end
% % end
% % 
% % % map of measurements
% % % figure;
% % % m_proj('mercator','lon',[min(lon_focus) max(lon_focus)],'lat',[min(lat_focus) max(lat_focus)]);
% % % m_gshhs_h('save','gumby');
% % % m_usercoast('gumby','patch','w'); 
% % % m_grid('box','fancy','linestyle','-','gridcolor','k','backcolor','none');
% % % hold on
% % % m_scatter(sst_mat_ab_time(index_focus,5,1), sst_mat_ab_time(index_focus,4,1), 45, sst_mat_ab_time(index_focus,1,1), 'filled');
% % % colorbar;
% % % cbdate('yy-mm-dd')
% % % title('Focus area');
% % 
% % figure;
% % hold on
% % for i = 1:11
% %     mdl = fitlm(sst_ab(:),group_ab(i,:));
% %     subplot(3,4,i)
% %     plot(mdl)
% %     hold on
% %     scatter(sst_mat_ab_time(:,3,i),sst_mat_ab_time(:,2,i),15,sst_mat_ab_time(:,1,i),'filled'); %replace ":" by index_focus pour focus area
% %     ylabel('Abundance (cell.cm-3)')
% %     xlabel('SST (°C)')
% %     colorbar;
% %     cbdate('yy-mm-dd');
% %     title([file(1,i+1) mdl.Rsquared.Ordinary])
% %     subtitle('total area, no replacement')
% % end   
% % %%
% % mdl = fitlm(chl_ab(index_focus),sst_ab(index_focus));
% %     
% % figure;
% % hold on
% %     plot(mdl)
% %     hold on
% %     scatter(chl_ab(index_focus),sst_ab(index_focus),40,sst_mat_ab_time(index_focus,1,i),'filled'); %replace ":" by index_focus pour focus area
% %     xlabel('CHL (µg/L)')
% %     ylabel('SST (°C)')
% %     colorbar;
% %     cbdate('yy-mm-dd');
% %     title(['CHL/SST, focus area, no replacement'])
% %     subtitle(mdl.Rsquared.Ordinary)  
% %     
% %  
% % %% SST in function of T(tsg)
% % lat_sst=[];
% % lon_sst=[];
% % sst_ab=[];
% % delta_lat_sst=[];
% % delta_lon_sst=[];
% % sst_from_averaged=[];
% % lat_sst_from_averaged = [];
% % lon_sst_from_averaged = [];
% % averaged_sst_iter = [];
% % % figure;
% % % hold on
% % for n = 1:11
% % for i = 1:length(dataCYTONEW)
% % group_ab(n,i) = dataCYTONEW(i,6+n); 
% % date = dataCYTONEW(i,3);
% %     if date == 30 
% %         sst_data = sst_file(:,:,2);
% %     else
% %         sst_data = sst_file(:,:,date+2);
% %     end   
% %     
% %         [d, iy ] = min( abs( D.lat-dataCYTONEW(i,end)) );         
% %         lat_sst(i)=iy;%indice des latitudes
% %         %delta_lat_sst(i) = D.lat(lat_chl(i))-lat_abundance(i);
% %         [d, ix ] = min( abs( D.lon-dataCYTONEW(i,end-1)));
% %         lon_sst(i)=ix;%indice des longitudes
% %         %delta_lon_sst(i) = D.lon(lon_sst(i))-lon_abundance(i);
% %         sst_ab(i)=sst_data(lat_sst(i),lon_sst(i));
% %         if isnan(sst_data(lat_sst(i),lon_sst(i)))
% %             %veraged_sst_iter(length(averaged_sst_iter)+1) = i;
% %             %sst_ab(i) = averaged_sst(lat_sst(i),lon_sst(i));
% %             sst_ab(i) = NaN;
% %             sst_from_averaged(length(sst_from_averaged)+1) = averaged_sst(lat_sst(i),lon_sst(i));
% %             lat_sst_from_averaged(length(lat_sst_from_averaged)+1) = lat_sst(i);
% %             lon_sst_from_averaged(length(lon_sst_from_averaged)+1)= lon_sst(i);
% %         else
% %             sst_ab(i) = sst_data(lat_sst(i),lon_sst(i));
% %         end
% % end
% % end
%%
% mdl = fitlm(sst_ab(:),T(:));
% R2_matrix(n,1)=mdl.Rsquared.Ordinary;
% figure;
%     plotAdded(mdl);
%     hold on
%     %scatter(sst_ab(:),T(start:stop),30,datenum(dataCYTONEW(start:stop,1),dataCYTONEW(start:stop,2), dataCYTONEW(start:stop,3), dataCYTONEW(start:stop,4), dataCYTONEW(start:stop,5), dataCYTONEW(start:stop,6)),'filled');
%     grid on
% grid minor
%     %scatter(sst_ab(:),T(:),40,sst_mat_ab_time(:,1,1),'filled'); %replace ":" by index_focus pour focus area
%     ylabel('T tsg (°C)')
%     xlabel('SST (°C)')
%     %colorbar;
%     %cbdate('yy-mm-dd');
%     title(['SST vs. T tsg']);
%     subtitle(mdl.Rsquared.Ordinary)
% % 
disp('Abundance_relations : done')