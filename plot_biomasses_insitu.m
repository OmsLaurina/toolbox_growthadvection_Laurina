%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 PLOT BIOMASSES IN SITU                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Use function : ga_*.m                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Author : Laurina Oms                                       %%%
%%% Creation : 25/01/2023                                      %%%     
%%% Region : Med sea, Cruise : PROTEVS 2018                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all

index_hipp = 407:511;
dataCYTONEW = dlmread('inputs/data_CYTO_NEW.txt');

name_group = 'pico';
parameters = phytoP_settings(name_group);
ab_cyto = parameters.biomasse_insitu.biomasse;

year = dataCYTONEW(index_hipp,1);
month = dataCYTONEW(index_hipp,2);
day = dataCYTONEW(index_hipp,3);
dates = datenum([year,month,day]);

figure, hold on
	plot(dates,ab_cyto,'LineWidth',2, 'Color','k')
    datetick('x','mmmyy')
% 	if min(dates)>datenum(1900,1,1), datetick('x','keeplimits'), xlabel('Time')
% 	else, xlabel('Time (days)')
%     end
	title('Pico in situ')
    
% datetick('x','mmmyy')