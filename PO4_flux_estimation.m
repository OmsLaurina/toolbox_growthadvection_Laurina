%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% PO4 FLUX ESTIMATION (method from Pulido-Villena, 2021) %%%%%%%%%%%%%
%%%% N. KIENTZ, 2022 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
clear all;
close all;
%%%%% IMPORT DATA %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Salinity  (lon, lat, depth, time)
S_file = dir('S_T/*sal*.nc');
salinity = double(ncread(S_file.name,'so'));
depth_sal = double(ncread(S_file.name,'depth'));%m
depth_sal = depth_sal(1:125); %to match level of depth for S, T and Nut
time_sal = double(ncread(S_file.name,'time'));%minutes since 1900-01-01 00:00:00
lat_sal = double(ncread(S_file.name,'lat'));
lon_sal = double(ncread(S_file.name,'lon'));

%%% Potential temperature (lon, lat, depth, time)
T_file = dir('S_T/*tem*.nc');
temperature = double(ncread(T_file.name,'thetao'));%°C
depth_tem = double(ncread(T_file.name,'depth'));%m
depth_tem = depth_tem(1:125); 
time_tem = double(ncread(T_file.name,'time'));%minutes since 1900-01-01 00:00:00
lat_tem = double(ncread(T_file.name,'lat'));
lon_tem = double(ncread(T_file.name,'lon'));

%%% MLD depth (lon, lat, time)
MLD_file = dir('MLD/*mld*.nc');
mld = double(ncread(MLD_file.name,'mlotst'));%m
time_mld = double(ncread(MLD_file.name,'time'));%minutes since 1900-01-01 00:00:00
lat_mld = double(ncread(MLD_file.name,'lat'));
lon_mld = double(ncread(MLD_file.name,'lon'));

%%% PHOSPHATE CONCENTRATION (lon, lat, depth, time)
NUT_file = dir('PHOSPHATE/*nut*.nc');
phosphate = double(ncread(NUT_file.name,'po4'));%mmol.m⁻³
depth_nut = double(ncread(NUT_file.name,'depth'));%m
time_nut = double(ncread(NUT_file.name,'time'));%seconds since 1970-01-01 00:00:00
lat_nut = double(ncread(NUT_file.name,'latitude'));
lon_nut = double(ncread(NUT_file.name,'longitude'));


%%% FIXED PARAMETERS
epsilon = 8*10^(-9); %6-10.10⁻⁹ W.kg⁻¹ entre 20 et 100m pour la zone couverte par PEACETIME 2017 (Cuypers et al., 2012)
rho_0 = 1025; %1025kg.m⁻³ (Lynne Talley, 2000 (SIO 210))
g = 9.81;%9.81 m.s⁻²

%%% PRESSURE ESTIMATION (lon, lat, depth)
% Z = -7x10⁻²P + 2x10⁻³P² (Leroy & Parthiot, 1997); Z (m); P (MPa)
P = ones(size(temperature,1),size(temperature,2),length(depth_tem));
depth = depth_tem;
for id = 1:length(depth_tem)
P(:,:,id) = P(:,:,id).* (7*10^(-2) - ((-7*10^(-2))^2 - 4*(2*10^(-3))*-depth(id))^(-1/2))/(2*2*10^(-3))*100; %*100 to convert MPa to dBar
end

%%% DENSITY ESTIMATION (UNESCO algorithms)
density = ones(size(temperature,1),size(temperature,2),length(depth_tem),size(temperature,4));
for it = 1:length(time_tem)
S = salinity(:,:,1:125,it);
T = temperature(:,:,1:125,it);
[SVAN,SIGMA]=swstate(S,T,P);
density(:,:,:,it)=SIGMA+1000;%to have density instead of anomaly density
end

%%% PHOSPHATE CONCENTRATION ACROSS ISOPYCNALS
% i.e. linearly fitting phosphate concentrations versus density
%%
%R2_matrix = [];
pente = [];
%pvalue_pente = [];
F_PO4_C = ones(169,96,21);
for it = 1:length(time_tem)
    
for ix = 1:length(lon_tem)
    for iy = 1:length(lat_tem)
            mdl=fitlm(squeeze(density(ix,iy,:,it)),squeeze(phosphate(ix,iy,:,it)));
            %R2_matrix(ix,iy,it) = mdl.Rsquared.Ordinary;
            table = mdl.Coefficients;
            table = table2array(table);
            pente(ix,iy,it) = table(2,1);
            %pvalue_pente(ix,iy,it) = table(2,4);
    end
end

%%% PO4 FLUX
%%% F_PO4 = -0.2 x epsilon x (rho_0/g) x d_C/d_rho (Pulido-Villena, 2021)
F_PO4 = ones(169,96);

F_PO4 = F_PO4.*(-0.2 * epsilon * (rho_0/g) * pente(:,:,it)); % kg PO4.m⁻².s⁻¹
F_PO4 = F_PO4.* (86400 * 1000 * 1/94.971 * 10^6); %µmol PO4 .m⁻².j⁻¹
F_PO4_C(:,:,it) = F_PO4 .*106*10^(-3)./mld(:,:,it); %mmol C.m⁻3.j⁻1 ; divided by the mld depth
end

%% create struct for nut inputs for every position and every time
grid = 0.042; %grid resolution of CMEMS products
Nsupply_PROTEVS=struct();
Nsupply_PROTEVS.lat = lat_tem(1):grid:lat_tem(end);
Nsupply_PROTEVS.lon = lon_tem(1):grid:lon_tem(end);
Nsupply_PROTEVS.time = datenum(2018,04,30):datenum(2018,05,20);

F_PO4_C_rg = ones(length(Nsupply_PROTEVS.lon),length(Nsupply_PROTEVS.lat),length(Nsupply_PROTEVS.time));
for i = 1:size(F_PO4_C_rg,1)
    for j = 1:size(F_PO4_C_rg,2)
        F_PO4_C_rg(i,j,:)=F_PO4_C(floor(i)+1,floor(j)+1,:);
    end
end
Nsupply_PROTEVS.Nsupply = F_PO4_C_rg;
Nsupply_PROTEVS.Nsupply_ini = F_PO4_C_rg(:,:,1);

for it = 1:21
    F_PO4_it = F_PO4_C_rg(:,:,it)';
Nsupply_PROTEVS.lon_interp = min(lon_nut):0.02:max(lon_nut);
Nsupply_PROTEVS.lat_interp = min(lat_nut):0.02:max(lat_nut);
[Xq, Yq] = meshgrid(Nsupply_PROTEVS.lon_interp,Nsupply_PROTEVS.lat_interp);
Nsupply_PROTEVS.Nsupply_interp(:,:,it) = interp2(Nsupply_PROTEVS.lon,Nsupply_PROTEVS.lat,F_PO4_it,Xq, Yq);
end
save('inputs/Nsupply_PROTEVS','Nsupply_PROTEVS')
%%
% figure('DefaultAxesFontSize',22)
% m_proj('mercator','lon',[0 7],'lat',[36 40]);
% m_pcolor(lon_tem,lat_tem,F_PO4_C(:,:,1)')
% shading flat
% m_usercoast('gumby','patch','w');
% m_grid('box','fancy','linestyle','-','gridcolor','none','backcolor','none');
% hold on
% cbar = colorbar;
% colorTitleHandle = get(cbar,'Title');
% titleString = ({'Flux de PO_4 (mmolC.m⁻³.j⁻¹)',''});
% set(colorTitleHandle ,'String',titleString);
toc



