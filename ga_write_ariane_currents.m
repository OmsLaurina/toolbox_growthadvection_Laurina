function ga_write_ariane_currents(curr,name_curr)


%% GA_WRITE_ARIANE_CURRENTS: writes current netcdf files in Ariane format
%
% Use:
% ga_write_ariane_currents(curr,name_curr)
%
% Required inputs:
% 	curr 		structure containing: 	.curr (3D matrix of dimensions lat,lon,time; complex numbers with real=u, imag=v)
%										.lat, .lon, .time (lon, lat, time are 1D vectors)
% 	name_curr 	unique identifier used to name the current dataset, will be an argument of GA processing
%
% Monique Messié, May 2021


global dir_ariane_global

if nargin<2, error('Need to provide both curr and name_curr!!!'); end
if ~isa(curr,'struct'), error('curr needs to be a structure containing .curr, .lon, .lat, .time'), end
if size(curr.curr,1)~=length(curr.lat) || size(curr.curr,2)~=length(curr.lon) size(curr.curr,3)~=length(curr.time)
	error('.curr does not have the correct dimensions: should be lat, lon, time'), end
if isempty(name_curr), error('name_curr cannot be empty!!!'), end
if length(name_curr)>60, error('Ariane doesn''t work if name_curr is too long'), end


% ------------------------------------------------------------------------------------------------------------------------------------- %
%														WRITE CURRENTS NETCDF FILES														%
% ------------------------------------------------------------------------------------------------------------------------------------- %

% cleanup potential old files
unix(['cd ',dir_ariane_global,'currents_data/; \rm ',name_curr,'_0*; \rm meshmask_',name_curr,'.nc;']);

% Compute C-grid
grid_u=struct(); grid_v=struct(); grid_r=struct();
grid_u.lon=curr.lon+(curr.lon(2)-curr.lon(1))/2; grid_u.lat=curr.lat;
grid_v.lon=curr.lon; grid_v.lat=curr.lat+(curr.lat(2)-curr.lat(1))/2;
% Put everything on the same dimensions 
grid_u.lon=grid_u.lon(1:end); grid_u.lat=grid_u.lat(1:end);
grid_v.lon=grid_v.lon(1:end); grid_v.lat=grid_v.lat(1:end); 
grid_r.lon=curr.lon(1:end); grid_r.lat=curr.lat(1:end);
nb_lon=length(grid_r.lon); nb_lat=length(grid_r.lat);
% Compute 2D grids
[curr.lon2D,curr.lat2D]=meshgrid(curr.lon,curr.lat);
[grid_u.lon2D,grid_u.lat2D]=meshgrid(grid_u.lon,grid_u.lat);
[grid_v.lon2D,grid_v.lat2D]=meshgrid(grid_v.lon,grid_v.lat);
[grid_r.lon2D,grid_r.lat2D]=meshgrid(grid_r.lon,grid_r.lat);

% Regrid currents on C-grid
disp('Interpol ucurr on u-grid...')
curr_u=nan(length(grid_u.lat),length(grid_u.lon),length(curr.time));
for itime=1:length(curr.time)
	curr_u(:,:,itime)=interp2(curr.lon2D,curr.lat2D,real(curr.curr(:,:,itime)),grid_u.lon2D,grid_u.lat2D);
end
curr_u(isnan(curr_u))=0;
disp('Interpol vcurr on v-grid...')
curr_v=nan(length(grid_v.lat),length(grid_v.lon),length(curr.time));
for itime=1:length(curr.time)
	curr_v(:,:,itime)=interp2(curr.lon2D,curr.lat2D,imag(curr.curr(:,:,itime)),grid_v.lon2D,grid_v.lat2D);
end
curr_v(isnan(curr_v))=0;

% currents files
disp('Write current netcdf files...')
for itime=1:length(curr.time), itime_str=['000000',num2str(itime)]; itime_str=itime_str(end-6:end);

	ncfile=[dir_ariane_global,'currents_data/',name_curr,'_',itime_str,'.nc'];
	nccreate(ncfile,'u','Dimensions',{'imt',nb_lon,'jmt',nb_lat,'kmt',2,'lmt',1})
	nccreate(ncfile,'v','Dimensions',{'imt',nb_lon,'jmt',nb_lat,'kmt',2,'lmt',1})
	nccreate(ncfile,'time','Dimensions',{'lmt',1})

	ncwrite(ncfile,'u',permute(repmat(curr_u(:,:,itime),1,1,2),[2 1 3 4]))
	ncwrite(ncfile,'v',permute(repmat(curr_v(:,:,itime),1,1,2),[2 1 3 4]))
	ncwrite(ncfile,'time',curr.time(itime))
	ncwriteatt(ncfile,'/','creation_date',datestr(now));
	ncwriteatt(ncfile,'/','comment','Currents interpolated on C-grid for Ariane, time is matlab format');

end
disp('Done!')


% ------------------------------------------------------------------------------------------------------------------------------------- %
%													WRITE CURRENTS MESHGRID FILE														%
% ------------------------------------------------------------------------------------------------------------------------------------- %

% calcul des facteurs d'echelle
R=6380e3; fac=pi/180; 
dx=grid_r.lon(2)-grid_r.lon(1); dy=grid_r.lat(2)-grid_r.lat(1); dz=10;
e1t=dx*fac*R*cos(grid_r.lat2D*fac);
e1u=dx*fac*R*cos(grid_u.lat2D*fac);
e1v=dx*fac*R*cos(grid_v.lat2D*fac);
e1f=dx*fac*R*cos(grid_v.lat2D*fac);
e2t=zeros(nb_lat,nb_lon); e2t(:,:)=dy*fac*R; e2u=e2t; e2v=e2t; e2f=e2t;
e3t=zeros(2,1); e3t(:)=dz; e3w=e3t;
tmask=repmat((1:2)',1,nb_lat,nb_lon);

% mask file
ncfile=[dir_ariane_global,'currents_data/meshmask_',name_curr,'.nc'];
for varname={'xt','xu','xv','xf','yt','yu','yv','yf','zt','zw','e1t','e1u','e1v','e1f','e2t','e2u','e2v','e2f','e3t','tmask'}, varname=varname{:};
	switch varname
	case 'tmask', dim_variable={'imt',nb_lon,'jmt',nb_lat,'kmt',2};
	case {'zt','zx'}, dim_variable={'kmt',2};
	otherwise, dim_variable={'imt',nb_lon,'jmt',nb_lat};
	end
	nccreate(ncfile,varname,'Dimensions',dim_variable)
end

ncwrite(ncfile,'xt',permute(grid_r.lon2D,[2 1]));
ncwrite(ncfile,'xu',permute(grid_u.lon2D,[2 1]));
ncwrite(ncfile,'xv',permute(grid_v.lon2D,[2 1]));
ncwrite(ncfile,'xf',permute(grid_u.lon2D,[2 1]));
ncwrite(ncfile,'yt',permute(grid_r.lat2D,[2 1]));
ncwrite(ncfile,'yu',permute(grid_u.lat2D,[2 1]));
ncwrite(ncfile,'yv',permute(grid_v.lat2D,[2 1])); 
ncwrite(ncfile,'yf',permute(grid_v.lat2D,[2 1]));
ncwrite(ncfile,'zt',[0 dz]);
ncwrite(ncfile,'zw',[dz/2 1.5*dz]);
for varname={'e1t','e1u','e1v','e1f','e2t','e2u','e2v','e2f','e3t'}, varname=varname{:}; ncwrite(ncfile,varname,permute(eval(varname),[2 1])); end
ncwrite(ncfile,'tmask',permute(tmask,[3 2 1]));
ncwriteatt(ncfile,'/','creation_date',datestr(now));
ncwriteatt(ncfile,'/','comment','Meshgrid for currents interpolated on C-grid for Ariane');



return
