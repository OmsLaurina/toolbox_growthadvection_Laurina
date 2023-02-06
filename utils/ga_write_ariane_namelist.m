function ga_write_ariane_namelist(name_curr,varargin)


%% GA_WRITE_ARIANE_NAMELIST: write the namelist file used by Ariane
%
% ga_write_ariane_namelist(name_curr,varargin)
%
% 'dt' 				time step for calculation, in days
% 'nbdays_advec' 	number of days to run the advection
% 'dtcurr'			current time step, in days
% 'itime_curr'		[1 curr.nb_time] par défaut, ITIME FOR CURRENTS (can be different from nbdays_advec if currents not daily)
% 'backwards'
%
% Monique Messié, May 2021


global dir_ariane_global

[arg,flag]=ga_read_varargin(varargin,{'dt',[],'nbdays_advec',[],'dtcurr',[],'itime_curr',[]},{'backwards'});
if isempty(arg.dt), error('Need to provide dt'), end
if isempty(arg.nbdays_advec), error('Need to provide nbdays_advec'), end
if isempty(arg.dtcurr), error('Need to provide dtcurr'), end
if isempty(arg.itime_curr), error('Need to provide itime_curr'), end
if flag.backwards, forback='backward'; else, forback='forward'; end


% getting grid size
infocurr=ncinfo([dir_ariane_global,'currents_data/meshmask_',name_curr,'.nc'],'xt');
gridsize=infocurr.Size;



% namelist file
filename=[dir_ariane_global,'namelist'];
fid=fopen(filename,'w');

fprintf(fid,'&ARIANE\n');

fprintf(fid,'	key_alltracers =.FALSE.,\n');
fprintf(fid,'	key_sequential =.TRUE.,\n');
fprintf(fid,'	key_ascii_outputs =.TRUE.,\n');
fprintf(fid,'	mode =''qualitative'',\n');
fprintf(fid,['	forback =''',forback,''',\n']);
fprintf(fid,'	bin =''nobin'',\n');
fprintf(fid,'	init_final =''init'',\n');
fprintf(fid,'	nmax =300000,\n');						% max number of particles
fprintf(fid,'	tunit =86400.,\n');						% time unit = daily (works with matlab time)
fprintf(fid,['	ntfic =',num2str(arg.dtcurr),',\n']);	% currents time resolution
fprintf(fid,'	tcyc =86400.,\n');
fprintf(fid,'/\n');

fprintf(fid,'&OPAPARAM\n');
fprintf(fid,['	imt =',num2str(gridsize(1)),'\n']);
fprintf(fid,['	jmt =',num2str(gridsize(2)),',\n']);
fprintf(fid,'	kmt =2,\n');
fprintf(fid,['	lmt =',num2str(arg.itime_curr(2)-arg.itime_curr(1)+1),',\n']);
fprintf(fid,'	key_periodic =.FALSE.,\n');
fprintf(fid,'	key_jfold =.FALSE.,\n');
fprintf(fid,'	key_computew =.TRUE.,\n');
fprintf(fid,'	key_partialsteps =.FALSE.,\n');
fprintf(fid,'/\n');
fprintf(fid,'&SEQUENTIAL\n');
fprintf(fid,'	maxcycles =1,\n');
fprintf(fid,'/\n');
fprintf(fid,'&QUALITATIVE\n');
fprintf(fid,['	delta_t =',num2str(round(3600*24*arg.dt)),'.,\n']);
fprintf(fid,'	frequency =1,\n');	
fprintf(fid,['	nb_output =',num2str(arg.nbdays_advec/arg.dt),',\n']);
fprintf(fid,'	key_region =.FALSE.,\n');
fprintf(fid,'/\n');

fprintf(fid,'&ZONALCRT\n');
fprintf(fid,['	c_dir_zo =''./currents_data'',\n']);
fprintf(fid,['	c_prefix_zo =''',name_curr,'_'',\n']);
fprintf(fid,['	ind0_zo =',num2str(arg.itime_curr(1)),',\n']);
fprintf(fid,['	indn_zo =',num2str(arg.itime_curr(2)),',\n']);
fprintf(fid,'	maxsize_zo =7,\n');								% time stamps in currents name (here 7 b/c 0000001)
fprintf(fid,'	c_suffix_zo =''.nc'',\n');
fprintf(fid,'	nc_var_zo =''u'',\n');
fprintf(fid,'	nc_var_eivu =''NONE'',\n');
fprintf(fid,'	nc_att_mask_zo =''NONE'',\n');
fprintf(fid,'/\n');

fprintf(fid,'&MERIDCRT\n');
fprintf(fid,['	c_dir_me =''./currents_data'',\n']);
fprintf(fid,['	c_prefix_me =''',name_curr,'_'',\n']);
fprintf(fid,['	ind0_me =',num2str(arg.itime_curr(1)),',\n']);
fprintf(fid,['	indn_me =',num2str(arg.itime_curr(2)),',\n']);
fprintf(fid,'	maxsize_me =7,\n');
fprintf(fid,'	c_suffix_me =''.nc'',\n');
fprintf(fid,'	nc_var_me =''v'',\n');
fprintf(fid,'	nc_var_eivv =''NONE'',\n');
fprintf(fid,'	nc_att_mask_me =''NONE'',\n');
fprintf(fid,'/\n');


fprintf(fid,'&MESH\n');
fprintf(fid,['	dir_mesh =''./currents_data'',\n']);
fprintf(fid,['	fn_mesh =''meshmask_',name_curr,'.nc'',\n']);
fprintf(fid,'	nc_var_xx_tt =''xt'',\n');
fprintf(fid,'	nc_var_xx_uu =''xu'',\n');
fprintf(fid,'	nc_var_yy_tt =''yt'',\n');
fprintf(fid,'	nc_var_yy_vv =''yv'',\n');
fprintf(fid,'	nc_var_zz_ww =''zw'',\n');
fprintf(fid,'	nc_var_e2u =''e2u'',\n');
fprintf(fid,'	nc_var_e1v =''e1v'',\n');
fprintf(fid,'	nc_var_e1t =''e1t'',\n');
fprintf(fid,'	nc_var_e2t =''e2t'',\n');
fprintf(fid,'	nc_var_e3t =''e3t'',\n');
fprintf(fid,'	nc_var_tmask =''tmask'',\n');
fprintf(fid,'	nc_mask_val =0.,\n');
fprintf(fid,'/\n');

fclose(fid);



return
