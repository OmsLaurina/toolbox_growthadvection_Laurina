function varargout=ga_read_varargin(varargin_input,default_parameters,flag_list)

%% GA_READ_VARARGIN: gets parameters given by a cell (e.g. varargin) and returns them as arguments or flags
% Varargin_input should contain both pairs of arguments (eg, 'arg_name',value) and flags as a single string.
% The function compares parameters and flags given in default_parameters and flag_list against varargin_inputs, and returns the corresponding values
% (true/false for flags, and custom/default parameters for arguments,  depending on whether they are found)
% The function also returns isfound (says if the parameter was found in varargin_input or if the default value was found, only for arguments) 
% and remaining_varargin which is the part of varargin_input that wasn't used.
%
% varargout=ga_read_varargin(varargin_input,default_parameters,flag_list);
%	All inputs are cells, with default_parameters containing pairs, flag_list a list of strings, and varargin_input a combination of both.
%
% Examples of use:
% [arg,flag,isfound,remaining_varargin]=ga_read_varargin(varargin,default_parameters,flag_list);
% [arg,isfound,remaining_varargin]=ga_read_varargin(varargin,default_parameters);
%
% Monique Messié, 2021 for public version


if nargin<3, flag_list={}; end
ikeep=true(size(varargin_input));

% Looking for arguments
arg=struct(); isfound=struct();
for ivar=1:2:length(default_parameters)-1, varname=default_parameters{ivar};
	index=find(strcmp(varargin_input,varname),1);
	if isempty(index)
		arg.(varname)=default_parameters{ivar+1}; isfound.(varname)=false; 
	elseif index==length(varargin_input), arg.(varname)=default_parameters{ivar+1}; isfound.(varname)=true; ikeep(index)=0; 
	else, arg.(varname)=varargin_input{index+1}; isfound.(varname)=true; ikeep(index:index+1)=0; 
	end
end

% Looking for flags
flag=struct();
for varname=flag_list, varname=varname{:}; 
	if max(strcmp(varargin_input,varname))
		flag.(varname)=true; ikeep(index)=0;
	else, flag.(varname)=false; 
	end
end

% Output
if nargin<3, varargout={arg,isfound,varargin_input(ikeep)}; else, varargout={arg,flag,isfound,varargin_input(ikeep)}; end
if nargout==0, varargout=varargout(1); else, varargout=varargout(1:nargout); end


return

