function list_filenames=ga_dir2filenames(rep_data,pattern)


%% GA_DIR2FILENAMES: gets the list of filenames in a directory, following a pattern
% 
% list_filenames=ga_dir2filenames(rep_data,pattern);
%
% Monique Messié, mai 2021

if nargin<2, pattern='*'; end

if ~isempty(rep_data) && rep_data(end)~='/'; rep_data=[rep_data,'/']; end
if ~contains(pattern,'*'), list_filenames={pattern}; 
else, list_dir=dir([rep_data,pattern]);
	list_filenames={}; 
	for iname=1:size(list_dir,1), filename=list_dir(iname).name;
		if ~max(strcmp(filename,{'.','..'}))
			list_filenames=[list_filenames,{filename}]; 
		end
	end
end


return
