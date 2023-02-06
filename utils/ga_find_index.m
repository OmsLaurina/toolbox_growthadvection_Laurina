function index=ga_find_index(series,target_value)


%% GA_FIND_INDEX: returns the positions of a value in a series, as decimal points such that interp1(1:length(series),index)=target_value.
% The function iteratively goes through all data points and, when crossing the target value, determines its position.
%
% index=ga_find_index(series,target_value)
%
% Monique Messié, 2021 for public version


index=[];
if isnan(target_value), return, end

ipts=find(~isnan(series),1,'first');
while ipts<length(series)
	if series(ipts)==target_value, index=[index,ipts]; ipts=ipts+1; 
	else
		if series(ipts)>target_value, series=-series; target_value=-target_value; end
		while (series(ipts)<target_value || isnan(series(ipts))) && ipts<length(series), ipts=ipts+1; end
		if min(~isnan(series(ipts-1:ipts)))
			coef=interp1(series(ipts-1:ipts),0:1,target_value); 
			index=[index,coef+ipts-1]; 
		end
	end
end
index=unique(index);



return
