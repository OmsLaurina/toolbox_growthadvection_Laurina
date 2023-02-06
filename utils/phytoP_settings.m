%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             FUNCTION - parameter for phytoP group          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs : name phytoP group                                 %%%
%%% Parameters from                                            %%%
%https://people.mio.osupytheas.fr/~doglioli/rapport_M2_Kientz.pdf%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Author : Laurina Oms                                       %%%
%%% Creation : 14/11/2022                                      %%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function param = phytoP_settings(name_group)
% load('inputs/mat_micro');
% load('inputs/mat_pico');
% load('inputs/iok');
% load('inputs/time2D');
% load('inputs/index_hipp');
dataCYTONEW = dlmread('inputs/data_CYTO_NEW.txt');
index_hipp = 407:511;

if strcmp(name_group, 'pico') 
    param.pos_group = 10;
    param.caxis_lim = [0.2 0.3];
    param.size_magnitude = 0.9;
%     mat_pico = mat_pico(iok);
%     abundance_matrix_pico_extp = reshape(mat_pico,[size(mat_pico,1)*size(mat_pico,2),1]);
%     for it = 1:length(time2D(1,:))-1
%         abundance_matrix_pico_extp = [abundance_matrix_pico_extp reshape(mat_pico,[size(mat_pico,1)*size(mat_pico,2),1])];
%     end
%     param.biomasse_extp = ab2biomasse((abundance_matrix_pico_extp(:,end-15:end)*1E6), param.size_magnitude); %convert cm3 to m3 (*1E6)
    param.biomasse_insitu = ab2biomasse((dataCYTONEW(index_hipp,param.pos_group)*1E6), param.size_magnitude);
    
elseif strcmp(name_group, 'micro')
    param.pos_group = 14;
    param.caxis_lim = [0.15 0.27];
    param.size_magnitude = 90;
%     mat_micro = mat_micro(iok);
%     abundance_matrix_micro_extp = reshape(mat_micro,[size(mat_micro,1)*size(mat_micro,2),1]);
%     for it = 1:length(time2D(1,:))-1
%         abundance_matrix_micro_extp = [abundance_matrix_micro_extp reshape(mat_micro,[size(mat_micro,1)*size(mat_micro,2),1])];
%     end
%     param.biomasse_extp = ab2biomasse((abundance_matrix_micro_extp(:,end-15:end)*1E6), param.size_magnitude);
    param.biomasse_insitu = ab2biomasse((dataCYTONEW(index_hipp,param.pos_group)*1E6), param.size_magnitude);
end
       
end
