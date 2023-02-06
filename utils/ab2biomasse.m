%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    FUNCTION - convert abundance in biomasse (mmol C/m3)    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs : cytometry abundance, size magnitude of groups     %%%
%%% Parameters from                                            %%%
%https://people.mio.osupytheas.fr/~doglioli/rapport_M2_Kientz.pdf%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Author : Laurina Oms                                       %%%
%%% Creation : 14/11/2022                                      %%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bioV_QC_biom = ab2biomasse(abundances, size_magnitude)

%%
%coef beta0 et beta1 de la power law function = obtenu à partir d'une regression log log entre
%le FWS et la taille des billes SPECIFIQUE AU CYTOMETRE UTILISÉ 
beta0 = -5.8702;
beta1 = 0.9228;

a = 0.26;
b = 0.86;
%%

bioV_QC_biom = struct();

%%
% size_magnitude = 0.9 (Pico); 90 (Micro) [um]
bioV_QC_biom.biovolume = exp(beta0)*(size_magnitude)^beta1; %bioV formule "toute faite" : power law function (Olson,2003) [unité?]
bioV_QC_biom.biovolume = (4/3)*pi*(size_magnitude)^3; %bioV formule d'une sphère [um^3]
%%
bioV_QC_biom.Q_C = (((a*bioV_QC_biom.biovolume)*b))*1E-12/12.106;
bioV_QC_biom.biomasse = abundances.*bioV_QC_biom.Q_C; %mmol C/m3

return
