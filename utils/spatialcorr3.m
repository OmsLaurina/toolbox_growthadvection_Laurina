function [ C ] = spatialcorr3(varargin)
%%Gives the Spatial Correlation between two variables for the specific time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%spatialcorr3 is a part of TimeSaving ToolBox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Created by Ankur Kumar                        Wednesday; January 03, 2018
%                                              Version: 1.0
%
% National Institute of Techonology, Rourkela
% Rourkela, Odisa - 769008, India
% Department of Earth and Atmospheric Sciences
% Email: ankurk017@gmail.com
%        416AS2025@nitrkl.ac.in
%
%
% Function:
%           For the spatial correlation between two variables, spatialcorr3
%           helps us to get on the desired result. Both the input matrices
%           should be in three dimension. One of the dimensions of a matrix
%           should reflect the time. If you have taken the mean of all time
%           steps (suggestion is not to take the mean) or if the data is
%           available in monthly mean format or daily mean format, then use
%           spatialcorr2.
%
% Syntax:
%           C1=spatialcorr3(A,B);
%           C2=spatialcorr3(A,B,2);
%
%Inputs:
%           First and second input should be the matrices between which you
%           want to find correlation.
%           Third dimension should be of dimension which reflects the time.
%           If you keep this argument empty, then by default it set to 3.
%
%
% Example:
%           A=randi(25,10,20,365);
%           B=randi(55,10,20,365);
%           C1=spatialcorr3(A,B,1);
%           C2=spatialcorr3(A,B,2);
%           C3=spatialcorr3(A,B,3);
%
%
% Please send your suggestions to the email id: ankurk017@gmail.com or
%                                               416AS2025@nitrkl.ac.in
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%spatialcorr3 is a part of TimeSaving ToolBox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin>=4
    error('Too many inputs. Only two inputs required.')
end
a=varargin{1};
if length(size(a))~=3
    error('Input matrices should be of dimension 3. Error in first input matrix. One of the dimensions should reflect the time. If you have taken the mean of all timesteps (suggestion is not to take the mean) or if the data is available in monthly mean format or daily mean format, then use spatialcorr2.')
end
b=varargin{2};
if length(size(b))~=3
    error('Input matrices should be of dimension 3. Error in first input matrix. One of the dimensions should reflect the time. If you have taken the mean of all timesteps (suggestion is not to take the mean) or if the data is available in monthly mean format or daily mean format, then use spatialcorr2.')
end
if nargin==3
    dim=varargin{3};
    if size(dim,1)==1 && size(dim,2)==1 && dim<=3
        dim=dim;
    else
        error('Third input should be singleton (Dimension must be 1*1) and must be less than 3. It should be the dimension of time series, along which you want to find correlation.')
    end
    if dim==2
        a=permute(a,[1,3,2]);
        b=permute(b,[1,3,2]);
    elseif dim==1
        a=permute(a,[2,3,1]);
        b=permute(b,[2,3,1]);
    elseif dim==3
        a=permute(a,[1,2,3]);
        b=permute(b,[1,2,3]);
    else
        error('Error in program')
    end
end
for i=1:size(a,1)
    for j=1:size(a,2)
        q=a(i,j,:);
        w=b(i,j,:);
        id=union(find(isnan(q)),find(isnan(w)));
        q(id)=[]; w(id)=[];
        C(i,j)=corr2(q,w);
    end
end
end

