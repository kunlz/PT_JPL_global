
function [f_APARA_max] = get_fapara_max(EVI_3D)
% ---------------
% function input:
% EVI_3D -- Enhanced Vegetation Index, range in [0~1]
%           all EVI data for each year listed in a 3D matrix
%           z-axis is time series
% ---------------

% parameters:
m1 = 1.2*1.136;
b1 = 1.2*(-0.04);

EVImax = max(EVI_3D,[],3);
f_APARA_max = m1 .* EVImax + b1;
f_APARA_max(m1 == 0) = 0;
end