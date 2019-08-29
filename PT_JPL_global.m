% ------------ %
% PT_JPL model %
% ------------ %

% Code: Kun Zhang, Lanzhou University, China
% Questions to : zhangkun322@foxmail.com
% First version: 2015.3.11
% Last Modified: 2016.10.9

% Reference:
% Fisher et al.(2008), Global estimates of the land-atmosphere water flux based on...
% monthly AVHRR and ISLSCP-II data, validated at 16 FLUXNET sites, Remote Sens Environ, 112(3), 901-919...
% doi:10.1016/j.rse.2007.06.025.

%----------------------------------------------
% Forcing data:
% Rn --------- Net Radiation, w/m^2
% Ta --------- Air Temperature, C
% RH --------- Relative Humidity, 1
% G  --------- Soil Heat Flux, w/m^2
% NDVI ------- Normalized Difference Vegetation Index
% EVI --------
% f_APARmax -- the max f_APAR during the study period
%----------------------------------------------

function [E,Et,Es,Ei] = PT_JPL_global(Rn,Ta,RH,G,NDVI,EVI,f_APARA_max)
% -------------------------------------
% Parameters
beta = 1; % 1.0 kPa
m1   = 1.2*1.136; % 1.2*1.136 Gao et al,2000; Huete,1988,2006
b1   = 1.2*(-0.04); % Gao et al,2000; Huete,1988,2006
m2   = 1; % 1.0 Fisher,2008
b2   = -0.05; % -0.05 Fisher,2008
k_Rn = 0.6; % 0.6 Imens&Lemur,1969
k_PAR= 0.5; % 0.5 Ross,1976
Topt = 25; % 25 Fisher,2008 ???
alfa = 1.26; % 1.26 Priestley&Taylor,1972
gamma= 0.066; % psychrometric constant
% -------------------------------------
% Main Programe

es    = 0.6108 .* exp(17.27 .* Ta ./ (Ta + 237.3));
delta = 4098 .* es ./ (Ta + 237.3) .^ 2;
ea    = RH .* es;
VPD   = es -ea;

% Fraction of PAR absorbed by green vegetation cover
f_APAR = m1 .* EVI + b1;

% Fraction of PAR intercepted by total vegetation cover
f_IPAR = m2 .* NDVI + b2;

f_APAR(f_APAR<0) = 0;
f_IPAR(f_IPAR<0) = 0;

% Fraction total vegetation cover
fc          = f_IPAR;
fc(fc >= 1) = 0.95;

% Leaf area index
LAI = -log(1 - fc) ./ k_PAR;

% Net radiation to the soil
Rns = Rn .* exp(-k_Rn .* LAI);

% Net radiation to the canopy
Rnc = Rn - Rns;

% Green canopy fraction
fg  = f_APAR ./ f_IPAR;
fg(f_APAR>f_IPAR) = 1;

% Plant temperature constraint
ft  = exp(-((Ta - Topt) ./ Topt) .^ 2);

% Plant moisture constraint
fm  = f_APAR ./ f_APARA_max;

% Soil moisture constraint
fsm = RH .^ (VPD ./ beta);

% Relative surface wetness
fwet   = RH .^ 10;

% cal
[LEc, LEs, LEi] = cal(alfa,delta,gamma,fwet,fg,ft,fm,fsm,Rnc,Rns,G);
LE = LEc + LEs + LEi;

% convert to mm
E  = w2mm(LE, Ta, 24);
Et = w2mm(LEc, Ta, 24);
Es = w2mm(LEs, Ta, 24);
Ei = w2mm(LEi, Ta, 24);

end

function [LEc, LEs, LEi] = cal(alfa,delta,gamma,fwet,ft,fg,fm,fsm,Rnc,Rns,G)
PJ  = alfa.*delta./(delta+gamma);
LEc = (1-fwet).*fg.*ft.*fm.*PJ.*Rnc;
LEs = (fwet+fsm.*(1-fwet)).*PJ.*(Rns-G);
LEi = fwet.*PJ.*Rnc;
end

function [out] = w2mm(inn, Ta, conts)
lambda = 1e6.*(2.501-0.002361.*Ta);
Aa     = inn./lambda; % kg m-2 s-1
out    = Aa .*3600 .*conts; % kg m-2 day-1 ==> mm day-1
end