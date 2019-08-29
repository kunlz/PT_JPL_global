%==============%
% PT_JPL model %
%==============%
% This program is developed for calculating daily or (monthly) actual ET

% Code: Kun Zhang, Lanzhou University, China
% Questions to: zhangk12@lzu.edu.cn
% First version: 2015.3.11 
% Last Modified: 2016.10.9

% Reference: 
% Fisher et al.(2008), Global estimates of the land-atmosphere water flux based on...
% monthly AVHRR and ISLSCP-II data, validated at 16 FLUXNET sites, Remote Sens Environ, 112(3), 901-919...
% doi:10.1016/j.rse.2007.06.025.
%---------------------------------------------------------
function LE = PT_JPL(Rn, Ta, RH, G, NDVI, EVI)
%-----------
% Input:
% Rn -- Net Radiation (W m-2)
% Ta -- Air temperature (Celsius)
% RH -- Relative humidity (0.01)
% G  -- Soil heat flux (W m-2)
% NDVI--Normalized differential vegetation index
% EVI --Enhanced vegetation index
%-----------
% Constant Parameter:
    k_Rn = 0.6;       % Fisher,2008; Imens&Lemur,1969
    k_PAR = 0.5;      % Fisher,2008; Ross,1976
    m1 = 1.2*1.136;   % Fisher,2008; Gao etal,2000; Huete,1988,2006
    m2 = 1;           % Fisher,2008
    b1 = 1.2*(-0.04); % Fisher,2008; Gao etal,2000; Huete,1988,2006
    b2 = -0.05;       % Fisher,2008
    alfa = 1.26;      % Fisher,2008; Priestley&Taylor,1972
    beta = 1;         % Fisher,2008; (kPa)
    gamma = 0.066;    % psychrometric constant,(kPa C-1)

    % Intermediate variable
    es(:,1) = 0.6108.*exp(17.27.*Ta./(Ta+237.3)); % satureation vapur pressure, kPa
    delta = 4098.*es(:,1)./(Ta+237.3).^2; 
    ea(:,1) = RH.*es(:,1); % actual vapor pressure, kPa
    VPD(:,1) = es(:,1)-ea(:,1); % vapor pressure deficit. kPa

    % Cal f_APAR during the study time (Yearly or..)
    f_APARmax = max(m1.*EVI+b1);

    % Cal Topt
    TTT = Rn.*Ta.*EVI./VPD;
    TTT(VPD==0,:) = [];
    Topt = Ta(TTT == max(TTT),1); % Optmum growth temperature, Celsius
    %Topt=25;         % Garcia,2013; Yao,2013 

    % Main program
    f_APAR = m1.*EVI+b1;       % Fraction of PAR absorbed by green vegetation cover
    f_IPAR = m2.*NDVI+b2;      % Fraction of PAR intercepted by total vegetation cover
    fc = f_IPAR;               % Fraction total vegetation cover
    LAI = -log(1-fc)./k_PAR;   % Leaf area index
    Rns = Rn.*exp(-k_Rn.*LAI); % Net radiation to the soil
    Rnc = Rn-Rns;              % Net radiation to the canopy
    fg = f_APAR./f_IPAR;       % Green canopy fraction
    ft = exp(-((Ta-Topt)./Topt).^2); % Plant temperature constraint
    fm = f_APAR./f_APARmax;    % Plant moisture constraint
    fsm = RH.^((VPD(:,1))./beta); % Soil moisture constraint
    fwet = RH.^4;              % Relative surface wetness

    % function cal
    [LEc, LEs, LEi] = cal(alfa,delta,gamma,fwet,fg,ft,fm,fsm,Rnc,Rns,G);
    LE = LEc+LEs+LEi;
end

function [LEc, LEs, LEi] = cal(alfa,delta,gamma,fwet,fg,ft,fm,fsm,Rnc,Rns,G)
    PJ = alfa.*delta./(delta+gamma);
    LEc = (1-fwet).*fg.*ft.*fm.*PJ.*Rnc;
    LEs = (fwet+fsm.*(1-fwet)).*PJ.*(Rns-G);
    LEi = fwet.*PJ.*Rnc;
end