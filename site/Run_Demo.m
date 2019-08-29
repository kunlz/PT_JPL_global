clc;clear
load data

data = Data{1,1};
data = data(:,7:14);
Rn = data(:,1);
Ta = data(:,2);
RH = data(:,3)/100;
G = data(:,4);
obs = data(:,6);
NDVI = data(:,7);
EVI = data(:,8);

result1 = PT_JPL(Rn, Ta, RH, G, NDVI, EVI);
result2 = PT_JPL_v1(Rn, Ta, RH, G, NDVI, EVI);


plot([obs,result1,result2])
xlabel('Days')
ylabel('Latent heat flux (W/m^2)')
legend('obs','sim')
box on
grid on

R1 = corr(obs,result1)
R2 = corr(obs,result2)