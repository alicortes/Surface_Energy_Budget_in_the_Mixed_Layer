[num txt raw] = xlsread('Temp_pt5m_HOBO.xlsx',1);
date = (raw(2:end,1));date2 = dnum2doy(datenum(date,'dd/mm/yyyy'));
dt = [15/24/60];ref=[0:dt:(1-dt)];m = length(ref);
date3 = unique(date2); n=length(date3);
hours = repmat(ref,1,n);time1 = date2 + hours';
T1 = num(:,2); %diffT1 = diff(T1);
clear num txt raw date date2 date3 dt hours m n ref
load Carioca_T_2013_corrected.mat
figure
plot(time,T1,'r');hold on; plot(time,T); ylabel('Temp (^oC)')
ylim([18 33])
legend('1.0m','3.0m','5.0m','7.0m','9.0m')
title('Corr data')
hold on; 