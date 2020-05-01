load Carioca_T_2013_corrected.mat

%I need to get the max and the minimum temperature value at each day????

decimals = time - floor(time);
use = find(decimals == 0.0);

for i = 1:(length(use)-1)
    day_s = use(i);
    day_e = use(i+1);
    day_T = T(day_s:day_e,1);
    day_Max = max(day_T)
    day_Min = min(day_T)
    deltaT(i) = day_Max - day_Min
    clear day_s day_e day_max day_min
end

de

figure
plot(deltaT)