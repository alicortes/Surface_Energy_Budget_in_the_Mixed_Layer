function Ananmean = nanmean2(A)

Ananmean = mean(A(~isnan(A)));

end