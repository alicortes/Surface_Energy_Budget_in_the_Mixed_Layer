function Ananstd = nanstd2(A)

Ananstd = std(A(~isnan(A)));

end