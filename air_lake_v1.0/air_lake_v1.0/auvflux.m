%% Heat flux calculations

i=1;

testT = [];
testD = [];
testt = [];

while (i<=13)
    S = find(MDS == 0);
    testT = [testT, MTS(S,:)];
    testD = [testD, MTD(S,:)];
    testt = [testt, MtS(S,:)];
    i+1;
end

    