%% figure(1)
figure('name','BasicMet')
set(gcf,'position',goodfigsize(gcf))
ax = subplot(5,1,1);
h = plot(doy,SWin,'k');
ylabel('SW [W/m^2]')
set(ax,'ylim',[-25 1400]) %chad 1500
dynamicDateTicks([], [], 'mm/dd');

ax = subplot(5,1,2);
h = plot(doy,AirT,'k',doy,SrfT,'k:');
ylabel('T [\circC]')
legend('Air','Water @ 0.1 m','Location','BestOutside')
dynamicDateTicks([], [], 'mm/dd');

ax = subplot(5,1,3);
h = plot(doy,RH,'k');
ylabel('RH [%]')
set(ax,'ylim',[0 100])
dynamicDateTicks([], [], 'mm/dd');

ax = subplot(5,1,4);
h = plot(doy,WS,'k');
ylabel('U [m/s]')
xlabel('Day')
set(ax,'ylim',[0 max(WS)])
dynamicDateTicks([], [], 'mm/dd');

dum = find(WDir < 90);
WDir2 = WDir;
WDir2(dum) = WDir(dum)+360;
ax = subplot(5,1,5);
h = plot(doy,WDir2,'k');
ylabel('WDir [^o]')
xlabel('Day')
set(ax,'ylim',[90 450])
dynamicDateTicks([], [], 'mm/dd');

%% figure(2)
figure('name','HeatFluxes')
set(gcf,'position',goodfigsize(gcf))
ax = subplot(2,1,1);
h = plot(doy,[SE LE LWNet],doy,zeros(size(doy)),'k');
ylabel('Heat Flux [W/m^2]')
set(h,{'linestyle'},{'-'; '--'; '-.'; ':'})
legend('SE','LE','LWNet','Location','BestOutside')
set(ax,'ylim',[min([-300,min([SE LE LWNet])]) 100]) %chad -400
dynamicDateTicks([], [], 'mm/dd');

ax = subplot(2,1,2);
h = plot(doy,SRFFLX,doy,HETFLX,doy,HETFLX_tot,doy,zeros(size(doy)),'k:');
ylabel('Heat Flux [W/m^2]')
xlabel('Day')
legend('Surface','Effective','Location','BestOutside')
set(ax,'ylim',[min(min([SRFFLX HETFLX]))-50 max(max([SRFFLX HETFLX]))+50]);
dynamicDateTicks([], [], 'mm/dd');


%% figure(3)
figure('name','Atmosphere')
set(gcf,'position',goodfigsize(gcf))
ax = subplot(3,1,1);
h = plot(doy,L,'k');
ylabel('L_a')
set(ax,'ylim',[-20 20])
dynamicDateTicks([], [], 'mm/dd');

ax = subplot(3,1,2);
h = plot(doy,Sinfo.sheight./L,'k');
ylabel('z/L_a')
set(ax,'ylim',[-20 20])
dynamicDateTicks([], [], 'mm/dd');

ax = subplot(3,1,3);
h = plot(doy,[CD CHE]);
ylabel('Transfer Coeff.')
set(h,{'linestyle'},{':'; '-'})
legend('C_D', 'C_H, C_E','Location','BestOutside')
dynamicDateTicks([], [], 'mm/dd');

%% figure(4)
figure('name','MLD&Ln')
set(gcf,'position',goodfigsize(gcf))
ax = subplot(3,1,1);
h = plot(doy,MLD,'-');
set(ax,'ydir','rev')
legend('MLD','Location','BestOutside')
ylabel('h [m]')
dynamicDateTicks([], [], 'mm/dd');

ax = subplot(3,1,2);
h = semilogy(binave(doy,binsize),binave(Ln,binsize),'k',doy,ones(size(doy)),'k:'); % binsize is the average window
ylabel('Lake Number')
if exist('Ln','var')
set(ax,'ylim',[.00005 50000])
set(gca,'YTick',[10.^-4 10.^-3 10.^-2 10.^-1 10.^0 10.^1 10.^2 10.^3])
end
legend('L_N','Location','BestOutside')
dynamicDateTicks([], [], 'mm/dd');

ax = subplot(3,1,3);
h = plot(doy,atn,'k');
ylabel('atn [m^-^1]')
xlabel('Day')
legend('atn','Location','BestOutside')
dynamicDateTicks([], [], 'mm/dd');

%% Set common limits in all figures

FHS = getfighdls(1:4);

n = 0;
for i = 1:length(FHS)
    for j = 1:size(FHS(i).handles.Axes)
    n = n+1;
    allax(n) = FHS(i).handles.Axes(j);
    xlimall(n,:) = get(FHS(i).handles.Axes(j),'xlim');
    end
end

xlim = [min(xlimall(:)) max(xlimall(:))];

for i = 1:length(FHS)
    %FRAM: added if statement.  Can't find position variable
    if exist('figpos','var')
        set(FHS(i).handles.Fig,'position',figpos)
    end
    set(FHS(i).handles.Axes,'yminortick','on','xminortick','on','xlim',xlim)
end

eqax(1:4);
linkaxes(allax,'x');

