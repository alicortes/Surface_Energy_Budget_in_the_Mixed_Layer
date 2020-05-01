% Feed this function hourly data that you want 5 min data from and a 5 
% minute data whose spectral properties you want the created 5 minute data
% to look like.
% Output is 5 minute wind data for the 1 hour time step in the form of a
% structure named out (out.WS, out.WDir, out.WN, and out.WE)
%
% FRAM March 2007
%
function out=windHrTo5min(doyH,WSH,WDirH,doy5,WS5,WDir5)

verbose = 0; % Plot spectra as opposed to just produce stuff

land.doy=doyH; % Lake and land refer to toolik where we always have hourly
lake.doy=doy5; % land data and we sometimes have 5 minute lake data.
land.WN=WSH.*cos(WDirH*pi/180);
land.WE=WSH.*sin(WDirH*pi/180);
lake.WN=WS5.*cos(WDir5*pi/180);
lake.WE=WS5.*sin(WDir5*pi/180);

for s=1:2
    if s==1; ss='WN'; else ss='WE'; end
    %% plot spectra
    y=eval(['lake.',ss,';']);
    if verbose
        hh=figure('name',[ss,' wind']);
        [PSD,CLV] =OpPSD(nanmedian(24*3600*diff(lake.doy)),y,1,1,64,1.08,1,0);
        loglog(PSD(:,1),PSD(:,2),'b','linewidth',2); hold on

        [PSD,CLV] =OpPSD(nanmedian(24*3600*diff(land.doy)),eval(['land.',ss]),1,1,32,1.08,1,0);
        loglog(PSD(:,1),PSD(:,2),'r','linewidth',2);

        lake.W_hr=kron(binave(y,12),ones(1,12));
        [PSD,CLV] =OpPSD(nanmedian(24*3600*diff(lake.doy)),lake.W_hr,1,1,64,1.08,1,0);
        loglog(PSD(:,1),PSD(:,2),'c','linewidth',2);
        loglog(PSD(:,1),PSD(:,2),'m','linewidth',2);
        loglog(PSD(:,1),PSD(:,2),'c','linewidth',2);

        legend('lake 5min','land hrly','lake hrly avg');
        vert_line([1/(24*3600) 1/(3600) 1/(60*5)]) ;
        title({'The hourly lake and land data match at high f,',' so I should be able to make approriate transfer function'});
    end
    % It's the same, so fill in high frequency wind to make land look like lake

    %% covariance. define semivariogram
    % We do this on WN and WE not WS because WS is not normally distributed.
    % Oddly, WS looks log normal to me at least in 2004.

    dt=median(diff(lake.doy))*24*60; %minutes
    maxlags=60/5;
    [C,lags]=xcov(y,y,maxlags,'coeff');
    if size(C,1)==size(lags,1);
        dex=find(lags>=0 & C>0); lags=lags(dex); C=C(dex);
    else
        dex=find(lags>=0 & C'>0); lags=lags(dex)'; C=C(dex);
    end
    p=polyfit(lags*dt,log(C),1);
    p(1)=p(1)*.5;
    % Use this to define correlation scale over an hour
    % A different correlation function might work better, such as a spherical
    % one, but exponential is good enough for now.  See Y. Rubin's book for
    % details on other appropriate correlation functions.  Not all functions
    % that fit the correlogram are usable.
    if verbose
        figure('name','correlation');
        semilogy(lags*dt/60,C,lags*dt/60,exp(p(1)*lags*dt));
%         title(['5min lake data: ',ss]);
        ylabel('cross-correlation coefficient');
        xlabel('lag in hours');
    end

    %% compute random variables with proper covariance
    % Method follows HWK#5 from CE202b taught by Y. Rubin
    y=eval(['land.',ss,';']);

    var_Y=var(y)*2; % variance x 1.1 to go from hourly to appx 5min variance
    meanabs_all=mean(abs(y));
    e_number=length(y);
    Q_size=12+1; % integer, # of nodes

    % Make Q matrix
    [i,j]=ndgrid(1:1:Q_size,1:1:Q_size);
    Q=var_Y*exp(p(1)*dt*abs(i-j)); % exponential correlation
    L_sym=chol(Q)';  % Matlab LU decomposition (actual "lu" function not symmetric)

    ii=randn(Q_size,e_number);  % random numbers generated.
    Y_est = L_sym*ii;   % matrix of estimated values formed.
    % mean_var=mean(var(Y_est))  matlab "var" command calc.'s a mean, thus biasing result
    % mean_mean=mean(mean(Y_est));

    if verbose
        % Compute Ensemble mean... of the covariances
        Ensemble_Y=ones(Q_size,Q_size);
        for i = 1:1:Q_size
            for j=1:1:Q_size
                Ensemble_Y(i,j)=1/e_number * sum(Y_est(i,:).*Y_est(j,:));
            end
        end
        % Compute covariances
        Y_cov_est=ones(1,Q_size)*NaN; node_dist=Y_cov_est;
        for i=0:1:Q_size-1
            Y_cov_est(i+1)= sum(diag(Ensemble_Y,i))/length(diag(Ensemble_Y,i));
            node_dist(i+1)=i*1000;
        end
        hhh=0:.1:Q_size-1;
        Y_cov_model=var_Y.*exp(p(1)*dt.*hhh);

        figure('name','do results fit the model?')
        plot(hhh*1000,Y_cov_model,'r',node_dist,Y_cov_est,'b*');
        legend('Exponential Covariance (Given)', 'Ensemble averaged Covariance')
        ylabel('Covariance'), xlabel('Distance');
    end

    % make "kriged" mean space
    out.doy=land.doy(1):5/24/60:land.doy(end)+1/24*55/60;
    kr=interp1(land.doy,y,out.doy,'nearest','extrap');

    % fill in data over each hour
    eval(['out.',ss,'=kr*NaN;']);
    for i=1:ceil(length(y)/12) %length(y)-1
        % disp(int2str([i length(y) size(Y_est,1) size(Y_est,2)]));
        % find line through end points
        tmp=Y_est(:,i);
        pp=polyfit([1,13]',tmp([1,13]),1);
        tmp=tmp-polyval(pp,(1:13)');
        %disp([i length(kr)])
        tmp=tmp(1:12)'+kr((i-1)*12+(1:12));
        eval(['out.',ss,'((i-1)*12+(1:12))=tmp;']);
        %   lake.WN_new((i-1)*12+(1:12))=tmp(1:12)'*mean(abs(y(i+(0:1))))/meanabs_all+lake.WN_kr((i-1)*12+(1:12));
    end
    tmp=Y_est(:,end);
    pp=polyfit([1,13]',tmp([1,13]),1);
    tmp=tmp-polyval(pp,(1:13)');
    tmp=tmp(1:12)'+kr(end+(-11:0));
    eval(['out.',ss,'(end+(-11:0))=tmp;']);
    %lake.WN_new(end+(-11:0))=tmp(1:12)'*abs(y(end))/meanabs_all+lake.WN_kr(end+(-11:0));

    if verbose
        figure(hh);
        [PSD,CLV] =OpPSD(nanmedian(24*3600*diff(lake.doy)),eval(['out.',ss]),1,1,64,1.08,1,0);
        loglog(PSD(:,1),PSD(:,2),'m');
        legend({'lake 5min','land hrly','lake hrly avg','land 5min created'});
    end
end
out.WS=sqrt(out.WN.^2+out.WE.^2);
th=atan(out.WE./out.WN)*180/pi;
dex=find(out.WN>0 & out.WE<0);
th(dex)=th(dex)+360;
dex=find(out.WN<=0 & out.WE<0);
th(dex)=th(dex)+180;
dex=find(out.WN<0 & out.WE>0);
th(dex)=th(dex)+180;
out.WDir=th;