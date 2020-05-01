function seml(Sinfo)
% 
% VERSION:
%      The original version of this code was developed by Sally MacIntyre,
%      University of California, Santa Barbara (1983).  Subsequent
%      improvements  and alterations have been made under her supervision.  
%      1.X -- S. MacIntyre, 1983 (APL and Lotus)
%      2.X -- J. Romero, 1998 (FORTRAN)
%      3.0 -- C. Helmle, 8 MAR 2005
%      3.1 -- C. Helmle, 11 May 2005: Updated Sinfo structure; Changed
%      CDN10 from 1.0e-3 to 1.3e-3;  Implemented Charnock's WS10
%      calculation
%      3.2 -- J. Fram, 30 May 2007: Correcdcted HETFLX algorithm based on
%      finding error in Imberger (1985). Note, all runs of seml prior to
%      this date are incorrect.\
%      3.3 -- A. Cortes, Nov 3, 2014: Corrected CHE and LW and SW
%      calculated
%      3.4 -- A. Cortes, October 2016: Corrected kSR and kSR_600
%      3.5 -- A.Cortes, January 2018: Made multiple changes to make the
%      code avaible to other users. Search for ACC Jan 2018 to see specific
%      comments. Main changes: 
%           * Compute density as a function of temperature and salinitiy;
%           * The user can upload a specific conductance matrix  (SC) along with T in the raw input file, 
%           or specify a mean spcond value in the GUI
%           * Let the user choose the factor to compute: Sal(g/kg) = f * SC(mS/cm)
%           * Add buoyancy frequency variable to outputs;
%           * Ask the user the type of temperature loggers used, and select
%           the method to compute Mixed Layer Depth accordingly
%           * Define the filtering window of Lake numbers
%           * Modify some figures in outputs
%
% INPUT: (All fields shown are required; Sinfo is returned from seml_gui.m)
%
%      Sinfo = 
%           name: 'Lake Name'
%            lat: <value> Latitude (degrees)
%           elev: <value> Lake elevation (m)
%        sheight: <value> Wind sensor height (m)
%            met: 'Path to met data'
%         tchain: 'Path to tchain data' (raw)
%           bath: 'Path to hypsographic data'
%         spcond: <value> mean specific conductance value in the water in
%                   uS/cm. Default = 80
%           fsal: <value> factor to convert spcond (mS/cm) to salinity
%                   (g/kg). Default = 0.9
%     loggertype: 'RBR' or 'HOBO'
%             dT: <value> if RBR. Default = 0.02 oC
%         drhodz: <value> if HOBO. Default = 0.003 kg/m3/m
%
% NB! atn: Must be a time series in the met.mat input file
%
% 
% MAT file variable requirements:
%      File                    Variable     Size          Comments
%      -------------------------------------------------------------------
%      <METEOROLOGY.mat>:      doy          Vertical
%                              AirT         size(doy)
%                              RH           size(doy)
%                              SWin         size(doy)
%                              SWout        size(doy)     Optional - If none, it will be NaN                             
%                              LWin         size(doy)
%                              LWout        size(doy)     Optional - If none, it will be NaN 
%                              WS           size(doy)
%                              WDir         size(doy)     Optional - If none, it will be NaN 
%                              atn          size(doy)     Optional - If none, ask a mean value
%      <THERMISTORS.mat>:      time         Vertical      May be different than doy
%                              depth        Horizontal
%                              T            [time x depth]
%   OPTIONAL FILES:
%      <HYPSOGRPHY.mat>:       A            Vertical      Area
%                              d            size(A)       depth - May be different than depth
%
% OUTPUT: 
%      Name      Units     Description
%      -------------------------------------------------------------------
%      doy       Days      Day of year
%      AirT      deg C     Air Temperature
%      RH        %         Relative Humidity
%      WS        m/s       Wind Speed
%      WDir      o         Wind Direction
%      SWin      W/m^2     Shortwave Radiation In
%      SWout     W/m^2     Shortwave Radiation Out
%      SW0       W/m2      Shortwave Radiation Net
%      SWML      W/m2      Shortwave Radiation in the Mixed Layer
%      LWin      W/m^2     Longwave Radiation In
%      LWout     W/m^2     Longwave Radiation Out
%      LWnet     W/m^2     Longwave Radiation Net
%      MLD       m         Mixed Layer Depth
%      SrfT      deg C     Surface Temperature
%      atn       m^-1      Attenuation Coefficient
%      rhoair    kg/m3     Air density
%      ustair    m/s       Horizontal Frictional Velocity of the Air
%      usth2o    m/s       Horizontal Frictional Velocity of the Water
%      WS10      m/s       Velocity at 10 m (Law of the Wall)
%      wstar     m/s       Vertical Frictional Velocity of the Water
%      SE        W/m^2     Sensible Heat Exchange at the Air/Water Interface
%      LE        W/m^2     Latent Heat Exchange at the Air/Water Interface         
%      L         m         Monin-Obukhov Length Scale - Atmosphere          
%      Tau       kg/s^2/m  Shear Stress at the Air/Water Interface
%      CD        --        Momentum transfer coefficient
%      CHE       --        Heat and mass transfer coefficient
%      Ln        --        Lake Number
%      St        kg m      Overall Stability
%      z_t       m         Center of 'density gradient'
%      SRFFLX    W/m^2     Surface Energy Flux
%      HETFLX    W/m^2     Effective Heat Flux into the Active Mixing Layer
%      HETFLX_tot W/m^2    Total Effective Heat Flux into the water colum
%      T         deg C     Thermistor String Data
%      depth     m         Depth of Each Thermistor
%      SC        mS/cm     Matrix of specific conductance - Optional
%      Sal       g/kg      Matrix of Salinity - Optional
%      fsal                Factor of conversion: Sal (g/kg) = f * SC(mS/cm) (included in Sinfo)
%      RHO       kg/m3     Water density, computed with freshwater_density(T,Sal)
%      nsq       s^(-2)    Buoyancy frequency
%      n_cph     cph       Buoyancy frequency in cycles per hour
%      Sinfo     --        Updated Information Structure
%
       
% COORDINATE SYSTEM:
%      Vertical (+) upwards, hence heat gain positive, heat loss negative.
%
% SERVICE FUNCTIONS:
%      1) ATPRES   - atmospheric pressure calculation
%      2) SATVAP   - saturation vapor pressure calculation
%      3) VAPPRS   - vapor pressure calculation
%      3) SPCHUM   - specific humidity calculation
%      4) AIRDEN   - air density calculation
%      5) H2ODEN   - water density calculation
%      6) ALBEDO   - calculates albedo and refraction angle
%      7) TRNSFR   - calculates air column stability correction
%                    of neutral surface transfer coefficients
%evalin('base','clear');

h = waitbar(0);
% % Load meteorological data %%%
waitbar(.33,h,'Loading Meteorological Data')
load(Sinfo.met);
doy = vert(doy);
close(h);

%% ... Create missing vectors if necessary %%%
% -------------------------------------------------------------------------

% ... ACC Dec 2014: added LWin statement. We should include LWin
% calculation from cloud cover. 

if ~exist('LWin','var')
    LWin = nan*ones(size(doy));
end
if ~exist('LWout','var')
    LWout = nan*ones(size(doy));
end
if ~exist('WDir','var')
    WDir = nan*ones(size(doy));
end
if ~exist('WDstd','var')
    WDstd = nan*ones(size(doy));
end
if ~exist('atn','var')                               
   atn = repmat(str2num(Sinfo.atn),length(doy),1);  
end

%% ... Load T and average it. Load also SC matrix (if available)
%...... Otherwise, estimate from a mean value. Compute salinity as Sal (g/kg) = fsal * SC (mS/cm)
% -------------------------------------------------------------------------

% % Load Temperature %%
h = waitbar(0);
waitbar(.67,h,'Loading Thermistor Data')
load(Sinfo.tchain); T=double(T); 

%isotemps=double(isotemps); isotherm=double(isotherm); ACC Jan 2018 - Raw
%temperature data is accepted as an input. Isotherms are not needed
if any(isnan(T(:)))
    pernan = sum(any(isnan(T')))/length(T)*100;
    qda = questdlg([{'NaNs in thermistor data will cause NaNs in calculated results.'};...
        {''}; {['(' num2str(pernan,'%3.0f') '% NaN rate in current thermistor array)']};...
        {''}; {'Press ''OK'' to continue or ''Cancel'' to stop SEML'};],'Warning','OK','Cancel','OK');
    switch qda
        case 'Cancel'
            return
    end
end
time = vert(time);
depth = vert(depth)';
dtime = nanmedian(diff(time)); %diff(time(1:2)); %FRAM change
ddoy  = diff(doy(1:2));
close(h);

% ----- ACC Jan 2018
% Load values requested to the user 
spcond = (Sinfo.spcond)./1000; %mS/cm
fsal = (Sinfo.fsal);

% ----- ACC Jan 2018
% % Compute Salinity from SC (available or mean value)
if exist('SC','var')  
   Sal = fsal.* SC;
else % uniform matrix of SC using the mean input value (in uS/cm)
   SC = repmat(spcond,length(time),length(depth));  
   Sal = fsal.*SC;
end

% % Block average data if necessary and compute salinity from spcond
%
Sinfo.binave = 'T not block averaged';
if ddoy/dtime > 5                                   % If met interval
%     T = BANDPASS(T,[0 1/ddoy],dtime);             % is > 5x t-chain
  for j = 1:length(depth)                           % ACC Jan 2015, added to binaverage all depths
    Tn(:,j) = binave(T(:,j),floor(ddoy/dtime));                % interval
  end
  clear T
  T = Tn; clear Tn
  
  % ----- ACC Jan 2018 - If SC is available, average it as well
  for j = 1:length(depth)                           
      SCn(:,j) = binave(SC(:,j),floor(ddoy/dtime));
      Saln(:,j) = binave(Sal(:,j),floor(ddoy/dtime));
  end
    clear SC Sal
    SC = SCn; clear SCn
    Sal = Saln; clear Saln

    Sinfo.binave = ['T block averaged from ' num2str(dtime*24*60*60,'%4.0f')...
      ' seconds to ' num2str(floor(ddoy/dtime)*dtime*24*60*60, ...
      '%4.0f') ' seconds'];
    time = binave(time,floor(ddoy/dtime));
end

% % Interpolate to doy %%
% ACC Jan 2018 added SC and Sal variables
if time(1)>doy(1) && time(2)>doy(1)
    T = interp1([time(1)*2-time(2);time],[T(1,:);T],doy); %added Aug 5,2008
    SC = interp1([time(1)*2-time(2);time],[SC(1,:);SC],doy);
    Sal = interp1([time(1)*2-time(2);time],[Sal(1,:);Sal],doy);
else
    T = interp1(time,T,doy);
    SC = interp1(time,SC,doy);
    Sal = interp1(time,Sal,doy);
end


%% Compute density = f(T,Sal), and buoyancy frequency (nsq) % ACC Jan 2018
% -------------------------------------------------------------------------

% % Compute density as a function of temperature and Salinity
RHO = freshwater_density(T,Sal);%% Salinity in g/kg or psu, computed from SC in mS/cm and fsal

% % Compute the buoyancy frequency (nsq) and save it into the outputs
mean_rho = nanmean(nanmean(RHO));
[drhodz] = gradient(RHO,depth,0.3);
nsq = 9.81*drhodz/mean_rho; %s^-2
n_cph = real(sqrt(nsq))*3600/(2*pi);%cph


%% ... Calculate Mixed Layer Depth %%% % ACC Jan 2018
% -------------------------------------------------------------------------

disp(['T-Logger type = ',(Sinfo.loggertype)])
loggertype = Sinfo.loggertype;
SrfT = T(:,1);
SrfSC = SC(:,1);

% % Select the method to compute Mixed Layer Depth depending on the
% temperature logger type.

% NB. The user must look at the temperature and density data and determine
% the suitable temperature difference (dT) or density gradient (drhodz) to
% compute mixed layer depth
% Here, the GUI sets some default values of the two variables, BUT THOSE VALUES MUST BE CHECKED BY THE
% USERS THEIR OWN DATA

switch loggertype
    case 'RBR' % more accurated T data
        % .... Method: Temperature Difference from the surface approach (dT = 0.02oC, MacIntyre et al 2002)
        MLD = vert(ml_depth(doy,T,depth,Sinfo.dT)); 
        %%..Anavilhanas: MLD = [5,10,15]cm
        %if Sinfo.name == 'Anavilhanas'
        %   MLD = 0.15.*ones(size(doy));
        %   disp('Anavilhanas MLD = 5cm')
        %end
        disp(['Temperature difference (oC) = ',num2str(Sinfo.dT)])
    case 'HOBO' % less accurated T data. 
        % .... Method: Density gradient approach as in Rueda et al 2007.
        % The drhodz is lake specific (e.g. 0.003kg/m3/m)
        % RECOMMENDATION: Calibrate HOBO sensors to 0.05oC and
        % intercalibrate time series when the water column should be mixed.

        MLD = vert(ml_rueda(doy,RHO,depth,Sinfo.drhodz));
        disp(['Density gradient (kg/m3/m) = ',num2str(Sinfo.drhodz)])

end

%% ... Assign Parameter Values %%%
% -------------------------------------------------------------------------
Cn     = 1.33;       CDN10m = 1.3e-3;     CHEN10 = 1.35e-3;
emsvty = 0.97;       StfBlt = 5.67e-8;    LWa1   = 0.17;          
LWa2   = 0.642;      LWAlb  = 0.03;       ToK    = 273.15;         
TvCnst = 0.6078;     k      = 0.4;        g      = 9.81;
Lta1   = 2.5008e6;   Lta2   = 2.3e3;      Lcnst  = 0.61;
Cp_air = 1004.0;     Cp_h2o = 4184.0;     nuair  = 1.5e-5; 
ac     = .018;       z10    = 10;
z      = Sinfo.sheight;
Lat    = Sinfo.lat;
LkElev = Sinfo.elev;
%%% Cp_h20=4184 for ~15degC water, 4181 for ~20 degC 

%% ... Assign attenuation data to specific wavelengths %%%
% -------------------------------------------------------------------------
[SWFr Atn] = vardist([...
    0.04  26;     % UV band (280-400 nm)
    0.45  NaN;    % PAR band (400-700 nm) -- Inserted below
    0.13  1.1;    % 700-800 nm
    0.09  3.4;    % 800-900 nm
    0.04  26;     % 900-1000 nm
    0.25  870;    % 1200-1800 nm
    0.02  7800]); % 1800-2800 nm
%Sum(SWFr) = 1.02 -RLS 6/19/07
Atn = repmat(Atn',length(doy),1);
SWFr = repmat(SWFr',length(doy),1);
Atn(:,2) = atn;

%% .... Set a Wind Speed Threshold %%%
% -------------------------------------------------------------------------

% ACC Jan 2018. Treshold condition due to the resolution of the wind speed
% sensor (~0.5 m/s)

WSthresh = 0.5; % wind speed resolution
dumWS = find(WS <= WSthresh);
if ~isempty(dumWS)
    WS(dumWS) = 0.1; % minimum values. This number can change.
end


%% .... Calculate Atmospheric Condition Values %%%
% -------------------------------------------------------------------------
AtPr   = ATPRES(LkElev);
sphair = SPCHUM(AirT,AtPr,RH); %specific humidity.  Agrees with air-sea toolbox
vpair  = VAPPRS(AtPr,sphair);%Vapor pressure in air (mbar)
rhoair = AIRDEN(AtPr,AirT,sphair); % air_densin air-sea and 
Tvair = (AirT+ToK).*(1+0.6078*sphair); %VIRTUAL TEMPERATURE OF AIR

sphh2o = SPCHUM(SrfT,AtPr,100.0); %specific humidity at saturation.  Agrees with air-sea toolbox qsat
svp = SATVAP(SrfT,AtPr);%Saturation Vapor pressure over the plane air-water (mbar)
Tvwater = (SrfT+ToK).*(1+0.6078*sphh2o); %VIRTUAL TEMPERATURE OF SURFACE WATER
rhoh2o = freshwater_density(SrfT,SrfSC.*0.9);
[Albdo RAng] = Albedo(doy,Lat);

%% .... Radiation calculations %%%
% -------------------------------------------------------------------------
% Longwave Radiation Calculations 
% ACC, March 2017 - We force the code to compute LWout as a function of SurfT, to avoid problems
% with condensation
LWout = StfBlt*(SrfT+ToK).^4;
LWout(isnan(LWin)) = nan;  % Don't calculate LWout if no LWin

% Calculate SWout if not provided %%
if ~exist('SWout','var'); 
    % FRAM swapped Albdo for 0.05 5/31/07
    Albdo(Albdo > 0.35) = 0.35; % near horizon reflection is probably overestimated
    SWout = SWin.*Albdo;
    
   elseif any(isnan(SWout))
     %SWout(isnan(SWout)) = 0.05*SWin(isnan(SWout)); %Jose
     SWout(isnan(SWout)) = Albdo(isnan(SWout)).*SWin(isnan(SWout));
end

% Computed values  of LW and LWnet assumptions when no data is available
LWOut  = - LWout; 
LWNet  = LWin + LWOut; 

if any(isnan(LWNet))
    disp('Assuming LWnet = 100 W/m2');%ACC 2015
    dum = find(isnan(LWNet));
    LWNet(dum) = -100; % ACC 2015
end


%% .... Air Column Stability Loop [Hicks, 1975] %%%
% -------------------------------------------------------------------------
% - Account for sensor height deviation from 10m standard with formula from
%   Amorocho and Devries (1980: JGR 85 (NC1): 433-442)

LtHeat = Lta1-Lta2*AirT; %heat to vaporize. One needs to get it to 100oC and then vaporize. From vapor.m in air-sea toolbox
Tv     = (AirT+ToK).*(1+TvCnst*sphair);
count  = zeros(length(AirT),1);
L      = zeros(length(AirT),1);

CDN = CDN10m/(1-CDN10m^0.5/k*log(10.0/z))^2.0;  % FRAM confirmed that this is right 6/4/07. 
CHEN = CHEN10/(1-CHEN10^0.5/k*log(10.0/z))^2.0;  % ACC 11/03/14. 

[Tau ustair usth2o SE LE L CD CHE] = deal(NaN*ones(length(doy),1));
h = waitbar(0);
for i = 1:length(AirT)
    waitbar(i/length(AirT),h,['Processing Energy Budget Data: '...
            num2str(round(i/length(AirT)*100)) ' %'])
    init = true;
    cont = true;
    CHE(i) = CHEN;                         % Mass and heat transfer coefficients
    CD(i)  = CDN;                          % Momentum transfer coefficients
    while cont
        Tau(i)    = rhoair(i)*CD(i)*WS(i)^2;
        ustair(i) = (Tau(i)/rhoair(i))^0.5;
        usth2o(i) = (Tau(i)/rhoh2o(i))^0.5;
        SE(i)     = rhoair(i)*Cp_air*CHE(i)*WS(i)*(AirT(i)-SrfT(i));
        LE(i)     = rhoair(i)*LtHeat(i)*CHE(i)*WS(i)*(sphair(i)-sphh2o(i));
        L(i)      = rhoair(i)*ustair(i)^3*Tv(i)/k/g; % SE, LE, and L are from Hicks 1975
        L(i)      = L(i)/((SE(i)/Cp_air)+(Lcnst*(AirT(i)+ToK)*LE(i)/LtHeat(i))); %Monin-Obukhov Length
        [CD(i) CHE(i)]  = TRNSFR(z,L(i),CDN,CHEN); % from Paulson 1970, ACC 11/03/14
        
        if init
            OldL = L(i);
            init = false;
            continue;
        end
        count(i)  = count(i) + 1;
        if (OldL - L(i))/OldL > 1e-5
            OldL = L(i);
        else
            cont = false;
        end
    end
end
close(h);
[count Tau ustair usth2o SE LE L CD CHE] = ...
    vardist(real([count Tau ustair usth2o SE LE L CD CHE]));

%% ... Calculate Lake Number %%%
% % Only if hypsographic data is available % %
% -------------------------------------------------------------------------

load(Sinfo.bath)

d  = vert(d); A = vert(A);
[d ix] = sort(d); A  = A(ix);
A  = interp1(d,A,vert(depth))';
% ACC Jan 2018 changed to mylakenumber_vic. Now St is an output
[Ln,z_t,St,z_g,drhodz] = mylakenumber_vic(depth,A,RHO,Tau); % St = stability

% % Define the binsize to average Lake  numbers % ACC Jan 2018 
dt = floor(nanmean(diff(doy)*24*60)); % min
Tiw = 24; % Default => 24 h. Period of the dominant internal wave mode. In small lakes, it can be ~4h
Tiw4th = Tiw/4; % Default => 6h. In small lakes ~ 1h
binsize = floor(Tiw4th*60/dt);
Sinfo.Tiw4th = Tiw4th;
disp(['Ln binsize = 1/4 of the period of the dominant internal wave mode (h) = ',num2str(Tiw4th)])

%% ... Compute the heat flux of the mixed-layer [Imberger, 1985] % % %
% -------------------------------------------------------------------------
SW0    = SWin - SWout;
SW1    = SW0.*cos(RAng).*sum(SWFr./Atn,2);
SW2    = SW0.*cos(RAng).*...
         sum(SWFr./Atn.*exp(-Atn.*repmat(MLD./cos(RAng),1,size(Atn,2))),2);
SWML   = SW0.*sum(SWFr.*exp(-Atn.*repmat(MLD./cos(RAng),1,size(Atn,2))),2);

SRFFLX = LE + SE + LWNet;

%% ... Effective surface heat flux [Kim, 1976] & [Rayner, 1980] % % %
% -------------------------------------------------------------------------

HETFLX = SRFFLX + SW0 - SWML;
HETFLX_tot = SRFFLX + SW0; % ACC on March 2015

%% ...Calculate WS10 from Smith [1988, JGR-Oceans, pg15467]
% -------------------------------------------------------------------------

zc   = ac*ustair.^2/g;          %Charnock. sm. 5aug2018. We are using 0.018, for coastal waters. Open is 0.011.
zs   = 0.11*nuair./ustair;      %Businger
    %0.11 is correct for smallish lakes.  See Ian Jones p.151 -FRAM 2/1/06
zo   = zc+zs;
% stable. see TRNSFR.
stab=(z./L>=0 & z./L<.5);
if ~isempty(stab); phi(stab) = -5*z10./L(stab); end
stab=(z./L>0.5 & z./L<10);
if ~isempty(stab); phi(stab)=0.5./(z./L(stab)).^2.0-4.25./(z./L(stab))-7*log(z./L(stab))-0.852; end
stab=(z./L>10 & z./L<15);
if ~isempty(stab); phi(stab)=log(z./L(stab))-0.76*(z./L(stab))-12.093; end
% unstable
% x    = (1-16*z10./L(~stab)).^(.25); was this
stab=(z./L>=0);  % need to choose all stable values
x    = abs((1-16*z10./L(~stab)).^(.25)); % i^4==1 and 1^4==1, so need abs
phi(~stab) = 2*log((1+x)/2) + log((1+x.^2)/2)-2*atan(x)+pi/2;
phi  = vert(phi);
WS10 = ustair/k.*(log(z10./zo)-phi);

%% Print Figures % % % 
% -------------------------------------------------------------------------
assignin('base','Sinfo',Sinfo);
plot_seml

%% Write Variables to matfile % % %
% -------------------------------------------------------------------------
% ACC Jan 2018 increased the relevant outputs

Sinfo.rundate = datestr(now);
savvars = {'doy'     'AirT'    'RH'      'WS'        'SWin'    'SWout'   ...
           'LWin'    'LWout'   'LWNet'   'MLD'       'atn'     ...
           'WDir'    'SW0'     'SWML'    'WS10'      'rhoair'            ...
           'SrfT'    'SE'      'LE'      ...
           'L'       'Tau'     'SRFFLX'  'HETFLX'   'HETFLX_tot' ...
           'CD'      'CHE'     ...
           'Ln'      'z_t'     'St'      'A'         'd'   ...
           'T'       'depth'   'Sinfo'   ...
           'RHO'     'nsq'     'n_cph'   'SC'        'Sal'     'binsize' ...
           'AtPr'};

for i = 1:length(savvars)
eval(['assignin(''base'',''' savvars{i} ''',' savvars{i} ');']);
end

%cd seml
Sinfo.filename = uisave_ch(savvars,Sinfo.name,'Save SEML results');
%cd ..

disp(Sinfo)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function atpres = ATPRES(LkElev)
% CALCULATE ATMOSPHERIC PRESSURE
% This calculation is based on Berry et al (1945) Handbook of Meteorology
%     P = 1013.25x(1-0.0065xelev/288)^5.25587
% where P is in mb and elev is in m

a1=1013.25; a2=0.0065; a3=288; a4=5.25587;
atpres = a1*(1-a2*LkElev/a3)^a4;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function satvap = SATVAP(T,P)
% CALCULATE SATURATION VAPOR PRESSURE
% This calculation of saturation vapor pressure over a plane water surface
% is based on SMITHOSONIAN METEOROLOGICAL TABLES for temperatures +/- 40øC
% (TABLE 94), and correction for atmospheric conditions (TABLE 89)
%     satvap = fw*es
% where es = 10^[(0.7859+0.03477T)/(1+0.00412T)] = ideal sat vap
%       fw = 1+10E-6xPx(4.5+0.0006T = sat vap correction
%       satvap = corrected sat vap
% Units millibars where es correct to 1 part in 500, and fw 2 parts in 10E4
% See Gill (1982) Atmosphere-Ocean Dynamics
% VARIABLES:
%     T   Temperature (øC)
%     P   Atmospheric Pressure (mb)

a1=0.7859;  a2=0.03477;  a3=0.00412;
b1=1.0e-6;  b2=4.5;   b3=0.0006;

es = 10.^((a1+a2*T)./(1+a3*T));
fw = 1+(b1*P*(b2+b3*T.^2));
satvap = es.*fw;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function spchum = SPCHUM(T,P,RH)
% CALCULATE SPECIFIC HUMIDITY
% This calculation of specific humidity is based on the following:
% #1  vp/P=sph/(0.622+(1-0.622)xsph)
% rearranging
% #2  sph = (0.622xvp)/(P-0.378xvp)
% #2 used to calcualte saturation specific humidity (sphsat) from
%     saturation vapor pressure (svp)
% specific humidity related to relative humidity (RH) as,
% #3  RH = r/rw
% where r is the mixing ratio equal to sph/(1-sph) and rw is the
%     saturation mixing ratio which yields
% #4  RH = sph x (1 - sphsat) / ((sphsat x (1 - sph))
% rearranging
% #5  spchum = RH x sphsat / (1 - sphsat + RH x sphsat)
% #5 used to calculated specific humidity (sph) from saturation
%     specific humidity (sphsat)
% INPUTS and UNITS:
%     svp     mb
%     P       mb
%     RH      %

a1=0.62197;

svp = SATVAP(T,P);
sphsat = a1*svp./(P-(1-a1)*svp);
spchum = RH/100.*sphsat./(1-sphsat+RH/100.*sphsat);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function vapprs = VAPPRS(P,sph)
% COMPUTE VAPOR PRESSURE (from sph)

epslon = 0.62197;
vapprs = P*sph./(epslon+(1.0-epslon)*sph);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function airden = AIRDEN(P,T,sph)
% CALCULATE AIR DENSITY
% This calculation is based on equation of state for air:
%     airden = P/(R*Tv)
% where R is the universal gas constant = 287.04 J/kg/øK
%     and Tv is the virtual temperature given as
%     Tv = Tx(1+0.6078xsph)
% INPUTS and UNITS:
%     sph     no units
%     T       C  NOTE: must convert to units of øK
%     P       mb  NOTE: must convert to units of Pa(N/m2) by multiplying by 100
%                       (i.e. 1 mb = 100 Pa)

a1=0.6078;	R=287.04;	ToK=273.15;
Tv = (T+ToK).*(1+0.6078*sph);
airden = (P*100)./(R*Tv);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function h2oden = H2ODEN(T)
% CALCULATE WATER DENSITY
% Compute density of water based on temperature only, IES 80
% Millero and Poisson(1981)

a0=999.842594;  a1=6.793952e-2; a2=-9.09529e-3;
a3=1.001685e-4;   a4=-1.120083e-6;   a5=6.536332e-9;

h2oden = a0 + a1*T + a2*T.^2 + a3*T.^3 + a4*T.^4 + a5*T.^5;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [CD,CE] = TRNSFR(z,L,CDN,CEN)
% CALCULATE TRANSFER COEFFICIENTS
% Compute transfer coefficients based on the integrated similarity 
% functions of Paulson (1970)
% Paulson, C. A.: 1970, The Mathematical Representation of Wind Speed
% and Temperature Profiles in the Unstable Atmospheric Surface Layer,
% J. Appl. Meteorol. 9, 857–861.
% Also, see Panofsky and Dutton 1984 book "Atmospheric Turbulence"
% and Cheng, Y and Brutsaert, W: 2005, Flux-profile relationshipts for wind 
% speed in the stable atmospheric boundary layer. BL Met. 114, 519-538.
% and Cheng, Y., M.B. Parlange, and W. Brutsaert: 2005. JGR, 110, D06101.
k=0.4;  % von Karman number

% compute integrated similarity functions
if (z/L) > 15.0; L = z/15.0; end
if (z/L) < -15.0; L = -z/15.0; end
if (z/L) < 0 %Paulson 1970 eqn 4. Unstable
    x = (1-16*(z/L))^0.25;
    sm  = 2.0*log((1.0+x)/2.0)+log((1+x^2.0)/2.0)-2.0*atan(x);
    sm  = sm + pi/2.0;
    se  = 2.0*log((1+x^2.0)/2.0);
elseif (z/L) < 0.5 % Stable.  Track down Hicks 1976 from the library.
    sm = -5*(z/L); % section 6.5.2 Panofsky eqn 9.
    se = sm;
elseif (z/L) < 10.0
    sm = 0.5/(z/L)^2.0-4.25/(z/L)-7*log(z/L)-0.852; %This is from Chad's
    se = sm;  % class notes via Todd Cowen. Cheng and Brutsaert give 
    % comparable equations, but let's stick with Chad's for now. -FRAM 6/5/07
else
    sm = log(z/L)-0.76*(z/L)-12.093; % Ditto
    se = sm;
end

% compute transfer coefficients
% Where are these from? - FRAM 6/4/07
CD = CDN/(1+CDN/k^2.0*(sm*sm-k*sm/CDN^0.5-k*sm*CDN^0.5/CDN));
CE = CEN/(1+CEN/k^2.0*(sm*se-k*se/CDN^0.5-k*sm*CDN^0.5/CEN));

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [Albdo RefAng] = Albedo(doy,Lt)
% This program calculates the albedo on a flat lake from Fresnal's
% Law (Neumann & Pierson, 1966).  The angle of refraction is given
% by Snell's Law.  The index of refraction is assumed to be
% independent of wavelength.  The Zenith angle is calculated from
% Milankovich (1930).

% Assign Parameter Values %
Tropic = 23.45; YrDays = 365.0;       DclDay = 284.0;
DgCrcl = 360.0; DgToRd = 57.29577951; RefInd = 1.34;

Lat = Lt/DgToRd;
doy = doy-1;                                    % Start year at 0.0
z = DclDay + floor(doy);
x = DgCrcl*z/YrDays/DgToRd;
y = sin(x);
Decl = Tropic/DgToRd*y;

% Hour-angle calculation where time is hours from midnight  %
HrAng = ((doy-floor(doy))*24 -12)*15.0/DgToRd;

% Zenith angle calculation %
Zenith = acos(sin(Decl)*sin(Lat)+cos(Decl)*cos(Lat).*cos(HrAng));

% Angle of Refraction calculation based on Snell's Law %
RefAng = asin(sin(Zenith)/RefInd);

% Albedo Calculation %
A1 = tan(Zenith-RefAng).^2;
A2 = tan(Zenith+RefAng).^2;
A3 = sin(Zenith-RefAng).^2;
A4 = sin(Zenith+RefAng).^2;
Albdo = 0.5 * (A1./A2 + A3./A4);

% Set albedo to 1 if greater than 1 %
Albdo(Albdo>1) = 1;

% Index doy back to original value %
doy = doy +1;




