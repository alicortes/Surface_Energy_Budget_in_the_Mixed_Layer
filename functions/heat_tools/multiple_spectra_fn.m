function multipl_spectra(time1,time2,
% multiple_spectra.m
% calculates spectral densities
% for multiple periods
% (linear detrending and hanning windows)
%
% Modified: BL, 04/2005 - to be called by other script for multiple panels

disp ('This program requires the following variables: ');
disp (' ');
disp ('   time       time vector          unit: days');
disp ('   isotherm   isotherm matrix      unit: m');
disp ('   isotemps   temperature vector   unit: C');
disp (' ');
file=input('Enter mat-file to load, or press enter: ','s');
if ~isempty(file)
	eval(['load ',file]);
end

% plot all data ----------------------------------------------------
 
figure(1);
hold off;
t=[];
delta_t=(max(time)-min(time))/length(time);	% interval between data points
s=size(isotherm);
for i=1:s(2)				% create time matrix            
	t=[t time];
end
plot(t,isotherm);
set (gca,'ydir','rev');
leg=num2str(isotemps);
legend(leg);
label=['Time (Julian Day)'];
xlabel(label);
ylabel('Depth (m)');
title ('Isotherms (ºC)');

% create signature for plot

xl=get (gca,'xlim');
yl=get (gca,'ylim');
x=(xl(1)+(xl(2)-xl(1))*.7);
y=(yl (2)+(yl(2)-yl(1))*.08);
txt=([pwd,'\',file,'.mat']);
   for f=1:length(txt)
      if txt(f)=='\'
         txt (f)='/';
      end
      if txt(f)=='_'
         txt(f)='-'
      end   
   end   
text (x,y,txt);
orient landscape;

% select intervals ----------------------------------------------- 

new_plot='y';
while new_plot=='y';
	time1=input('Enter beginning of interval (day): ');
	time2=input('Enter end of interval (day): ');
	t1=floor((time1-t(1))/delta_t);
	t2=floor((time2-t(1))/delta_t);

% plot interval --------------------------------------------------

	figure(2);
	hold off;
	plot(t(t1:t2,:),isotherm(t1:t2,:));
	set (gca,'ydir','rev');
  leg=num2str(isotemps);
	legend(leg);
 	xlabel(label);
	ylabel('depth (m)');
	title ('isotherms (C)');
   xl=get (gca,'xlim');
   yl=get (gca,'ylim');
   x=(xl(1)+(xl(2)-xl(1))*.7);
   y=(yl (2)+(yl(2)-yl(1))*.08);
   text (x,y,txt);
   orient landscape
   disp (['This interval contains ',num2str(t2-t1+1),' data points.']);
   window=input ('Enter size of interval for PSD (e.g. 1024, 2048...): ');
   windows=window;
   isoth_deg=input('Enter isotherm: ');	% in degrees C
   freq_lowers=.0001;
   freq_uppers=100000;
   
% change isotherm from degrees C to integer --------------------

	isoth=round((isoth_deg-isotemps(1))/(isotemps(2)-isotemps(1)))+1;

   conf_interval=input ('Enter confidence interval (0.0 - 1.0): ');
   overlap=input ('Enter degree of overlap: ');
   N1=input('Enter BV-Frequency (rad/s) or 1 if not available: '); 
   if N1~=1
      N=N1/2/pi;
   else
      N=N1;
   end   
   
   number_of_windows=0;  
   choice=' ';
   while choice~='n'
      number_of_windows=number_of_windows+1;
      
% calculate spectral density without overlap ----------------------------- 
% we do this to get confidence intervals (Pconf).

	detrend='linear';  % detrending
	[P,Pconf,freq]=psd(isotherm(t1:t2,isoth),window,1/delta_t,hanning(window),0,conf_interval,detrend);
  Pconf_N=Pconf./[P P];
 
% spectral density with overlap ------------------------------------------   
   
   [P,Pconf,freq]=psd(isotherm(t1:t2,isoth),window,1/delta_t,hanning(window),floor (overlap*window),conf_interval,detrend);
   disp ('Here, the confidence intervals for non-overlap are used.');
   P_N=P*N*86400;	% spectral density * N
   
% combine old and new PSD data
  
  if choice~=' ' 
     lower=max(find(freq_old<freq_lower));
     if ~isnan (lower)
        lower1=min(find(freq>=freq_lower));
        freq_new=freq_old(1:lower);
        P_N_new=P_N_old(1:lower);
        Pconf_N_new=Pconf_N_old(1:lower,:);
        freq_new=[freq_new' freq(lower1:length(freq))']';
        P_N_new=[P_N_new' P_N(lower1:length(freq))']';
        Pconf_N_new=[Pconf_N_new' Pconf_N(lower1:length(freq),:)']';
     else   
        freq_new=freq;
        P_N_new=P_N;
        Pconf_N_new=Pconf_N;
     end   
     upper=min(find(freq_new>=freq_upper));
     if ~isnan (upper)
        upper1=min(find(freq_old>freq_upper));
        freq_new=freq_new(1:upper);
        P_N_new=P_N_new(1:upper,:);
        Pconf_N_new=Pconf_N_new(1:upper,:);
        freq_new=[freq_new' freq_old(upper1:length(freq_old))']';
        P_N_new=[P_N_new' P_N_old(upper1:length(freq_old))']';
        Pconf_N_new=[Pconf_N_new' Pconf_N_old(upper1:length(freq_old),:)']';
     end   
     freq=freq_new;
     P_N=P_N_new;
     Pconf_N=Pconf_N_new;
  end       
  
% plot spectra and confidence intervals ----------------------------------- 

% BL - multiple panels;
if exist('ca','var'),
    figure(ca.fig);
    subplot(ca.panels,1,ca.ax);
else
    figure(3);
end

   hold off
   h1=loglog (freq,P_N,'b');
   hold on
   logticks;
   orient landscape
   all_P_N=P_N;
   
% plot confidence intervals at the bottom of the graph   
   
   ylimits=get (gca,'ylim');
   Pconf_N1=Pconf_N*ylimits(1)*10;
   loglog ([freq freq],Pconf_N1,'k');
   loglog (freq,ones(1,length(freq))*ylimits(1)*10,'k');
   
% plot slope of -2 and N ---------------------------------------------------

   loglog([max(freq) max(freq)/50]',[P_N(4)/2.5 P_N(4)*1000]'/100,'k');
	loglog([N*86400 N*86400]',[ymin P_N(4)/5]','k:');
	%txt=['N = ',num2str(N),' Hz'];
	%text ((N*86400)*1.1,ymin*3,txt);
	xlabel ('Frequency (cycles per day)');
   if N~=1
      ylabel ('Spectral Density * N (m^2)');
   else
      ylabel ('Spectral Density (m^2/Hz)');
   end   
  xl=get (gca,'xlim');
  yl=get (gca,'ylim');
  x=(xl(1)*(xl(2)/xl(1))^.7); 
  y=(yl (1)/(yl(2)/yl(1))^.08);
  text (x,y,txt);
 
   choice=input ('Calculate PSD for smaller windows? (y/n): ','s');   
   if choice=='y'
       window=input ('Enter size of window (e.g. 1024, 2048...): ');
       disp ('The results of this PSD will be plotted in a selected frequency range.');
       freq_lower=input ('Select lower frequency: ');
       freq_upper=input ('Select upper frequency: ');
       freq_old=freq;
       Pconf_N_old=Pconf_N;
       P_N_old=P_N;
       windows=[windows;window];
       freq_lowers=[freq_lowers freq_lower];
       freq_uppers=[freq_uppers freq_upper];
   end    
 end
 
 leg='                                  ';
 ttl=[num2str(isoth_deg),' C Isotherm, Day ',num2str(time1),' - ',num2str(time2)];
 leg(1:length(ttl))=ttl;
 legend(leg);
 
 %add additional spectra if desired
 
 new_spectra=input('Add another plot (y/n) ','s');
 number_of_spectra=1
 while new_spectra=='y'
    number_of_spectra=number_of_spectra+1;
	time1=input('Enter beginning of interval (day): ');
	time2=input('Enter end of interval (day): ');
	t1=floor((time1-t(1))/delta_t);
	t2=floor((time2-t(1))/delta_t);
   isoth_deg=input('Enter isotherm: ');	% in degrees C
	isoth=round((isoth_deg-isotemps(1))/(isotemps(2)-isotemps(1)))+1;
   
      N1=input('Enter BV-Frequency (rad/s) or 1 if not available: '); 
   if N1~=1
      N=N1/2/pi;
   else
      N=N1;
   end   

   for i=1:number_of_windows
      
   % spectral density with overlap ------------------------------------------   
   
   [P,Pconf,freq]=psd(isotherm(t1:t2,isoth),windows(i),1/delta_t,hanning(windows(i)),floor (overlap*windows(i)),conf_interval,detrend);
   disp ('Here, the confidence intervals for non-overlap are used.');
   P_N=P*N*86400;	% spectral density * N
   
% combine old and new PSD data
  
  if i>1 
     lower=max(find(freq_old<freq_lowers(i)));
     if ~isnan (lower)
        lower1=min(find(freq>=freq_lowers(i)));
        freq_new=freq_old(1:lower);
        P_N_new=P_N_old(1:lower);
        freq_new=[freq_new' freq(lower1:length(freq))']';
        P_N_new=[P_N_new' P_N(lower1:length(freq))']';
     else   
        freq_new=freq;
        P_N_new=P_N;
     end   
     upper=min(find(freq_new>=freq_uppers(i)));
     if ~isnan (upper)
        upper1=min(find(freq_old>freq_uppers(i)));
        freq_new=freq_new(1:upper);
        P_N_new=P_N_new(1:upper,:);
        freq_new=[freq_new' freq_old(upper1:length(freq_old))']';
        P_N_new=[P_N_new' P_N_old(upper1:length(freq_old))']';
     end   
     freq=freq_new;
     P_N=P_N_new;
  else
     
  end
 P_N_old=P_N;
 freq_old=freq;
end 
  
% plot spectra  ----------------------------------- 

hold off
all_P_N=[all_P_N P_N];
h1=loglog (freq,all_P_N);
hold on
   logticks;
   orient landscape

   ttl=[num2str(isoth_deg),' C Isotherm, Day ',num2str(time1),' - ',num2str(time2)];
   leg=[leg;'                                  '];
   leg(end,1:length(ttl))=ttl;
   legend(leg);
   
   % plot confidence intervals at the bottom of the graph   
   
   ylimits=get (gca,'ylim');
   Pconf_N1=Pconf_N*ylimits(1)*10;
   loglog ([freq freq],Pconf_N1,'k');
   loglog (freq,ones(1,length(freq))*ylimits(1)*10,'k');
   
% plot slope of -2 and N ---------------------------------------------------

   loglog([max(freq) max(freq)/50]',[P_N(4)/2.5 P_N(4)*1000]'/100,'k');
	loglog([N*86400 N*86400]',[ymin P_N(4)/5]','k:');
	%txt=['N = ',num2str(N),' Hz'];
	%text ((N*86400)*1.1,ymin*3,txt);
	xlabel ('Frequency (cycles per day)');
   if N~=1
      ylabel ('Spectral Density * N (m^2)');
   else
      ylabel ('Spectral Density (m^2/Hz)');
   end   
  xl=get (gca,'xlim');
  yl=get (gca,'ylim');
  x=(xl(1)*(xl(2)/xl(1))^.7); 
  y=(yl (1)/(yl(2)/yl(1))^.08);
  text (x,y,txt);

   
   new_spectra=input('Add another plot (y/n) ','s'); 
end
 new_plot=input('New plot (y/n): ','s');
end;
