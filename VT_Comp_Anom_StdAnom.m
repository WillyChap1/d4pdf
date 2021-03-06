%Composite maps 0f vapor transport. 

%First make Anomaly Map...
clear all

% =========================================================================
% ============================ Switches ===================================
% =========================================================================

Ninomonths1 = [12];  %months in first year el nino [Dec]
Ninomonths2 = [1 2];  %months in second year el nino [JAN, FEB] +1 -

%months to test for map.
months1 = [12];  %months in first year ex. [NOV, DEC]  0
months2 = [1 2];  %months in second year ex. [JAN, FEB] +1 -

anomaNINO = [0.5,4];    %min and max anomaly nino
anomaNINA = [0.5,4];    %min and max anomaly nina 

daysinMonth = [31 28 31 30 31 30 31 31 30 31 30 31];

Level = 925;%[1000,925,850,700,600,500,400,300];
for RRa = 1:length(Level)
close all
% =========================================================================
% =========================================================================
% =========================================================================

%==========================================================================
    %======================= Get Mean Response  =====================
%==========================================================================

cd /Users/wchapman/SIO/research/100Amip/d4pdf/Q/ARmasked/

%find the leap years: 
leapstart = 1952;
leapyears = zeros(1,31);
leapyears(1) = leapstart;

for jj = 2:30
    leapyears(jj) = leapyears(jj-1)+4; 
end

findfilesQ = strcat('Q_ARmask',num2str(Level(RRa)),'*.mat');
pp = dir(findfilesQ);

cd /Users/wchapman/SIO/research/100Amip/d4pdf/UV/ARmasked/

findfilesU = strcat('U_ARmask',num2str(Level(RRa)),'*.mat'); 
qq = dir(findfilesU);
findfilesV = strcat('V_ARmask',num2str(Level(RRa)),'*.mat'); 
rr = dir(findfilesV);

cd /Users/wchapman/SIO/research/100Amip/d4pdf/Q/ARmasked/
load(pp(1).name);

datestart = '1951-01-01 00:00';
formatin = 'yyyy-MM-dd HH:mm';
date1 = datetime(datestart,'InputFormat',formatin);
timetot = date1+minutes(time);

Qcat = zeros(size(Qmask,1),size(Qmask,2));
timekeep = [];
timekeepMonths =[];
QkeepAll = [];
QkeepMonth = [];
Qperiod = [];
Uperiod = [];
Vperiod = [];
timekeepPeriod = [];

for hh = 1:(length(pp)-2)
    Qperiodtemp = [];
    Uperiodtemp = [];
    Vperiodtemp = [];
    timePeriod = [];
    
    if isempty(months1) == 0
    disp(strcat('Im in year:',{' '},pp(hh).name(end-11:end-8)))
    fileQ = pp(hh).name;
    disp(fileQ)
    load(fileQ);
    fileU = strcat('/Users/wchapman/SIO/research/100Amip/d4pdf/UV/ARmasked/',qq(hh).name);
    disp(fileU)
    load(fileU);
    fileV = strcat('/Users/wchapman/SIO/research/100Amip/d4pdf/UV/ARmasked/',rr(hh).name);
    load(fileV);
    disp(fileV)
    timetot = date1+minutes(time);
        
    for jr = min(months1):max(months1)    %Dec Nov. 
        [~,month,~]=datevec(timetot);
        indexdo = find(month==jr);
        Qtemp = Qmask(:,:,:,indexdo);
        Utemp = Umask(:,:,:,indexdo);
        Vtemp = Vmask(:,:,:,indexdo);
        
        QkeepMonth = cat(4,QkeepMonth,sum(Qtemp,4));
        %QkeepAll = cat(4,QkeepAll,Qtemp);
        timekeepMonths = cat(1,timekeepMonths,timetot(indexdo(1)));
        timekeep = cat(1,timekeep,timetot(indexdo));
        
        Qperiodtemp = cat(4,Qperiodtemp,Qtemp);
        Uperiodtemp = cat(4,Uperiodtemp,Utemp);
        Vperiodtemp = cat(4,Vperiodtemp,Vtemp);
        
        
        timePeriod = cat(1,timePeriod,timetot(indexdo(1)));
        
        
    end

    end
    
    if isempty(months2)==0
  
    fileQ = pp(hh+1).name;
    disp(fileQ)
    load(fileQ);
    fileU = strcat('/Users/wchapman/SIO/research/100Amip/d4pdf/UV/ARmasked/',qq(hh+1).name);
    disp(fileU)
    load(fileU);
    fileV = strcat('/Users/wchapman/SIO/research/100Amip/d4pdf/UV/ARmasked/',rr(hh+1).name);
    load(fileV);
    disp(fileV)
    
    
    timetot = date1+minutes(time);
        
    for jr = min(months2):max(months2)    %Dec Nov. 
        [~,month,~]=datevec(timetot);
        indexdo = find(month==jr);
        Qtemp = Qmask(:,:,:,indexdo);
        Utemp = Umask(:,:,:,indexdo);
        Vtemp = Vmask(:,:,:,indexdo);
        
        QkeepMonth = cat(4,QkeepMonth,sum(Qtemp,4)); 
        %QkeepAll = cat(4,QkeepAll,Qtemp);
        timekeepMonths = cat(1,timekeepMonths,timetot(indexdo(1)));
        timekeep = cat(1,timekeep,timetot(indexdo));
        timePeriod = cat(1,timePeriod,timetot(indexdo(1)));
        
        Qperiodtemp = cat(4,Qperiodtemp,Qtemp);
        Uperiodtemp = cat(4,Uperiodtemp,Utemp);
        Vperiodtemp = cat(4,Vperiodtemp,Vtemp);
    end
        
    end
    
    Qperiod = cat(4,Qperiod,nansum(Qperiodtemp,4));
    Uperiod = cat(4,Uperiod,nansum(Uperiodtemp,4));
    Vperiod = cat(4,Vperiod,nansum(Vperiodtemp,4));
    timekeepPeriod = cat(1,timekeepPeriod,timePeriod(1));
      
end 



VTperiod = Qperiod.*((Uperiod.^2 + Vperiod.^2).^.5);



VTMean = nanmean(VTperiod,3);
VTMean = squeeze(VTMean);
VTMean = nanmean(VTMean,3);

%%
%==========================================================================
%=========================== Get Nino/Nina Composites =====================
%==========================================================================

%MAKE COMPOSITES FROM YEARS---  across ensembles. 
%el nino years:strong: 57-58,65-66,72-73,87-88,91-92...1997-1998, 1972-1973, 1982-1983, 1997-1998.

if isempty(Ninomonths2)
    Ninomonthsall = Ninomonths1;
else 
    Ninomonthsall = cat(2,Ninomonths1,Ninomonths2);   %combined months for plotting.
end

if isempty(months2)
    monthsall = months1;
else 
    monthsall = cat(2,months1,months2);   %combined months for plotting.
end

Persinmonth = sum(daysinMonth(monthsall))*4;

nino = nino_34_get(1952,2008,Ninomonthsall,1);  %function to get our nino index shit..
avenino = mean(nino);
anomnino = nino-avenino;
figs1= figure;
timeDJF = 1952:1:2008;
plot(timeDJF,anomnino,'Linewidth',1.5)
hold on
ylabel('NINO 3.4 Anomoly')
plot([1952,2009],[anomaNINO(2),anomaNINO(2)],'r--')
plot([1952,2009],[anomaNINO(1),anomaNINO(1)],'r--')

plot([1952,2009],[-anomaNINA(2),-anomaNINA(2)],'b--')
plot([1952,2009],[-anomaNINA(1),-anomaNINA(1)],'b--')

title('Nino 3.4 Anomaly');
xticks(1950:5:2009);
ax = gca;
grid on
ax.FontSize = 25;
pos = get(figs1,'position');
set(figs1,'position',[pos(1)/7 pos(2) pos(3)*4 pos(4)])

Ninoyears = timeDJF(anomnino<=anomaNINO(2) & anomnino>=anomaNINO(1));
Ninayears = timeDJF(anomnino>=-anomaNINA(2) & anomnino<=-anomaNINA(1));
NoNinoyears = timeDJF(anomnino<=anomaNINO(1) & anomnino>=-anomaNINO(1));

hold on 
plot(Ninoyears,anomnino(anomnino<=anomaNINO(2) & anomnino>=anomaNINO(1)),'b*')
plot(Ninayears,anomnino(anomnino>=-anomaNINA(2) & anomnino<=-anomaNINA(1)),'k*')

%==========================================================================
%========================       Nino .       ==============================
%==========================================================================

[yearsIn,~] = datevec(timekeepPeriod);

for ff = 1:length(Ninoyears)
    
    NinoIndex = find(yearsIn==Ninoyears(ff));
    
    VTcatMemNino(:,:,:,ff) = VTperiod(:,:,:,NinoIndex);     
    timeNino(ff) = timekeepPeriod(NinoIndex);
    
end

VTCompositeNino = nanmean(VTcatMemNino,3);
VTCompositeNino = squeeze(VTCompositeNino);
VTCompositeNino = nanmean(VTCompositeNino,3);
VTCompositeNino = squeeze(VTCompositeNino);


%==========================================================================
%========================       Nina .       ==============================
%==========================================================================

for ff = 1:length(Ninayears)
    
    NinaIndex = find(yearsIn==Ninayears(ff));
    
    VTcatMemNina(:,:,:,ff) = VTperiod(:,:,:,NinaIndex); 
    timeNina(ff) = timekeepPeriod(NinaIndex);
    
end

VTCompositeNina = nanmean(VTcatMemNina,3);
VTCompositeNina = squeeze(VTCompositeNina);
VTCompositeNina = nanmean(VTCompositeNina,3);
VTCompositeNina = squeeze(VTCompositeNina);


%==========================================================================
%========================     No  Nino .       ============================
%==========================================================================

for ff = 1:length(NoNinoyears)
    
    NoNinaIndex = find(yearsIn==NoNinoyears(ff));
    
    VTcatMemNoNino(:,:,:,ff) = VTperiod(:,:,:,NoNinaIndex); 
    timeNina(ff) = timekeepPeriod(NoNinaIndex);
    
end

VTCompositeNoNino = nanmean(VTcatMemNoNino,3);
VTCompositeNoNino = squeeze(VTCompositeNoNino);
VTCompositeNoNino = nanmean(VTCompositeNoNino,3);
VTCompositeNoNino = squeeze(VTCompositeNoNino);

%==========================================================================
%==========================================================================
%==========================================================================



%% Make Plotting Variables && plot 
VTmeanYear = nanmean(VTperiod,3);
VTmeanNino = nanmean(VTcatMemNino,3);
VTmeanNina = nanmean(VTcatMemNina,3);
VTmeanNoNino = nanmean(VTcatMemNoNino,3);

VTCompInt = VTperiod - VTmeanYear;
VTCompIntNino = VTcatMemNino - VTmeanNino;
VTCompIntNina = VTcatMemNina - VTmeanNina;
VTCompIntNoNino = VTcatMemNoNino - VTmeanNoNino;

VTVarCompInt = nanvar(VTCompInt,0,3);
VTVarCompInt = (squeeze(VTVarCompInt)).^0.5;
VTStdCompInt = mean(squeeze(VTVarCompInt),3);
    
VTVarCompIntNino = nanvar(VTCompIntNino,0,3);
VTVarCompIntNino = (squeeze(VTVarCompIntNino)).^0.5;
VTStdCompIntNino = mean(squeeze(VTVarCompIntNino),3);
    
VTVarCompIntNina = nanvar(VTCompIntNina,0,3);
VTVarCompIntNina = (squeeze(VTVarCompIntNina)).^0.5;
VTStdCompIntNina = mean(squeeze(VTVarCompIntNina),3);

VTVarCompIntNoNino = nanvar(VTCompIntNoNino,0,3);
VTVarCompIntNoNino = (squeeze(VTVarCompIntNoNino)).^0.5;
VTStdCompIntNoNino = mean(squeeze(VTVarCompIntNoNino),3);

VTplotInt = VTStdCompInt;
VTplotIntNino = VTStdCompIntNino;
VTplotIntNina = VTStdCompIntNina;
VTplotIntNoNino = VTStdCompIntNoNino;

%%
% %%%%%%%%%%%%%%%%%%%%%% .    BOOTSTRAP IT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
disp('Entering BootStrap')
[QlowNINO,QhighNINO] = Bootstrap4D(VTperiod,size(VTcatMemNino,4)/size(VTperiod,4),10000,[1,99]);
[QlowNINA,QhighNINA] = Bootstrap4D(VTperiod,size(VTcatMemNina,4)/size(VTperiod,4),10000,[1,99]);
[QlowNoNINO,QhighNoNINO] = Bootstrap4D(VTperiod,size(VTcatMemNina,4)/size(VTperiod,4),10000,[1,99]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

                                %options: [gap]   %hght [Top Bot]  %Wid  [L R]                          
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
if ~make_it_tight,  clear subplot;  end


load('Q_ARmask300____19520101_19521231.mat');
titMon = 'o';
mons =['J','F','M','A','M','J','J','A','S','O','N','D'];
for yy = 1:length(monthsall)
    
    mon = monthsall(yy);
    
    monC = mons(mon);
    
    titMon = strcat(titMon,monC);      
end 
ax11=figure;
padvar =4;
lonpadf = zeros(padvar,1);
latpadf = zeros(padvar,1);
lonpadb = zeros(padvar,1);
latpadb = zeros(padvar,1);
counter = 1;

for kk = padvar:-1:1 
    lonpadf(counter) = lon(1) - (lon(2)-lon(1))*kk;
    latpadf(counter) = lat(1) - (lat(2)-lat(1))*kk;
    counter = counter+1; 
end 

for kk = 1:padvar 
    lonpadb(kk) = lon(end) + (lon(2)-lon(1))*kk;
    latpadb(kk) = lat(end) + (lat(2)-lat(1))*kk;
end 

lon = cat(1,lon,lonpadb);
lat = cat(1,lat,latpadb);

lon = cat(1,lonpadf,lon);
lat = cat(1,latpadf,lat);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find color scale%
fmax = cat(3,VTplotIntNino,VTplotIntNina,VTplotIntNoNino);
fmax(fmax==Inf)=NaN;
fmax(fmax==-Inf)=NaN;
cv12 = max(max(max(abs(fmax))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


aa2=subplot(3,4,2);
padvar =4;
padCOMPmean = padarray(VTCompositeNino,[padvar,padvar],nan,'both');
VTplot = padCOMPmean;
VTplotInt =  padarray((VTplotIntNino),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,VTplot','Color',[1 1 1],'LineWidth',3);
hold on
m_contourf(lon,lat,VTplotInt',24);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Blues',24);
colormap(aa2,RuBU);
title(strcat('Nino',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR VT'))
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
caxis([0 cv12])
ylabel(h1,'Std Dev of VT (kg/kg*m/s)');

aa3=subplot(3,4,3);
padvar =4;
padCOMPmean = padarray(VTCompositeNina,[padvar,padvar],nan,'both');
VTplot = padCOMPmean;
VTplotInt =  padarray((VTplotIntNina),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,VTplot','Color',[1,1,1],'LineWidth',3);
hold on
m_contourf(lon,lat,VTplotInt',24);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Blues',24);
colormap(aa3,RuBU);
title(strcat('Nina',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR  VT'))
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
caxis([0 cv12])
ylabel(h1,'Std Dev of VT (kg/kg*m/s)');

aa4=subplot(3,4,4);
padvar =4;
padCOMPmean = padarray(VTCompositeNoNino,[padvar,padvar],nan,'both');
VTplot = padCOMPmean;
VTplotInt =  padarray((VTplotIntNoNino),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,VTplot','Color',[1,1,1],'LineWidth',3);
hold on
m_contourf(lon,lat,VTplotInt',24);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Blues',24);
colormap(aa4,RuBU);
title(strcat('Nuetral',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR VT'))
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
caxis([0 cv12])
ylabel(h1,'Std Dev of VT (kg/kg*m/s)');



% fig panel 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Plotting Variables && plot 
VTmeanYear = nanmean(VTperiod,3);
VTmeanNino = nanmean(VTcatMemNino,3);
VTmeanNina = nanmean(VTcatMemNina,3);
VTmeanNoNino = nanmean(VTcatMemNoNino,3);

VTCompInt = VTperiod - VTmeanYear;
VTCompIntNino = VTcatMemNino - VTmeanNino;
VTCompIntNina = VTcatMemNina - VTmeanNina;
VTCompIntNoNino = VTcatMemNoNino - VTmeanNoNino;

VTVarCompInt = nanvar(VTCompInt,0,3);
VTVarCompInt = (squeeze(VTVarCompInt)).^0.5;
VTStdCompInt = mean(squeeze(VTVarCompInt),3);
    
VTVarCompIntNino = nanvar(VTCompIntNino,0,3);
VTVarCompIntNino = (squeeze(VTVarCompIntNino)).^0.5;
VTStdCompIntNino = mean(squeeze(VTVarCompIntNino),3);
    
VTVarCompIntNina = nanvar(VTCompIntNina,0,3);
VTVarCompIntNina = (squeeze(VTVarCompIntNina)).^0.5;
VTStdCompIntNina = mean(squeeze(VTVarCompIntNina),3);

VTVarCompIntNoNino = nanvar(VTCompIntNoNino,0,3);
VTVarCompIntNoNino = (squeeze(VTVarCompIntNoNino)).^0.5;
VTStdCompIntNoNino = mean(squeeze(VTVarCompIntNoNino),3);

VTplotInt = VTStdCompInt;
VTplotIntNino = VTStdCompIntNino;
VTplotIntNina = VTStdCompIntNina;
VTplotIntNoNino = VTStdCompIntNoNino;

load('Q_ARmask300____19520101_19521231.mat');
titMon = 'o';
mons =['J','F','M','A','M','J','J','A','S','O','N','D'];
for yy = 1:length(monthsall)
    
    mon = monthsall(yy);
    
    monC = mons(mon);
    
    titMon = strcat(titMon,monC);      
end 
aa1=subplot(3,4,5);
padvar =4;
padCOMPmean = padarray(VTMean,[padvar,padvar],nan,'both');
lonpadf = zeros(padvar,1);
latpadf = zeros(padvar,1);
lonpadb = zeros(padvar,1);
latpadb = zeros(padvar,1);
counter = 1;

for kk = padvar:-1:1 
    lonpadf(counter) = lon(1) - (lon(2)-lon(1))*kk;
    latpadf(counter) = lat(1) - (lat(2)-lat(1))*kk;
    counter = counter+1; 
end 

for kk = 1:padvar 
    lonpadb(kk) = lon(end) + (lon(2)-lon(1))*kk;
    latpadb(kk) = lat(end) + (lat(2)-lat(1))*kk;
end 

lon = cat(1,lon,lonpadb);
lat = cat(1,lat,latpadb);

lon = cat(1,lonpadf,lon);
lat = cat(1,latpadf,lat);

VTplot = padCOMPmean;
VTdivide = VTplotInt;
VTplotInt =  padarray((VTplotInt),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,VTplot','Color',[1 1 1],'LineWidth',3);
hold on
m_contourf(lon,lat,VTplotInt',24);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Blues',24);
colormap(aa1,RuBU);
title(strcat(titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR VT'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([0 cv12])
ylabel(h1,'Std Dev of Q (kg/kg)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find color scale%

fmax = cat(3,VTCompositeNino-VTMean,VTCompositeNina-VTMean,VTCompositeNoNino-VTMean);
fmax(fmax==Inf)=NaN;
fmax(fmax==-Inf)=NaN;
cv1 = max(max(max(abs(fmax))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aa2=subplot(3,4,6);
padvar =4;
padCOMPmean = padarray(VTCompositeNino-VTMean,[padvar,padvar],nan,'both');

moreNINO = padCOMPmean(5:end-4,5:end-4)>QhighNINO;
lessNINO = padCOMPmean(5:end-4,5:end-4)<QlowNINO;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

VTplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
m_contourf(lon,lat,VTplot','Color',[1 1 1],'LineWidth',3);
hold on
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',3,'markerfacecolor','w');
RuBU = cbrewer('div','RdBu',24);
colormap(aa2,RuBU);
title(strcat('Nino',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR VT'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv1 cv1])
ylabel(h1,'VT (kg/kg * m/s)');

aa3=subplot(3,4,7);
padvar =4;
padCOMPmean = padarray(VTCompositeNina-VTMean,[padvar,padvar],nan,'both');

moreNINO = padCOMPmean(5:end-4,5:end-4)>QhighNINA;
lessNINO = padCOMPmean(5:end-4,5:end-4)<QlowNINA;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

VTplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
m_contourf(lon,lat,VTplot','Color',[1,1,1],'LineWidth',3);
hold on
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',3,'markerfacecolor','w');
RuBU = cbrewer('div','RdBu',24);
colormap(aa3,RuBU);
title(strcat('Nina',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR VT'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv1 cv1])
ylabel(h1,'VT (kg/kg * m/s)');

aa4=subplot(3,4,8);
padvar =4;
padCOMPmean = padarray(VTCompositeNoNino-VTMean,[padvar,padvar],nan,'both');

moreNINO = padCOMPmean(5:end-4,5:end-4)>QhighNoNINO;
lessNINO = padCOMPmean(5:end-4,5:end-4)<QlowNoNINO;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);


VTplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
m_contourf(lon,lat,VTplot','Color',[1,1,1],'LineWidth',3);
hold on
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',3,'markerfacecolor','w');
RuBU = cbrewer('div','RdBu',24);
colormap(aa4,RuBU);
title(strcat('Nuetral',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR VT'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv1 cv1])
ylabel(h1,'VT (kg/kg * m/s)');


% Fig Panel 3
% Make Plotting Variables && plot 
VTmeanYear = nanmean(VTperiod,3);
VTmeanNino = nanmean(VTcatMemNino,3);
VTmeanNina = nanmean(VTcatMemNina,3);
VTmeanNoNino = nanmean(VTcatMemNoNino,3);

VTCompInt = VTperiod - VTmeanYear;
VTCompIntNino = VTcatMemNino - VTmeanNino;
VTCompIntNina = VTcatMemNina - VTmeanNina;
VTCompIntNoNino = VTcatMemNoNino - VTmeanNoNino;

VTVarCompInt = nanvar(VTCompInt,0,3);
VTVarCompInt = (squeeze(VTVarCompInt)).^0.5;
VTStdCompInt = mean(squeeze(VTVarCompInt),3);
    
VTVarCompIntNino = nanvar(VTCompIntNino,0,3);
VTVarCompIntNino = (squeeze(VTVarCompIntNino)).^0.5;
VTStdCompIntNino = mean(squeeze(VTVarCompIntNino),3);
    
VTVarCompIntNina = nanvar(VTCompIntNina,0,3);
VTVarCompIntNina = (squeeze(VTVarCompIntNina)).^0.5;
VTStdCompIntNina = mean(squeeze(VTVarCompIntNina),3);

VTVarCompIntNoNino = nanvar(VTCompIntNoNino,0,3);
VTVarCompIntNoNino = (squeeze(VTVarCompIntNoNino)).^0.5;
VTStdCompIntNoNino = mean(squeeze(VTVarCompIntNoNino),3);

VTplotInt = VTStdCompInt;
VTplotIntNino = VTStdCompIntNino;
VTplotIntNina = VTStdCompIntNina;
VTplotIntNoNino = VTStdCompIntNoNino;


load('Q_ARmask300____19520101_19521231.mat');
titMon = 'o';
mons =['J','F','M','A','M','J','J','A','S','O','N','D'];
for yy = 1:length(monthsall)
    
    mon = monthsall(yy);
    
    monC = mons(mon);
    
    titMon = strcat(titMon,monC);      
end

padvar =4;
lonpadf = zeros(padvar,1);
latpadf = zeros(padvar,1);
lonpadb = zeros(padvar,1);
latpadb = zeros(padvar,1);
counter = 1;

for kk = padvar:-1:1 
    lonpadf(counter) = lon(1) - (lon(2)-lon(1))*kk;
    latpadf(counter) = lat(1) - (lat(2)-lat(1))*kk;
    counter = counter+1; 
end 

for kk = 1:padvar 
    lonpadb(kk) = lon(end) + (lon(2)-lon(1))*kk;
    latpadb(kk) = lat(end) + (lat(2)-lat(1))*kk;
end 

lon = cat(1,lon,lonpadb);
lat = cat(1,lat,latpadb);

lon = cat(1,lonpadf,lon);
lat = cat(1,latpadf,lat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find color scale%

fmax = cat(3,(VTCompositeNino-VTMean)./VTdivide,(VTCompositeNina-VTMean)./VTdivide,(VTCompositeNoNino-VTMean)./VTdivide);
fmax(fmax==Inf)=NaN;
fmax(fmax==-Inf)=NaN;
cv1 = max(max(max(abs(fmax))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aa2=subplot(3,4,10);
padvar =4;
padCOMPmean = padarray((VTCompositeNino-VTMean)./VTdivide,[padvar,padvar],nan,'both');
moreNINO = (VTCompositeNino-VTMean)>QhighNINO;
lessNINO = (VTCompositeNino-VTMean)<QlowNINO;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

VTplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,VTplot','Color',[1 1 1],'LineWidth',3);
hold on
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',2,'markerfacecolor','w');
RuBU = cbrewer('div','BrBG',24);
colormap(aa2,RuBU);
title(strcat('Nino',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR VT'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv1 cv1])

aa3=subplot(3,4,11);
padvar =4;
padCOMPmean = padarray((VTCompositeNina-VTMean)./VTdivide,[padvar,padvar],nan,'both');

moreNINO = (VTCompositeNina-VTMean)>QhighNINA;
lessNINO = (VTCompositeNina-VTMean)<QlowNINA;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

VTplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,VTplot','Color',[1,1,1],'LineWidth',3);
hold on
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',2,'markerfacecolor','w');
RuBU = cbrewer('div','BrBG',24);
colormap(aa3,RuBU);
title(strcat('Nina',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR VT'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv1 cv1])

aa4=subplot(3,4,12);
padvar =4;
padCOMPmean = padarray((VTCompositeNoNino-VTMean)./VTdivide,[padvar,padvar],nan,'both');

moreNINO = (VTCompositeNoNino-VTMean)>QhighNoNINO;
lessNINO = (VTCompositeNoNino-VTMean)<QlowNoNINO;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);


VTplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,VTplot','Color',[1,1,1],'LineWidth',3);
hold on
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',2,'markerfacecolor','w');
RuBU = cbrewer('div','BrBG',24);
colormap(aa4,RuBU);
title(strcat('Nuetral',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR VT'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv1 cv1])


%set position
set(ax11,'position',[ 229          19        2045        1326])

%annotate each row. 

% %1
assx = annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Mean', ...
    'EdgeColor', 'none', ...
    'Position', [0.59    0.9000    1.0000    0.1000],'FontSize',25,'FontWeight','bold');

annotation('line',[0.265 .97],[0.978 0.978],'color','k');


% %2
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Mean', ...
    'EdgeColor', 'none', ...
    'Position', [0.11    0.574000    1.0000    0.1000],'FontSize',25,'FontWeight','bold');

annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Anomaly', ...
    'EdgeColor', 'none', ...
    'Position', [0.58    0.574000    1.0000    0.1000],'FontSize',25,'FontWeight','bold');

annotation('line',[0.03 .258],[0.65 0.65],'color','k');
annotation('line',[0.265 .97],[0.65 0.65],'color','k');
% 
% 
% %3

annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Standard Anomaly', ...
    'EdgeColor', 'none', ...
    'Position', [0.55    0.255000    1.0000    0.1000],'FontSize',25,'FontWeight','bold');

annotation('line',[0.265 .97],[0.33 0.33],'color','k');

%vert line 


annotation('line',[0.2615 0.2615],[0.978 0.05],'color','k');



figsave=strcat('/Users/wchapman/Desktop/QUVfigs/',num2str(titMon(2:end)),'_',num2str(Level(RRa)),'_VTall.jpg');
saveas(ax11,figsave);
end







