%Composite maps 

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

Level = [925,850,700,600,500,400,300];
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


QMean = nanmean(Qperiod,3);
QMean = squeeze(QMean);
QMean = nanmean(QMean,3);

VMean = nanmean(Vperiod,3);
VMean = squeeze(VMean);
VMean = nanmean(VMean,3);

UMean = nanmean(Uperiod,3);
UMean = squeeze(UMean);
UMean = nanmean(UMean,3);
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
    
    QcatMemNino(:,:,:,ff) = Qperiod(:,:,:,NinoIndex); 
    UcatMemNino(:,:,:,ff) = Uperiod(:,:,:,NinoIndex);
    VcatMemNino(:,:,:,ff) = Vperiod(:,:,:,NinoIndex);
    timeNino(ff) = timekeepPeriod(NinoIndex);
    
end

QCompositeNino = nanmean(QcatMemNino,3);
QCompositeNino = squeeze(QCompositeNino);
QCompositeNino = nanmean(QCompositeNino,3);
QCompositeNino = squeeze(QCompositeNino);

UCompositeNino = nanmean(UcatMemNino,3);
UCompositeNino = squeeze(UCompositeNino);
UCompositeNino = nanmean(UCompositeNino,3);
UCompositeNino = squeeze(UCompositeNino);

VCompositeNino = nanmean(VcatMemNino,3);
VCompositeNino = squeeze(VCompositeNino);
VCompositeNino = nanmean(VCompositeNino,3);
VCompositeNino = squeeze(VCompositeNino);

%==========================================================================
%========================       Nina .       ==============================
%==========================================================================

for ff = 1:length(Ninayears)
    
    NinaIndex = find(yearsIn==Ninayears(ff));
    
    QcatMemNina(:,:,:,ff) = Qperiod(:,:,:,NinaIndex); 
    UcatMemNina(:,:,:,ff) = Uperiod(:,:,:,NinaIndex);
    VcatMemNina(:,:,:,ff) = Vperiod(:,:,:,NinaIndex);
    timeNina(ff) = timekeepPeriod(NinaIndex);
    
end

QCompositeNina = nanmean(QcatMemNina,3);
QCompositeNina = squeeze(QCompositeNina);
QCompositeNina = nanmean(QCompositeNina,3);
QCompositeNina = squeeze(QCompositeNina);

UCompositeNina = nanmean(UcatMemNina,3);
UCompositeNina = squeeze(UCompositeNina);
UCompositeNina = nanmean(UCompositeNina,3);
UCompositeNina = squeeze(UCompositeNina);

VCompositeNina = nanmean(VcatMemNina,3);
VCompositeNina = squeeze(VCompositeNina);
VCompositeNina = nanmean(VCompositeNina,3);
VCompositeNina = squeeze(VCompositeNina);

%==========================================================================
%========================     No  Nino .       ============================
%==========================================================================

for ff = 1:length(NoNinoyears)
    
    NoNinaIndex = find(yearsIn==NoNinoyears(ff));
    
    QcatMemNoNino(:,:,:,ff) = Qperiod(:,:,:,NoNinaIndex); 
    UcatMemNoNino(:,:,:,ff) = Uperiod(:,:,:,NoNinaIndex);
    VcatMemNoNino(:,:,:,ff) = Vperiod(:,:,:,NoNinaIndex);
    timeNina(ff) = timekeepPeriod(NoNinaIndex);
    
end
QCompositeNoNino = nanmean(QcatMemNoNino,3);
QCompositeNoNino = squeeze(QCompositeNoNino);
QCompositeNoNino = nanmean(QCompositeNoNino,3);
QCompositeNoNino = squeeze(QCompositeNoNino);

UCompositeNoNino = nanmean(UcatMemNoNino,3);
UCompositeNoNino = squeeze(UCompositeNoNino);
UCompositeNoNino = nanmean(UCompositeNoNino,3);
UCompositeNoNino = squeeze(UCompositeNoNino);

VCompositeNoNino = nanmean(VcatMemNoNino,3);
VCompositeNoNino = squeeze(VCompositeNoNino);
VCompositeNoNino = nanmean(VCompositeNoNino,3);
VCompositeNoNino = squeeze(VCompositeNoNino);

%==========================================================================
%==========================================================================
%==========================================================================



%% Make Plotting Variables && plot 
QmeanYear = nanmean(Qperiod,3);
QmeanNino = nanmean(QcatMemNino,3);
QmeanNina = nanmean(QcatMemNina,3);
QmeanNoNino = nanmean(QcatMemNoNino,3);

QCompInt = Qperiod - QmeanYear;
QCompIntNino = QcatMemNino - QmeanNino;
QCompIntNina = QcatMemNina - QmeanNina;
QCompIntNoNino = QcatMemNoNino - QmeanNoNino;

QVarCompInt = nanvar(QCompInt,0,3);
QVarCompInt = (squeeze(QVarCompInt)).^0.5;
QStdCompInt = mean(squeeze(QVarCompInt),3);
    
QVarCompIntNino = nanvar(QCompIntNino,0,3);
QVarCompIntNino = (squeeze(QVarCompIntNino)).^0.5;
QStdCompIntNino = mean(squeeze(QVarCompIntNino),3);
    
QVarCompIntNina = nanvar(QCompIntNina,0,3);
QVarCompIntNina = (squeeze(QVarCompIntNina)).^0.5;
QStdCompIntNina = mean(squeeze(QVarCompIntNina),3);

QVarCompIntNoNino = nanvar(QCompIntNoNino,0,3);
QVarCompIntNoNino = (squeeze(QVarCompIntNoNino)).^0.5;
QStdCompIntNoNino = mean(squeeze(QVarCompIntNoNino),3);

QplotInt = QStdCompInt;
QplotIntNino = QStdCompIntNino;
QplotIntNina = QStdCompIntNina;
QplotIntNoNino = QStdCompIntNoNino;


UmeanYear = nanmean(Uperiod,3);
UmeanNino = nanmean(UcatMemNino,3);
UmeanNina = nanmean(UcatMemNina,3);
UmeanNoNino = nanmean(UcatMemNoNino,3);


UCompInt = Uperiod - UmeanYear;
UCompIntNino = UcatMemNino - UmeanNino;
UCompIntNina = UcatMemNina - UmeanNina;
UCompIntNoNino = UcatMemNoNino - UmeanNoNino;


UVarCompInt = nanvar(UCompInt,0,3);
UVarCompInt = (squeeze(UVarCompInt)).^0.5;
UStdCompInt = mean(squeeze(UVarCompInt),3);
    
UVarCompIntNino = nanvar(UCompIntNino,0,3);
UVarCompIntNino = (squeeze(UVarCompIntNino)).^0.5;
UStdCompIntNino = mean(squeeze(UVarCompIntNino),3);
    
UVarCompIntNina = nanvar(UCompIntNina,0,3);
UVarCompIntNina = (squeeze(UVarCompIntNina)).^0.5;
UStdCompIntNina = mean(squeeze(UVarCompIntNina),3);

UVarCompIntNoNino = nanvar(UCompIntNoNino,0,3);
UVarCompIntNoNino = (squeeze(UVarCompIntNoNino)).^0.5;
UStdCompIntNoNino = mean(squeeze(UVarCompIntNoNino),3);

UplotInt = UStdCompInt;
UplotIntNino = UStdCompIntNino;
UplotIntNina = UStdCompIntNina;
UplotIntNoNino = UStdCompIntNoNino;


VmeanYear = nanmean(Vperiod,3);
VmeanNino = nanmean(VcatMemNino,3);
VmeanNina = nanmean(VcatMemNina,3);
VmeanNoNino = nanmean(VcatMemNoNino,3);


VCompInt = Vperiod - VmeanYear;
VCompIntNino = VcatMemNino - VmeanNino;
VCompIntNina = VcatMemNina - VmeanNina;
VCompIntNoNino = VcatMemNoNino - VmeanNoNino;


VVarCompInt = nanvar(VCompInt,0,3);
VVarCompInt = (squeeze(VVarCompInt)).^0.5;
VStdCompInt = mean(squeeze(VVarCompInt),3);
    
VVarCompIntNino = nanvar(VCompIntNino,0,3);
VVarCompIntNino = (squeeze(VVarCompIntNino)).^0.5;
VStdCompIntNino = mean(squeeze(VVarCompIntNino),3);
    
VVarCompIntNina = nanvar(VCompIntNina,0,3);
VVarCompIntNina = (squeeze(VVarCompIntNina)).^0.5;
VStdCompIntNina = mean(squeeze(VVarCompIntNina),3);

VVarCompIntNoNino = nanvar(VCompIntNoNino,0,3);
VVarCompIntNoNino = (squeeze(VVarCompIntNoNino)).^0.5;
VStdCompIntNoNino = mean(squeeze(VVarCompIntNoNino),3);

VplotInt = VStdCompInt;
VplotIntNino = VStdCompIntNino;
VplotIntNina = VStdCompIntNina;
VplotIntNoNino = VStdCompIntNoNino;

%%
% %%%%%%%%%%%%%%%%%%%%%% .    BOOTSTRAP IT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
disp('Entering BootStrap')
[QlowNINO,QhighNINO] = Bootstrap4D(Qperiod,size(QcatMemNino,4)/size(Qperiod,4),10000,[1,99]);
[QlowNINA,QhighNINA] = Bootstrap4D(Qperiod,size(QcatMemNina,4)/size(Qperiod,4),10000,[1,99]);
[QlowNoNINO,QhighNoNINO] = Bootstrap4D(Qperiod,size(QcatMemNina,4)/size(Qperiod,4),10000,[1,99]);

[UlowNINO,UhighNINO] = Bootstrap4D(Uperiod,size(UcatMemNino,4)/size(Uperiod,4),10000,[1,99]);
[UlowNINA,UhighNINA] = Bootstrap4D(Uperiod,size(UcatMemNina,4)/size(Uperiod,4),10000,[1,99]);
[UlowNoNINO,UhighNoNINO] = Bootstrap4D(Uperiod,size(UcatMemNina,4)/size(Uperiod,4),10000,[1,99]);

[VlowNINO,VhighNINO] = Bootstrap4D(Vperiod,size(VcatMemNino,4)/size(Vperiod,4),10000,[1,99]);
[VlowNINA,VhighNINA] = Bootstrap4D(Vperiod,size(VcatMemNina,4)/size(Vperiod,4),10000,[1,99]);
[VlowNoNINO,VhighNoNINO] = Bootstrap4D(Vperiod,size(VcatMemNina,4)/size(Vperiod,4),10000,[1,99]);
disp('Exiting BootStrap')
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
aa1=subplot(3,4,1);
padvar =4;
padCOMPmean = padarray(QMean,[padvar,padvar],nan,'both');
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

Qplot = padCOMPmean;
QplotInt =  padarray((QplotInt),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,Qplot','Color',[1 1 1],'LineWidth',3);
hold on
m_contourf(lon,lat,QplotInt',24);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Oranges',24);
colormap(aa1,RuBU);
title(strcat(titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR Q'))
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
cv1 =max(max(QplotInt));
caxis([0 cv1])
ylabel(h1,'Std Dev of Q (kg/kg)');


aa2=subplot(3,4,2);
padvar =4;
padCOMPmean = padarray(QCompositeNino,[padvar,padvar],nan,'both');
Qplot = padCOMPmean;
QplotInt =  padarray((QplotIntNino),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,Qplot','Color',[1 1 1],'LineWidth',3);
hold on
m_contourf(lon,lat,QplotInt',24);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Oranges',24);
colormap(aa2,RuBU);
title(strcat('Nino',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR Q'))
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
caxis([0 cv1])
ylabel(h1,'Std Dev of Q (kg/kg)');

aa3=subplot(3,4,3);
padvar =4;
padCOMPmean = padarray(QCompositeNina,[padvar,padvar],nan,'both');
Qplot = padCOMPmean;
QplotInt =  padarray((QplotIntNina),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,Qplot','Color',[1,1,1],'LineWidth',3);
hold on
m_contourf(lon,lat,QplotInt',24);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Oranges',24);
colormap(aa3,RuBU);
title(strcat('Nina',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR Q'))
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
caxis([0 cv1])
ylabel(h1,'Std Dev of Q (kg/kg)');

aa4=subplot(3,4,4);
padvar =4;
padCOMPmean = padarray(QCompositeNoNino,[padvar,padvar],nan,'both');
Qplot = padCOMPmean;
QplotInt =  padarray((QplotIntNoNino),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,Qplot','Color',[1,1,1],'LineWidth',3);
hold on
m_contourf(lon,lat,QplotInt',24);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Oranges',24);
colormap(aa4,RuBU);
title(strcat('Nuetral',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR Q'))
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
caxis([0 cv1])
ylabel(h1,'Std Dev of Q (kg/kg)');

aa5=subplot(3,4,5);
padvar =4;
padCOMPmean = padarray(UMean,[padvar,padvar],nan,'both');
Uplot = padCOMPmean;
UplotInt =  padarray((UplotInt),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,Uplot','Color',[1,1,1],'LineWidth',3);
hold on
m_contourf(lon,lat,UplotInt',24);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Blues',24);
colormap(aa5,RuBU);
title(strcat(titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR U'))
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
cv2 =max(max(UplotInt));
caxis([0 cv2])
ylabel(h1,'Std Dev of U wind (m/s)');

aa6=subplot(3,4,6);
padvar =4;
padCOMPmean = padarray(UCompositeNino,[padvar,padvar],nan,'both');
Uplot = padCOMPmean;
UplotInt =  padarray((UplotIntNino),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,Uplot','Color',[1,1,1],'LineWidth',3);
hold on
m_contourf(lon,lat,UplotInt',24);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Blues',24);
colormap(aa6,RuBU);
title(strcat('Nino',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR U'))
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
caxis([0 cv2])
ylabel(h1,'Std Dev of U wind (m/s)');

aa7 = subplot(3,4,7);
padvar =4;
padCOMPmean = padarray(UCompositeNina,[padvar,padvar],nan,'both');
Uplot = padCOMPmean;
UplotInt =  padarray((UplotIntNina),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,Uplot','Color',[1,1,1],'LineWidth',3);
hold on
m_contourf(lon,lat,UplotInt',24);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Blues',24);
colormap(aa7,RuBU);
title(strcat('Nina',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR U'))
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
caxis([0 cv2])
ylabel(h1,'Std Dev of U wind (m/s)');


aa8=subplot(3,4,8);
padvar =4;
padCOMPmean = padarray(UCompositeNoNino,[padvar,padvar],nan,'both');
Uplot = padCOMPmean;
UplotInt =  padarray((UplotIntNoNino),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,Uplot','Color',[1,1,1],'LineWidth',3);
hold on
m_contourf(lon,lat,UplotInt',24);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Blues',24);
colormap(aa8,RuBU);
title(strcat('Nuetral',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR U'))
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
caxis([0 cv2])
ylabel(h1,'Std Dev of U wind (m/s)');

aa9=subplot(3,4,9);
padvar =4;
padCOMPmean = padarray(VMean,[padvar,padvar],nan,'both');
Vplot = padCOMPmean;
VplotInt =  padarray((VplotInt),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,Vplot','Color',[1,1,1],'LineWidth',3);
hold on
m_contourf(lon,lat,VplotInt',24);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Greens',24);
colormap(aa9,RuBU);
title(strcat(titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR V'))
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
cv2 =max(max(VplotInt));
caxis([0 cv2])
ylabel(h1,'Std Dev of V wind (m/s)');

aa10=subplot(3,4,10);
padvar =4;
padCOMPmean = padarray(VCompositeNino,[padvar,padvar],nan,'both');
Vplot = padCOMPmean;
VplotInt =  padarray((VplotIntNino),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,Vplot','Color',[1,1,1],'LineWidth',3);
hold on
m_contourf(lon,lat,VplotInt',24);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Greens',24);
colormap(aa10,RuBU);
title(strcat('Nino',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR V'))
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
caxis([0 cv2])
ylabel(h1,'Std Dev of V wind (m/s)');

aa11 = subplot(3,4,11);
padvar =4;
padCOMPmean = padarray(VCompositeNina,[padvar,padvar],nan,'both');
Vplot = padCOMPmean;
VplotInt =  padarray((VplotIntNina),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,Vplot','Color',[1,1,1],'LineWidth',3);
hold on
m_contourf(lon,lat,VplotInt',24);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Greens',24);
colormap(aa11,RuBU);
title(strcat('Nina',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR V'))
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
caxis([0 cv2])
ylabel(h1,'Std Dev of V wind (m/s)');

aa12=subplot(3,4,12);
padvar =4;
padCOMPmean = padarray(VCompositeNoNino,[padvar,padvar],nan,'both');
Vplot = padCOMPmean;
VplotInt =  padarray((VplotIntNoNino),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,Vplot','Color',[1,1,1],'LineWidth',3);
hold on
m_contourf(lon,lat,VplotInt',24);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Greens',24);
colormap(aa12,RuBU);
title(strcat('Nuetral',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR V'))
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
caxis([0 cv2])
ylabel(h1,'Std Dev of V wind (m/s)');


set(ax11,'position',[ 203          92        1856        1253])


figsave=strcat('/Users/wchapman/Desktop/QUVfigs/',num2str(titMon(2:end)),'_',num2str(Level(RRa)),'_All.jpg');
saveas(ax11,figsave);
close all


%% fig panel 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make Plotting Variables && plot 
QmeanYear = nanmean(Qperiod,3);
QmeanNino = nanmean(QcatMemNino,3);
QmeanNina = nanmean(QcatMemNina,3);
QmeanNoNino = nanmean(QcatMemNoNino,3);

QCompInt = Qperiod - QmeanYear;
QCompIntNino = QcatMemNino - QmeanNino;
QCompIntNina = QcatMemNina - QmeanNina;
QCompIntNoNino = QcatMemNoNino - QmeanNoNino;

QVarCompInt = nanvar(QCompInt,0,3);
QVarCompInt = (squeeze(QVarCompInt)).^0.5;
QStdCompInt = mean(squeeze(QVarCompInt),3);
    
QVarCompIntNino = nanvar(QCompIntNino,0,3);
QVarCompIntNino = (squeeze(QVarCompIntNino)).^0.5;
QStdCompIntNino = mean(squeeze(QVarCompIntNino),3);
    
QVarCompIntNina = nanvar(QCompIntNina,0,3);
QVarCompIntNina = (squeeze(QVarCompIntNina)).^0.5;
QStdCompIntNina = mean(squeeze(QVarCompIntNina),3);

QVarCompIntNoNino = nanvar(QCompIntNoNino,0,3);
QVarCompIntNoNino = (squeeze(QVarCompIntNoNino)).^0.5;
QStdCompIntNoNino = mean(squeeze(QVarCompIntNoNino),3);

QplotInt = QStdCompInt;
QplotIntNino = QStdCompIntNino;
QplotIntNina = QStdCompIntNina;
QplotIntNoNino = QStdCompIntNoNino;


UmeanYear = nanmean(Uperiod,3);
UmeanNino = nanmean(UcatMemNino,3);
UmeanNina = nanmean(UcatMemNina,3);
UmeanNoNino = nanmean(UcatMemNoNino,3);


UCompInt = Uperiod - UmeanYear;
UCompIntNino = UcatMemNino - UmeanNino;
UCompIntNina = UcatMemNina - UmeanNina;
UCompIntNoNino = UcatMemNoNino - UmeanNoNino;


UVarCompInt = nanvar(UCompInt,0,3);
UVarCompInt = (squeeze(UVarCompInt)).^0.5;
UStdCompInt = mean(squeeze(UVarCompInt),3);
    
UVarCompIntNino = nanvar(UCompIntNino,0,3);
UVarCompIntNino = (squeeze(UVarCompIntNino)).^0.5;
UStdCompIntNino = mean(squeeze(UVarCompIntNino),3);
    
UVarCompIntNina = nanvar(UCompIntNina,0,3);
UVarCompIntNina = (squeeze(UVarCompIntNina)).^0.5;
UStdCompIntNina = mean(squeeze(UVarCompIntNina),3);

UVarCompIntNoNino = nanvar(UCompIntNoNino,0,3);
UVarCompIntNoNino = (squeeze(UVarCompIntNoNino)).^0.5;
UStdCompIntNoNino = mean(squeeze(UVarCompIntNoNino),3);

UplotInt = UStdCompInt;
UplotIntNino = UStdCompIntNino;
UplotIntNina = UStdCompIntNina;
UplotIntNoNino = UStdCompIntNoNino;


VmeanYear = nanmean(Vperiod,3);
VmeanNino = nanmean(VcatMemNino,3);
VmeanNina = nanmean(VcatMemNina,3);
VmeanNoNino = nanmean(VcatMemNoNino,3);


VCompInt = Vperiod - VmeanYear;
VCompIntNino = VcatMemNino - VmeanNino;
VCompIntNina = VcatMemNina - VmeanNina;
VCompIntNoNino = VcatMemNoNino - VmeanNoNino;


VVarCompInt = nanvar(VCompInt,0,3);
VVarCompInt = (squeeze(VVarCompInt)).^0.5;
VStdCompInt = mean(squeeze(VVarCompInt),3);
    
VVarCompIntNino = nanvar(VCompIntNino,0,3);
VVarCompIntNino = (squeeze(VVarCompIntNino)).^0.5;
VStdCompIntNino = mean(squeeze(VVarCompIntNino),3);
    
VVarCompIntNina = nanvar(VCompIntNina,0,3);
VVarCompIntNina = (squeeze(VVarCompIntNina)).^0.5;
VStdCompIntNina = mean(squeeze(VVarCompIntNina),3);

VVarCompIntNoNino = nanvar(VCompIntNoNino,0,3);
VVarCompIntNoNino = (squeeze(VVarCompIntNoNino)).^0.5;
VStdCompIntNoNino = mean(squeeze(VVarCompIntNoNino),3);

VplotInt = VStdCompInt;
VplotIntNino = VStdCompIntNino;
VplotIntNina = VStdCompIntNina;
VplotIntNoNino = VStdCompIntNoNino;


load('Q_ARmask300____19520101_19521231.mat');
titMon = 'o';
mons =['J','F','M','A','M','J','J','A','S','O','N','D'];
for yy = 1:length(monthsall)
    
    mon = monthsall(yy);
    
    monC = mons(mon);
    
    titMon = strcat(titMon,monC);      
end 
ax11=figure;
aa1=subplot(3,4,1);
padvar =4;
padCOMPmean = padarray(QMean,[padvar,padvar],nan,'both');
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

Qplot = padCOMPmean;
QplotInt =  padarray((QplotInt),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,Qplot','Color',[1 1 1],'LineWidth',3);
hold on
m_contourf(lon,lat,QplotInt',24);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Oranges',24);
colormap(aa1,RuBU);
title(strcat(titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR Q'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
cv1 =max(max(QplotInt));
caxis([0 cv1])
ylabel(h1,'Std Dev of Q (kg/kg)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find color scale%

fmax = cat(3,QCompositeNino-QMean,QCompositeNina-QMean,QCompositeNoNino-QMean);
fmax(fmax==Inf)=NaN;
fmax(fmax==-Inf)=NaN;
cv1 = max(max(max(abs(fmax))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aa2=subplot(3,4,2);
padvar =4;
padCOMPmean = padarray(QCompositeNino-QMean,[padvar,padvar],nan,'both');

moreNINO = padCOMPmean(5:end-4,5:end-4)>QhighNINO;
lessNINO = padCOMPmean(5:end-4,5:end-4)<QlowNINO;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

Qplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
m_contourf(lon,lat,Qplot','Color',[1 1 1],'LineWidth',3);
hold on
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',3,'markerfacecolor','w');
RuBU = cbrewer('div','RdBu',24);
colormap(aa2,RuBU);
title(strcat('Nino',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR Q'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv1 cv1])
ylabel(h1,'Q (kg/kg)');

aa3=subplot(3,4,3);
padvar =4;
padCOMPmean = padarray(QCompositeNina-QMean,[padvar,padvar],nan,'both');

moreNINO = padCOMPmean(5:end-4,5:end-4)>QhighNINA;
lessNINO = padCOMPmean(5:end-4,5:end-4)<QlowNINA;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

Qplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
m_contourf(lon,lat,Qplot','Color',[1,1,1],'LineWidth',3);
hold on
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',3,'markerfacecolor','w');
RuBU = cbrewer('div','RdBu',24);
colormap(aa3,RuBU);
title(strcat('Nina',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR Q'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv1 cv1])
ylabel(h1,'Q (kg/kg)');

aa4=subplot(3,4,4);
padvar =4;
padCOMPmean = padarray(QCompositeNoNino-QMean,[padvar,padvar],nan,'both');

moreNINO = padCOMPmean(5:end-4,5:end-4)>QhighNoNINO;
lessNINO = padCOMPmean(5:end-4,5:end-4)<QlowNoNINO;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);


Qplot = padCOMPmean;
QplotInt =  padarray((QplotIntNoNino),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
m_contourf(lon,lat,Qplot','Color',[1,1,1],'LineWidth',3);
hold on
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',3,'markerfacecolor','w');
RuBU = cbrewer('div','RdBu',24);
colormap(aa4,RuBU);
title(strcat('Nuetral',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR Q'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv1 cv1])
ylabel(h1,'Q (kg/kg)');


aa5=subplot(3,4,5);
padvar =4;
padCOMPmean = padarray(UMean,[padvar,padvar],nan,'both');
Uplot = padCOMPmean;
UplotInt =  padarray((UplotInt),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,Uplot','Color',[1,1,1],'LineWidth',3);
hold on
m_contourf(lon,lat,UplotInt',24);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Blues',24);
colormap(aa5,RuBU);
title(strcat(titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR U'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
cv2 =max(max(UplotInt));
caxis([0 cv2])
ylabel(h1,'Std Dev of U wind (m/s)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find color scale%

fmax = cat(3,UCompositeNino-UMean,UCompositeNina-UMean,UCompositeNoNino-UMean);
fmax(fmax==Inf)=NaN;
fmax(fmax==-Inf)=NaN;
cv2 = max(max(max(abs(fmax))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aa6=subplot(3,4,6);
padvar =4;
padCOMPmean = padarray(UCompositeNino-UMean,[padvar,padvar],nan,'both');

moreNINO = padCOMPmean(5:end-4,5:end-4)>UhighNINO;
lessNINO = padCOMPmean(5:end-4,5:end-4)<UlowNINO;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

Uplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
m_contourf(lon,lat,Uplot','Color',[1,1,1],'LineWidth',3);
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',3,'markerfacecolor','w');
hold on
RuBU = cbrewer('div','BrBG',24);
colormap(aa6,RuBU);
title(strcat('Nino',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR U'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv2, cv2])
ylabel(h1,'U wind (m/s)');


aa7 = subplot(3,4,7);
padvar =4;
padCOMPmean = padarray(UCompositeNina-UMean,[padvar,padvar],nan,'both');

moreNINO = padCOMPmean(5:end-4,5:end-4)>UhighNINA;
lessNINO = padCOMPmean(5:end-4,5:end-4)<UlowNINA;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

Uplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
m_contourf(lon,lat,Uplot','Color',[1,1,1],'LineWidth',3);
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',3,'markerfacecolor','w');
hold on
RuBU = cbrewer('div','BrBG',24);
colormap(aa7,RuBU);
title(strcat('Nina',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR U'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv2 cv2])
ylabel(h1,'U wind (m/s)');


aa8=subplot(3,4,8);
padvar =4;
padCOMPmean = padarray(UCompositeNoNino-UMean,[padvar,padvar],nan,'both');

moreNINO = padCOMPmean(5:end-4,5:end-4)>UhighNoNINO;
lessNINO = padCOMPmean(5:end-4,5:end-4)<UlowNoNINO;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

Uplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
m_contourf(lon,lat,Uplot','Color',[1,1,1],'LineWidth',3);
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',3,'markerfacecolor','w');
hold on
RuBU = cbrewer('div','BrBG',24);
colormap(aa8,RuBU);
title(strcat('Nuetral',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR U'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv2 cv2])
ylabel(h1,'U wind (m/s)');



aa9=subplot(3,4,9);
padvar =4;
padCOMPmean = padarray(VMean,[padvar,padvar],nan,'both');
Vplot = padCOMPmean;
VplotInt =  padarray((VplotInt),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,Vplot','Color',[1,1,1],'LineWidth',3);
hold on
m_contourf(lon,lat,VplotInt',24);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Greens',24);
colormap(aa9,RuBU);
title(strcat(titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR V'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
cv3 =max(max(VplotInt));
caxis([0 cv3])
ylabel(h1,'Std Dev of V wind (m/s)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find color scale%

fmax = cat(3,VCompositeNino-VMean,VCompositeNina-VMean,VCompositeNoNino-VMean);
fmax(fmax==Inf)=NaN;
fmax(fmax==-Inf)=NaN;
cv3 = max(max(max(abs(fmax))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


aa10=subplot(3,4,10);
padvar =4;
padCOMPmean = padarray(VCompositeNino-VMean,[padvar,padvar],nan,'both');

moreNINO = padCOMPmean(5:end-4,5:end-4)>VhighNINO;
lessNINO = padCOMPmean(5:end-4,5:end-4)<VlowNINO;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

Vplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
m_contourf(lon,lat,Vplot','Color',[1,1,1],'LineWidth',3);
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',3,'markerfacecolor','w');
hold on
RuBU = cbrewer('div','PRGn',24);
colormap(aa10,RuBU);
title(strcat('Nino',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR V'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv3 cv3])
ylabel(h1,'V wind (m/s)');

aa11 = subplot(3,4,11);
padvar =4;
padCOMPmean = padarray(VCompositeNina-VMean,[padvar,padvar],nan,'both');

moreNINO = padCOMPmean(5:end-4,5:end-4)>VhighNINA;
lessNINO = padCOMPmean(5:end-4,5:end-4)<VlowNINA;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

Vplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
m_contourf(lon,lat,Vplot','Color',[1,1,1],'LineWidth',3);
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',3,'markerfacecolor','w');
hold on
RuBU = cbrewer('div','PRGn',24);
colormap(aa11,RuBU);
title(strcat('Nina',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR V'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv3 cv3])
ylabel(h1,'V wind (m/s)');


aa12=subplot(3,4,12);
padvar =4;
padCOMPmean = padarray(VCompositeNoNino-VMean,[padvar,padvar],nan,'both');

moreNINO = padCOMPmean(5:end-4,5:end-4)>UhighNoNINO;
lessNINO = padCOMPmean(5:end-4,5:end-4)<UlowNoNINO;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

Vplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
m_contourf(lon,lat,Vplot','Color',[1,1,1],'LineWidth',3);
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',3,'markerfacecolor','w');
hold on
RuBU = cbrewer('div','PRGn',24);
colormap(aa12,RuBU);
title(strcat('Nuetral',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR V'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv3 cv3])
ylabel(h1,'V wind (m/s)');

set(ax11,'position',[ 203          92        1856        1253])


assx = annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Mean', ...
    'EdgeColor', 'none', ...
    'Position', [0.11    0.9000    1.0000    0.1000],'FontSize',25,'FontWeight','bold');

assx = annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Anomaly', ...
    'EdgeColor', 'none', ...
    'Position', [0.58    0.9000    1.0000    0.1000],'FontSize',25,'FontWeight','bold');

annotation('line',[0.03 .258],[0.978 0.978],'color','k');
annotation('line',[0.27 .97],[0.978 0.978],'color','k');

figsave=strcat('/Users/wchapman/Desktop/QUVfigs/',num2str(titMon(2:end)),'_',num2str(Level(RRa)),'_Aanomaly.jpg');
saveas(ax11,figsave);
close all


%% Fig Panel 3
% Make Plotting Variables && plot 
QmeanYear = nanmean(Qperiod,3);
QmeanNino = nanmean(QcatMemNino,3);
QmeanNina = nanmean(QcatMemNina,3);
QmeanNoNino = nanmean(QcatMemNoNino,3);

QCompInt = Qperiod - QmeanYear;
QCompIntNino = QcatMemNino - QmeanNino;
QCompIntNina = QcatMemNina - QmeanNina;
QCompIntNoNino = QcatMemNoNino - QmeanNoNino;

QVarCompInt = nanvar(QCompInt,0,3);
QVarCompInt = (squeeze(QVarCompInt)).^0.5;
QStdCompInt = mean(squeeze(QVarCompInt),3);
    
QVarCompIntNino = nanvar(QCompIntNino,0,3);
QVarCompIntNino = (squeeze(QVarCompIntNino)).^0.5;
QStdCompIntNino = mean(squeeze(QVarCompIntNino),3);
    
QVarCompIntNina = nanvar(QCompIntNina,0,3);
QVarCompIntNina = (squeeze(QVarCompIntNina)).^0.5;
QStdCompIntNina = mean(squeeze(QVarCompIntNina),3);

QVarCompIntNoNino = nanvar(QCompIntNoNino,0,3);
QVarCompIntNoNino = (squeeze(QVarCompIntNoNino)).^0.5;
QStdCompIntNoNino = mean(squeeze(QVarCompIntNoNino),3);

QplotInt = QStdCompInt;
QplotIntNino = QStdCompIntNino;
QplotIntNina = QStdCompIntNina;
QplotIntNoNino = QStdCompIntNoNino;


UmeanYear = nanmean(Uperiod,3);
UmeanNino = nanmean(UcatMemNino,3);
UmeanNina = nanmean(UcatMemNina,3);
UmeanNoNino = nanmean(UcatMemNoNino,3);


UCompInt = Uperiod - UmeanYear;
UCompIntNino = UcatMemNino - UmeanNino;
UCompIntNina = UcatMemNina - UmeanNina;
UCompIntNoNino = UcatMemNoNino - UmeanNoNino;


UVarCompInt = nanvar(UCompInt,0,3);
UVarCompInt = (squeeze(UVarCompInt)).^0.5;
UStdCompInt = mean(squeeze(UVarCompInt),3);
    
UVarCompIntNino = nanvar(UCompIntNino,0,3);
UVarCompIntNino = (squeeze(UVarCompIntNino)).^0.5;
UStdCompIntNino = mean(squeeze(UVarCompIntNino),3);
    
UVarCompIntNina = nanvar(UCompIntNina,0,3);
UVarCompIntNina = (squeeze(UVarCompIntNina)).^0.5;
UStdCompIntNina = mean(squeeze(UVarCompIntNina),3);

UVarCompIntNoNino = nanvar(UCompIntNoNino,0,3);
UVarCompIntNoNino = (squeeze(UVarCompIntNoNino)).^0.5;
UStdCompIntNoNino = mean(squeeze(UVarCompIntNoNino),3);

UplotInt = UStdCompInt;
UplotIntNino = UStdCompIntNino;
UplotIntNina = UStdCompIntNina;
UplotIntNoNino = UStdCompIntNoNino;


VmeanYear = nanmean(Vperiod,3);
VmeanNino = nanmean(VcatMemNino,3);
VmeanNina = nanmean(VcatMemNina,3);
VmeanNoNino = nanmean(VcatMemNoNino,3);


VCompInt = Vperiod - VmeanYear;
VCompIntNino = VcatMemNino - VmeanNino;
VCompIntNina = VcatMemNina - VmeanNina;
VCompIntNoNino = VcatMemNoNino - VmeanNoNino;


VVarCompInt = nanvar(VCompInt,0,3);
VVarCompInt = (squeeze(VVarCompInt)).^0.5;
VStdCompInt = mean(squeeze(VVarCompInt),3);
    
VVarCompIntNino = nanvar(VCompIntNino,0,3);
VVarCompIntNino = (squeeze(VVarCompIntNino)).^0.5;
VStdCompIntNino = mean(squeeze(VVarCompIntNino),3);
    
VVarCompIntNina = nanvar(VCompIntNina,0,3);
VVarCompIntNina = (squeeze(VVarCompIntNina)).^0.5;
VStdCompIntNina = mean(squeeze(VVarCompIntNina),3);

VVarCompIntNoNino = nanvar(VCompIntNoNino,0,3);
VVarCompIntNoNino = (squeeze(VVarCompIntNoNino)).^0.5;
VStdCompIntNoNino = mean(squeeze(VVarCompIntNoNino),3);

VplotInt = VStdCompInt;
VplotIntNino = VStdCompIntNino;
VplotIntNina = VStdCompIntNina;
VplotIntNoNino = VStdCompIntNoNino;

load('Q_ARmask300____19520101_19521231.mat');
titMon = 'o';
mons =['J','F','M','A','M','J','J','A','S','O','N','D'];
for yy = 1:length(monthsall)
    
    mon = monthsall(yy);
    
    monC = mons(mon);
    
    titMon = strcat(titMon,monC);      
end 
ax11=figure;
aa1=subplot(3,4,1);
padvar =4;
padCOMPmean = padarray(QMean,[padvar,padvar],nan,'both');
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

Qplot = padCOMPmean;
Qdivide=QplotInt;
QplotInt =  padarray((QplotInt),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
[Cmap,Ctext] = m_contour(lon,lat,Qplot','Color',[1 1 1],'LineWidth',3);
hold on
m_contourf(lon,lat,QplotInt',24);
clabel(Cmap,Ctext,'FontSize',15,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Oranges',24);
colormap(aa1,RuBU);
title(strcat(titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR Q'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
cv1 =max(max(QplotInt));
caxis([0 cv1])
ylabel(h1,'Std Dev of Q (kg/kg)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find color scale%

fmax = cat(3,(QCompositeNino-QMean)./Qdivide,(QCompositeNina-QMean)./Qdivide,(QCompositeNoNino-QMean)./Qdivide);
fmax(fmax==Inf)=NaN;
fmax(fmax==-Inf)=NaN;
cv1 = max(max(max(abs(fmax))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aa2=subplot(3,4,2);
padvar =4;
padCOMPmean = padarray((QCompositeNino-QMean)./Qdivide,[padvar,padvar],nan,'both');

moreNINO = (QCompositeNino-QMean)>QhighNINO;
lessNINO = (QCompositeNino-QMean)<QlowNINO;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

Qplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,Qplot','Color',[1 1 1],'LineWidth',3);
hold on
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',2,'markerfacecolor','w');
RuBU = cbrewer('div','RdBu',24);
colormap(aa2,RuBU);
title(strcat('Nino',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR Q'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv1 cv1])

aa3=subplot(3,4,3);
padvar =4;
padCOMPmean = padarray((QCompositeNina-QMean)./Qdivide,[padvar,padvar],nan,'both');

moreNINO = (QCompositeNina-QMean)>QhighNINA;
lessNINO = (QCompositeNina-QMean)<QlowNINA;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

Qplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,Qplot','Color',[1,1,1],'LineWidth',3);
hold on
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',2,'markerfacecolor','w');
RuBU = cbrewer('div','RdBu',24);
colormap(aa3,RuBU);
title(strcat('Nina',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR Q'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv1 cv1])

aa4=subplot(3,4,4);
padvar =4;
padCOMPmean = padarray((QCompositeNoNino-QMean)./Qdivide,[padvar,padvar],nan,'both');

moreNINO = (QCompositeNoNino-QMean)>QhighNoNINO;
lessNINO = (QCompositeNoNino-QMean)<QlowNoNINO;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);


Qplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,Qplot','Color',[1,1,1],'LineWidth',3);
hold on
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',2,'markerfacecolor','w');
RuBU = cbrewer('div','RdBu',24);
colormap(aa4,RuBU);
title(strcat('Nuetral',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR Q'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv1 cv1])


aa5=subplot(3,4,5);
padvar =4;
padCOMPmean = padarray(UMean,[padvar,padvar],nan,'both');
Uplot = padCOMPmean;
Udivide = UplotInt;
UplotInt =  padarray((UplotInt),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,Uplot','Color',[1,1,1],'LineWidth',3);
hold on
m_contourf(lon,lat,UplotInt',24);
clabel(Cmap,Ctext,'FontSize',15,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Blues',24);
colormap(aa5,RuBU);
title(strcat(titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR U'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
cv2 =max(max(UplotInt));
caxis([0 cv2])
ylabel(h1,'Std Dev of U wind (m/s)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find color scale%

fmax = cat(3,(UCompositeNino-UMean)./Udivide,(UCompositeNina-UMean)./Udivide,(UCompositeNoNino-UMean)./Udivide);
fmax(fmax==Inf)=NaN;
fmax(fmax==-Inf)=NaN;
cv2 = max(max(max(abs(fmax))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aa6=subplot(3,4,6);
padvar =4;
padCOMPmean = padarray((UCompositeNino-UMean)./Udivide,[padvar,padvar],nan,'both');

moreNINO = (UCompositeNino-UMean)>UhighNINO;
lessNINO = (UCompositeNino-UMean)<UlowNINO;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

Uplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
m_contourf(lon,lat,Uplot','Color',[1,1,1],'LineWidth',3);
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',2,'markerfacecolor','w');
hold on
RuBU = cbrewer('div','RdBu',24);
colormap(aa6,RuBU);
title(strcat('Nino',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR U'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv2, cv2])


aa7 = subplot(3,4,7);
padvar =4;
padCOMPmean = padarray((UCompositeNina-UMean)./Udivide,[padvar,padvar],nan,'both');

moreNINO = (UCompositeNina-UMean)>UhighNINA;
lessNINO = (UCompositeNina-UMean)<UlowNINA;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

Uplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,Uplot','Color',[1,1,1],'LineWidth',3);
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',2,'markerfacecolor','w');
hold on
RuBU = cbrewer('div','RdBu',24);
colormap(aa7,RuBU);
title(strcat('Nina',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR U'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv2 cv2])


aa8=subplot(3,4,8);
padvar =4;
padCOMPmean = padarray((UCompositeNoNino-UMean)./Udivide,[padvar,padvar],nan,'both');

moreNINO = (UCompositeNoNino-UMean)>UhighNoNINO;
lessNINO = (UCompositeNoNino-UMean)<UlowNoNINO;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

Uplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,Uplot','Color',[1,1,1],'LineWidth',3);
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',2,'markerfacecolor','w');
hold on
RuBU = cbrewer('div','RdBu',24);
colormap(aa8,RuBU);
title(strcat('Nuetral',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR U'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv2 cv2])



aa9=subplot(3,4,9);
padvar =4;
padCOMPmean = padarray(VMean,[padvar,padvar],nan,'both');
Vplot = padCOMPmean;
Vdivide = VplotInt;
VplotInt =  padarray((VplotInt),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
[Cmap,Ctext] = m_contour(lon,lat,Vplot','Color',[1,1,1],'LineWidth',3);
hold on
m_contourf(lon,lat,VplotInt',24);
clabel(Cmap,Ctext,'FontSize',15,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Greens',24);
colormap(aa9,RuBU);
title(strcat(titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR V'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
cv3 =max(max(VplotInt));
caxis([0 cv3])
ylabel(h1,'Std Dev of V wind (m/s)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find color scale%

fmax = cat(3,(VCompositeNino-VMean)./Vdivide,(VCompositeNina-VMean)./Vdivide,(VCompositeNoNino-VMean)./Vdivide);
fmax(fmax==Inf)=NaN;
fmax(fmax==-Inf)=NaN;
cv3 = max(max(max(abs(fmax))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


aa10=subplot(3,4,10);
padvar =4;
padCOMPmean = padarray((VCompositeNino-VMean)./Vdivide,[padvar,padvar],nan,'both');

moreNINO = (VCompositeNino-VMean)>VhighNINO;
lessNINO = (VCompositeNino-VMean)<VlowNINO;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

Vplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
m_contourf(lon,lat,Vplot','Color',[1,1,1],'LineWidth',3);
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',2,'markerfacecolor','w');
hold on
RuBU = cbrewer('div','RdBu',24);
colormap(aa10,RuBU);
title(strcat('Nino',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR V'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv3 cv3])

aa11 = subplot(3,4,11);
padvar =4;
padCOMPmean = padarray((VCompositeNina-VMean)./Vdivide,[padvar,padvar],nan,'both');

moreNINO = (VCompositeNina-VMean)>VhighNINA;
lessNINO = (VCompositeNina-VMean)<VlowNINA;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

Vplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,Vplot','Color',[1,1,1],'LineWidth',3);
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',2,'markerfacecolor','w');
hold on
RuBU = cbrewer('div','RdBu',24);
colormap(aa11,RuBU);
title(strcat('Nina',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR V'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv3 cv3])

aa12=subplot(3,4,12);
padvar =4;
padCOMPmean = padarray((VCompositeNoNino-VMean)./Vdivide,[padvar,padvar],nan,'both');

moreNINO = (VCompositeNoNino-VMean)>UhighNoNINO;
lessNINO = (VCompositeNoNino-VMean)<UlowNoNINO;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

Vplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,Vplot','Color',[1,1,1],'LineWidth',3);
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',2,'markerfacecolor','w');
hold on
RuBU = cbrewer('div','RdBu',24);
colormap(aa12,RuBU);
title(strcat('Nuetral',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR V'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv3 cv3])

set(ax11,'position',[ 203          92        1856        1253])


annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Mean', ...
    'EdgeColor', 'none', ...
    'Position', [0.11    0.9000    1.0000    0.1000],'FontSize',25,'FontWeight','bold');

annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Standard Anomaly', ...
    'EdgeColor', 'none', ...
    'Position', [0.55    0.9000    1.0000    0.1000],'FontSize',25,'FontWeight','bold');

annotation('line',[0.03 .258],[0.978 0.978],'color','k');
annotation('line',[0.27 .97],[0.978 0.978],'color','k');

figsave=strcat('/Users/wchapman/Desktop/QUVfigs/',num2str(titMon(2:end)),'_',num2str(Level(RRa)),'_StdAnomaly.jpg');
saveas(ax11,figsave);
end







