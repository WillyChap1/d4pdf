%Composite 200hpa heights. 

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

close all
% =========================================================================
% =========================================================================
% =========================================================================

%==========================================================================
    %======================= Get Mean Response  =====================
%==========================================================================

cd /Users/wchapman/SIO/research/100Amip/d4pdf/200Ens/

%find the leap years: 
leapstart = 1952;
leapyears = zeros(1,31);
leapyears(1) = leapstart;

for jj = 2:30
    leapyears(jj) = leapyears(jj-1)+4; 
end

findfilesGP = strcat('Geop200tot','*.mat');
findfilesU = strcat('/Volumes/HotSlop/d4pdf/U200/U200tot','*.mat');
findfilesV = strcat('/Volumes/HotSlop/d4pdf/V200/V200tot','*.mat');
pp = dir(findfilesGP);
uu = dir(findfilesU);
vv = dir(findfilesV);


load(pp(1).name);
timetot = timer;

timekeep = [];
timekeepMonths =[];
%GPkeepAll = [];
GPkeepMonth = [];
GPperiod = [];
Uperiod =[];
Vperiod =[];
timekeepPeriod = [];

for hh = 1:(length(pp)-2)
    GPperiodtemp = [];
    timePeriod = [];
    Vperiodtemp = [];
    Uperiodtemp = [];
    
    if isempty(months1) == 0
    disp(strcat('Im in year:',{' '},pp(hh).name(end-11:end-8)))
    fileGP = pp(hh).name;
    fileU = uu(hh).name;
    fileV = vv(hh).name;
    load(fileGP);
    load(strcat('/Volumes/HotSlop/d4pdf/U200/',fileU));
    load(strcat('/Volumes/HotSlop/d4pdf/V200/',fileV));
    
    timetot = timer;
        
    for jr = min(months1):max(months1)    %Dec Nov. 
        [~,month,~]=datevec(timetot);
        indexdo = find(month==jr);
        GPtemp = GPtot(:,:,:,indexdo);
        Utemp = U200tot(:,:,:,indexdo);
        Vtemp = v200tot(:,:,:,indexdo);
        
        GPkeepMonth = cat(4,GPkeepMonth,sum(GPtemp,4));
        %GPkeepAll = cat(4,GPkeepAll,GPtemp);
        timekeepMonths = cat(1,timekeepMonths,timetot(indexdo(1)));
        timekeep = cat(1,timekeep,timetot(indexdo));
        
        GPperiodtemp = cat(4,GPperiodtemp,GPtemp);   
        Vperiodtemp = cat(4,Vperiodtemp,Vtemp);   
        Uperiodtemp = cat(4,Uperiodtemp,Utemp);   
        timePeriod = cat(1,timePeriod,timetot(indexdo(1)));
        
    end

    end
    
    if isempty(months2)==0
  
    fileGP = pp(hh+1).name;    
    timetot = timer;
        
    for jr = min(months2):max(months2)    %Dec Nov. 
        [~,month,~]=datevec(timetot);
        indexdo = find(month==jr);
        GPtemp = GPtot(:,:,:,indexdo);
        Utemp = U200tot(:,:,:,indexdo);
        Vtemp = v200tot(:,:,:,indexdo);
        
        
        GPkeepMonth = cat(4,GPkeepMonth,sum(GPtemp,4)); 
        %GPkeepAll = cat(4,GPkeepAll,GPtemp);
        timekeepMonths = cat(1,timekeepMonths,timetot(indexdo(1)));
        timekeep = cat(1,timekeep,timetot(indexdo));
        timePeriod = cat(1,timePeriod,timetot(indexdo(1)));
        
        GPperiodtemp = cat(4,GPperiodtemp,GPtemp);
        Vperiodtemp = cat(4,Vperiodtemp,Vtemp);   
        Uperiodtemp = cat(4,Uperiodtemp,Utemp);   
    end
        
    end
    
    GPperiod = cat(4,GPperiod,nanmean(GPperiodtemp,4));
    Vperiod = cat(4,Vperiod,nanmean(Vperiodtemp,4));
    Uperiod = cat(4,Uperiod,nanmean(Uperiodtemp,4));
    timekeepPeriod = cat(1,timekeepPeriod,timePeriod(1));
      
end 


GPMean = nanmean(GPperiod,3);
GPMean = squeeze(GPMean);
GPMean = nanmean(GPMean,3);

UMean = nanmean(Uperiod,3);
UMean = squeeze(UMean);
UMean = nanmean(UMean,3);

VMean = nanmean(Vperiod,3);
VMean = squeeze(VMean);
VMean = nanmean(VMean,3);


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
    
    GPcatMemNino(:,:,:,ff) = GPperiod(:,:,:,NinoIndex); 
    VcatMemNino(:,:,:,ff) = Vperiod(:,:,:,NinoIndex); 
    UcatMemNino(:,:,:,ff) = Uperiod(:,:,:,NinoIndex); 
    timeNino(ff) = timekeepPeriod(NinoIndex);
    
end

GPCompositeNino = nanmean(GPcatMemNino,3);
GPCompositeNino = squeeze(GPCompositeNino);
GPCompositeNino = nanmean(GPCompositeNino,3);
GPCompositeNino = squeeze(GPCompositeNino);

VCompositeNino = nanmean(VcatMemNino,3);
VCompositeNino = squeeze(VCompositeNino);
VCompositeNino = nanmean(VCompositeNino,3);
VCompositeNino = squeeze(VCompositeNino);

UCompositeNino = nanmean(UcatMemNino,3);
UCompositeNino = squeeze(UCompositeNino);
UCompositeNino = nanmean(UCompositeNino,3);
UCompositeNino = squeeze(UCompositeNino);
%==========================================================================
%========================       Nina .       ==============================
%==========================================================================

for ff = 1:length(Ninayears)
    
    NinaIndex = find(yearsIn==Ninayears(ff));
    
    GPcatMemNina(:,:,:,ff) = GPperiod(:,:,:,NinaIndex); 
    VcatMemNina(:,:,:,ff) = Vperiod(:,:,:,NinaIndex); 
    UcatMemNina(:,:,:,ff) = Uperiod(:,:,:,NinaIndex);
    timeNina(ff) = timekeepPeriod(NinaIndex);
    
end

GPCompositeNina = nanmean(GPcatMemNina,3);
GPCompositeNina = squeeze(GPCompositeNina);
GPCompositeNina = nanmean(GPCompositeNina,3);
GPCompositeNina = squeeze(GPCompositeNina);

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
    
    GPcatMemNoNino(:,:,:,ff) = GPperiod(:,:,:,NoNinaIndex);
    VcatMemNoNino(:,:,:,ff) = Vperiod(:,:,:,NoNinaIndex); 
    UcatMemNoNino(:,:,:,ff) = Uperiod(:,:,:,NoNinaIndex); 
    
    timeNina(ff) = timekeepPeriod(NoNinaIndex);
    
end

GPCompositeNoNino = nanmean(GPcatMemNoNino,3);
GPCompositeNoNino = squeeze(GPCompositeNoNino);
GPCompositeNoNino = nanmean(GPCompositeNoNino,3);
GPCompositeNoNino = squeeze(GPCompositeNoNino);

VCompositeNoNino = nanmean(VcatMemNoNino,3);
VCompositeNoNino = squeeze(VCompositeNoNino);
VCompositeNoNino = nanmean(VCompositeNoNino,3);
VCompositeNoNino = squeeze(VCompositeNoNino);

UCompositeNoNino = nanmean(UcatMemNoNino,3);
UCompositeNoNino = squeeze(UCompositeNoNino);
UCompositeNoNino = nanmean(UCompositeNoNino,3);
UCompositeNoNino = squeeze(UCompositeNoNino);

%%
disp('Entering BootStrap')
[GPlowNINO,GPhighNINO] = Bootstrap4D(GPperiod,size(GPcatMemNino,4)/size(GPperiod,4),10000,[1,99]);
[GPlowNINA,GPhighNINA] = Bootstrap4D(GPperiod,size(GPcatMemNina,4)/size(GPperiod,4),10000,[1,99]);
[GPlowNoNINO,GPhighNoNINO] = Bootstrap4D(GPperiod,size(GPcatMemNina,4)/size(GPperiod,4),10000,[1,99]);
%%

GPmeanYear = nanmean(GPperiod,3);
GPmeanNino = nanmean(GPcatMemNino,3);
GPmeanNina = nanmean(GPcatMemNina,3);
GPmeanNoNino = nanmean(GPcatMemNoNino,3);

UmeanYear = nanmean(Uperiod,3);
UmeanNino = nanmean(UcatMemNino,3);
UmeanNina = nanmean(UcatMemNina,3);
UmeanNoNino = nanmean(UcatMemNoNino,3);

VmeanYear = nanmean(Vperiod,3);
VmeanNino = nanmean(VcatMemNino,3);
VmeanNina = nanmean(VcatMemNina,3);
VmeanNoNino = nanmean(VcatMemNoNino,3);

GPCompInt = GPperiod - GPmeanYear;
GPCompIntNino = GPcatMemNino - GPmeanNino;
GPCompIntNina = GPcatMemNina - GPmeanNina;
GPCompIntNoNino = GPcatMemNoNino - GPmeanNoNino;

GPVarCompInt = nanvar(GPCompInt,0,3);
GPVarCompInt = (squeeze(GPVarCompInt)).^0.5;
GPStdCompInt = mean(squeeze(GPVarCompInt),3);
    
GPVarCompIntNino = nanvar(GPCompIntNino,0,3);
GPVarCompIntNino = (squeeze(GPVarCompIntNino)).^0.5;
GPStdCompIntNino = mean(squeeze(GPVarCompIntNino),3);
    
GPVarCompIntNina = nanvar(GPCompIntNina,0,3);
GPVarCompIntNina = (squeeze(GPVarCompIntNina)).^0.5;
GPStdCompIntNina = mean(squeeze(GPVarCompIntNina),3);

GPVarCompIntNoNino = nanvar(GPCompIntNoNino,0,3);
GPVarCompIntNoNino = (squeeze(GPVarCompIntNoNino)).^0.5;
GPStdCompIntNoNino = mean(squeeze(GPVarCompIntNoNino),3);

GPplotInt = GPStdCompInt;
GPplotIntNino = GPStdCompIntNino;
GPplotIntNina = GPStdCompIntNina;
GPplotIntNoNino = GPStdCompIntNoNino;




                                %options: [gap]   %hght [Top Bot]  %Wid  [L R]                          
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
if ~make_it_tight,  clear subplot;  end


load('Geop200tot______19520101_19521231.mat');
titMon = [];
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
fmax = cat(3,GPplotIntNino,GPplotIntNina,GPplotIntNoNino);
fmax(fmax==Inf)=NaN;
fmax(fmax==-Inf)=NaN;
cv12 = max(max(max(abs(fmax))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


aa2=subplot(3,4,2);
padvar =4;
uPad = padarray(UCompositeNino,[padvar,padvar],nan,'both');
vPad = padarray(VCompositeNino,[padvar,padvar],nan,'both');
uPad=uPad.*3.6;
vPad=vPad.*3.6;
padCOMPmean = padarray(GPCompositeNino,[padvar,padvar],nan,'both');

GPplot = padCOMPmean;
GPplotInt =  padarray((GPplotIntNino),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
[Cmap,Ctext] = m_contour(lon,lat,GPplot','Color',[1 1 1],'LineWidth',3);
hold on
m_contourf(lon,lat,GPplotInt',24);
m_barbs(lon,lat,uPad,vPad,7,15,0.09,1,[1.0000    0.9294    0.6275]);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Blues',24);
colormap(aa2,RuBU);
title(strcat('Nino',{' '},titMon,{' '},'200 hPa GPH'))
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
caxis([0 cv12])
ylabel(h1,'Std Dev of GPH [m]');

aa3=subplot(3,4,3);
padvar =4;
uPad = padarray(UCompositeNina,[padvar,padvar],nan,'both');
vPad = padarray(VCompositeNina,[padvar,padvar],nan,'both');
uPad=uPad.*3.6;
vPad=vPad.*3.6;
padCOMPmean = padarray(GPCompositeNina,[padvar,padvar],nan,'both');
GPplot = padCOMPmean;
GPplotInt =  padarray((GPplotIntNina),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
[Cmap,Ctext] = m_contour(lon,lat,GPplot','Color',[1,1,1],'LineWidth',3);
hold on
m_contourf(lon,lat,GPplotInt',24);
m_barbs(lon,lat,uPad,vPad,7,15,0.09,1,[1.0000    0.9294    0.6275]);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Blues',24);
colormap(aa3,RuBU);
title(strcat('Nina',{' '},titMon,{' '},'200 hPa GPH'))
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
caxis([0 cv12])
ylabel(h1,'Std Dev of GPH [m]');

aa4=subplot(3,4,4);
padvar =4;
uPad = padarray(UCompositeNoNino,[padvar,padvar],nan,'both');
vPad = padarray(VCompositeNoNino,[padvar,padvar],nan,'both');
uPad=uPad.*3.6;
vPad=vPad.*3.6;
padCOMPmean = padarray(GPCompositeNoNino,[padvar,padvar],nan,'both');
GPplot = padCOMPmean;
GPplotInt =  padarray((GPplotIntNoNino),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
[Cmap,Ctext] = m_contour(lon,lat,GPplot','Color',[1,1,1],'LineWidth',3);
hold on
m_contourf(lon,lat,GPplotInt',24);
m_barbs(lon,lat,uPad,vPad,7,15,0.09,1,[1.0000    0.9294    0.6275]);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Blues',24);
colormap(aa4,RuBU);
title(strcat('Nuetral',{' '},titMon,{' '},'200 hPa GPH'))
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
caxis([0 cv12])
ylabel(h1,'Std Dev of GPH [m]');



% fig panel 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Plotting Variables && plot 
GPmeanYear = nanmean(GPperiod,3);
GPmeanNino = nanmean(GPcatMemNino,3);
GPmeanNina = nanmean(GPcatMemNina,3);
GPmeanNoNino = nanmean(GPcatMemNoNino,3);

GPCompInt = GPperiod - GPmeanYear;
GPCompIntNino = GPcatMemNino - GPmeanNino;
GPCompIntNina = GPcatMemNina - GPmeanNina;
GPCompIntNoNino = GPcatMemNoNino - GPmeanNoNino;

GPVarCompInt = nanvar(GPCompInt,0,3);
GPVarCompInt = (squeeze(GPVarCompInt)).^0.5;
GPStdCompInt = mean(squeeze(GPVarCompInt),3);
    
GPVarCompIntNino = nanvar(GPCompIntNino,0,3);
GPVarCompIntNino = (squeeze(GPVarCompIntNino)).^0.5;
GPStdCompIntNino = mean(squeeze(GPVarCompIntNino),3);
    
GPVarCompIntNina = nanvar(GPCompIntNina,0,3);
GPVarCompIntNina = (squeeze(GPVarCompIntNina)).^0.5;
GPStdCompIntNina = mean(squeeze(GPVarCompIntNina),3);

GPVarCompIntNoNino = nanvar(GPCompIntNoNino,0,3);
GPVarCompIntNoNino = (squeeze(GPVarCompIntNoNino)).^0.5;
GPStdCompIntNoNino = mean(squeeze(GPVarCompIntNoNino),3);

GPplotInt = GPStdCompInt;
GPplotIntNino = GPStdCompIntNino;
GPplotIntNina = GPStdCompIntNina;
GPplotIntNoNino = GPStdCompIntNoNino;

load('Geop200tot______19520101_19521231.mat');
titMon = [];
mons =['J','F','M','A','M','J','J','A','S','O','N','D'];
for yy = 1:length(monthsall)
    
    mon = monthsall(yy);
    
    monC = mons(mon);
    
    titMon = strcat(titMon,monC);      
end 
aa1=subplot(3,4,5);
padvar =4;
padCOMPmean = padarray(GPMean,[padvar,padvar],nan,'both');
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

GPplot = padCOMPmean;
GPdivide = GPplotInt;
GPplotInt =  padarray((GPplotInt),[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
[Cmap,Ctext] = m_contour(lon,lat,GPplot','Color',[1 1 1],'LineWidth',3);
hold on
m_contourf(lon,lat,GPplotInt',24);
clabel(Cmap,Ctext,'FontSize',20,'Color',[1,1,1],'LabelSpacing',500);
RuBU = cbrewer('seq','Blues',24);
colormap(aa1,RuBU);
title(strcat(titMon,{' '},'200 hPa GPH'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([0 cv12])
ylabel(h1,'Std Dev of GPH [m]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find color scale%

fmax = cat(3,GPCompositeNino-GPMean,GPCompositeNina-GPMean,GPCompositeNoNino-GPMean);
fmax(fmax==Inf)=NaN;
fmax(fmax==-Inf)=NaN;
cv1 = max(max(max(abs(fmax))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aa2=subplot(3,4,6);
padvar =4;
uPad = padarray(UCompositeNino - UMean,[padvar,padvar],nan,'both');
vPad = padarray(VCompositeNino - VMean,[padvar,padvar],nan,'both');
uPad=uPad.*3.6;
vPad=vPad.*3.6;
padCOMPmean = padarray(GPCompositeNino-GPMean,[padvar,padvar],nan,'both');

moreNINO = padCOMPmean(5:end-4,5:end-4)>GPhighNINO;
lessNINO = padCOMPmean(5:end-4,5:end-4)<GPlowNINO;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

GPplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,GPplot','Color',[1 1 1],'LineWidth',3);
hold on
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',3,'markerfacecolor','w');
m_barbs(lon,lat,uPad,vPad,9,15,0.07,2,'k');
RuBU = cbrewer('div','RdBu',24);
colormap(aa2,RuBU);
title(strcat('Nino',{' '},titMon,{' '},'200 hPa GPH'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv1 cv1])
ylabel(h1,'GPH [m]');

aa3=subplot(3,4,7);
padvar =4;
uPad = padarray(UCompositeNina - UMean,[padvar,padvar],nan,'both');
vPad = padarray(VCompositeNina - VMean,[padvar,padvar],nan,'both');
uPad=uPad.*3.6;
vPad=vPad.*3.6;
padCOMPmean = padarray(GPCompositeNina-GPMean,[padvar,padvar],nan,'both');

moreNINO = padCOMPmean(5:end-4,5:end-4)>GPhighNINA;
lessNINO = padCOMPmean(5:end-4,5:end-4)<GPlowNINA;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

GPplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
m_contourf(lon,lat,GPplot','Color',[1,1,1],'LineWidth',3);
hold on
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','w','linest','none','markersize',3,'markerfacecolor','w');
m_barbs(lon,lat,uPad,vPad,9,15,0.07,2,'k');
RuBU = cbrewer('div','RdBu',24);
colormap(aa3,RuBU);
title(strcat('Nina',{' '},titMon,{' '},'200 hPa GPH'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv1 cv1])
ylabel(h1,'GPH [m]');

aa4=subplot(3,4,8);
padvar =4;
uPad = padarray(UCompositeNoNino - UMean,[padvar,padvar],nan,'both');
vPad = padarray(VCompositeNoNino - VMean,[padvar,padvar],nan,'both');
uPad=uPad.*3.6;
vPad=vPad.*3.6;
padCOMPmean = padarray(GPCompositeNoNino-GPMean,[padvar,padvar],nan,'both');

moreNINO = padCOMPmean(5:end-4,5:end-4)>GPhighNoNINO;
lessNINO = padCOMPmean(5:end-4,5:end-4)<GPlowNoNINO;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);


GPplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,GPplot','Color',[1,1,1],'LineWidth',3);
hold on
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','w','linest','none','markersize',3,'markerfacecolor','w');
m_barbs(lon,lat,uPad,vPad,9,15,0.07,2,'k');
RuBU = cbrewer('div','RdBu',24);
colormap(aa4,RuBU);
title(strcat('Nuetral',{' '},titMon,{' '},'200 hPa GPH'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv1 cv1])
ylabel(h1,'GPH [m]');


% Fig Panel 3
% Make Plotting Variables && plot 
GPmeanYear = nanmean(GPperiod,3);
GPmeanNino = nanmean(GPcatMemNino,3);
GPmeanNina = nanmean(GPcatMemNina,3);
GPmeanNoNino = nanmean(GPcatMemNoNino,3);

GPCompInt = GPperiod - GPmeanYear;
GPCompIntNino = GPcatMemNino - GPmeanNino;
GPCompIntNina = GPcatMemNina - GPmeanNina;
GPCompIntNoNino = GPcatMemNoNino - GPmeanNoNino;

GPVarCompInt = nanvar(GPCompInt,0,3);
GPVarCompInt = (squeeze(GPVarCompInt)).^0.5;
GPStdCompInt = mean(squeeze(GPVarCompInt),3);
    
GPVarCompIntNino = nanvar(GPCompIntNino,0,3);
GPVarCompIntNino = (squeeze(GPVarCompIntNino)).^0.5;
GPStdCompIntNino = mean(squeeze(GPVarCompIntNino),3);
    
GPVarCompIntNina = nanvar(GPCompIntNina,0,3);
GPVarCompIntNina = (squeeze(GPVarCompIntNina)).^0.5;
GPStdCompIntNina = mean(squeeze(GPVarCompIntNina),3);

GPVarCompIntNoNino = nanvar(GPCompIntNoNino,0,3);
GPVarCompIntNoNino = (squeeze(GPVarCompIntNoNino)).^0.5;
GPStdCompIntNoNino = mean(squeeze(GPVarCompIntNoNino),3);

GPplotInt = GPStdCompInt;
GPplotIntNino = GPStdCompIntNino;
GPplotIntNina = GPStdCompIntNina;
GPplotIntNoNino = GPStdCompIntNoNino;


load('Geop200tot______19520101_19521231.mat');
titMon = [];
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

fmax = cat(3,(GPCompositeNino-GPMean)./GPdivide,(GPCompositeNina-GPMean)./GPdivide,(GPCompositeNoNino-GPMean)./GPdivide);
fmax(fmax==Inf)=NaN;
fmax(fmax==-Inf)=NaN;
cv1 = max(max(max(abs(fmax))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aa2=subplot(3,4,10);
padvar =4;
padCOMPmean = padarray((GPCompositeNino-GPMean)./GPdivide,[padvar,padvar],nan,'both');
moreNINO = (GPCompositeNino-GPMean)>GPhighNINO;
lessNINO = (GPCompositeNino-GPMean)<GPlowNINO;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

GPplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,GPplot','Color',[1 1 1],'LineWidth',3);
hold on
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',2,'markerfacecolor','w');
RuBU = cbrewer('div','BrBG',24);
colormap(aa2,RuBU);
title(strcat('Nino',{' '},titMon,{' '},'200 hPa GPH'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv1 cv1])

aa3=subplot(3,4,11);
padvar =4;
padCOMPmean = padarray((GPCompositeNina-GPMean)./GPdivide,[padvar,padvar],nan,'both');

moreNINO = (GPCompositeNina-GPMean)>GPhighNINA;
lessNINO = (GPCompositeNina-GPMean)<GPlowNINA;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);

GPplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,GPplot','Color',[1,1,1],'LineWidth',3);
hold on
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',2,'markerfacecolor','w');
RuBU = cbrewer('div','BrBG',24);
colormap(aa3,RuBU);
title(strcat('Nina',{' '},titMon,{' '},'200 hPa GPH'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([-cv1 cv1])

aa4=subplot(3,4,12);
padvar =4;
padCOMPmean = padarray((GPCompositeNoNino-GPMean)./GPdivide,[padvar,padvar],nan,'both');

moreNINO = (GPCompositeNoNino-GPMean)>GPhighNoNINO;
lessNINO = (GPCompositeNoNino-GPMean)<GPlowNoNINO;

sigNino = moreNINO+lessNINO;
[rNsigNino, cNsigNino] = find(sigNino == 1);


GPplot = padCOMPmean;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,GPplot','Color',[1,1,1],'LineWidth',3);
hold on
m_plot(lon(rNsigNino+4),lat(cNsigNino+4),'marker','x','color','k','linest','none','markersize',2,'markerfacecolor','w');
RuBU = cbrewer('div','BrBG',24);
colormap(aa4,RuBU);
title(strcat('Neutral',{' '},titMon,{' '},'200 hPa GPH'))
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
    'Position', [0.55    0.257000    1.0000    0.1000],'FontSize',25,'FontWeight','bold');

annotation('line',[0.265 .97],[0.332 0.332],'color','k');

%vert line 


annotation('line',[0.2615 0.2615],[0.978 0.05],'color','k');



% figsave=strcat('/Users/wchapman/Desktop/QUVfigs/',num2str(titMon),'_',num2str(Level(RRa)),'_GPall.jpg');
% saveas(ax11,figsave);





