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

Level =[1000,925,850,700,600,500,400,300];
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



VTperiod = ((Qperiod.*Uperiod).^2 + (Qperiod.*Vperiod).^2).^.5;

VTMean = nanmean(VTperiod,3);
VTMean = squeeze(VTMean);
VTMean = nanmean(VTMean,3);


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
    
    VTcatMemNino(:,:,:,ff) = VTperiod(:,:,:,NinoIndex); 
    QcatMemNino(:,:,:,ff) = Qperiod(:,:,:,NinoIndex); 
    UcatMemNino(:,:,:,ff) = Uperiod(:,:,:,NinoIndex);
    VcatMemNino(:,:,:,ff) = Vperiod(:,:,:,NinoIndex);
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
    QcatMemNina(:,:,:,ff) = Qperiod(:,:,:,NinaIndex); 
    UcatMemNina(:,:,:,ff) = Uperiod(:,:,:,NinaIndex);
    VcatMemNina(:,:,:,ff) = Vperiod(:,:,:,NinaIndex);
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
    QcatMemNoNino(:,:,:,ff) = Qperiod(:,:,:,NoNinaIndex); 
    UcatMemNoNino(:,:,:,ff) = Uperiod(:,:,:,NoNinaIndex);
    VcatMemNoNino(:,:,:,ff) = Vperiod(:,:,:,NoNinaIndex);
    timeNina(ff) = timekeepPeriod(NoNinaIndex);
    
end

VTCompositeNoNino = nanmean(VTcatMemNoNino,3);
VTCompositeNoNino = squeeze(VTCompositeNoNino);
VTCompositeNoNino = nanmean(VTCompositeNoNino,3);
VTCompositeNoNino = squeeze(VTCompositeNoNino);

%==========================================================================
%==========================================================================
%==========================================================================
% figure
% pUfit = pU(1)*VTts+pU(2);
% plot((VTts),(Uts),'*')
% hold on 
% plot(VTts,pUfit)
% title('U')
% 
% figure
% pVfit = pV(1)*VTts+pV(2);
% plot(VTts,Vts,'*')
% hold on 
% plot(VTts,pVfit)
% title('V')
% 
% 
% figure
% pQfit = pQ(1)*VTts+pQ(2);
% plot(VTts,Qts,'*')
% hold on 
% plot(VTts,pQfit)
% title('Q')
% 
% figure
% pWSfit = pWS(1)*VTts+pWS(2);
% plot(VTts,WSts,'*')
% hold on 
% plot(VTts,pWSfit)
% title('WS')

VTperR = reshape(VTperiod,size(VTperiod,1),size(VTperiod,2),size(VTperiod,3)*size(VTperiod,4));

%U R2 map
UperR = reshape(Uperiod,size(Uperiod,1),size(Uperiod,2),size(Uperiod,3)*size(Uperiod,4));

R2U = zeros(size(Uperiod,1),size(Uperiod,2));
for ii = 1:size(Uperiod,1)
    for jj = 1:size(Uperiod,2)
       
        
       Uts = UperR(ii,jj,:);
       VTts = VTperR(ii,jj,:);
       pU = polyfit(VTts,Uts,1);
       pUfit = pU(1)*VTts+pU(2);
       Uresid = Uts - pUfit;
       uSSresid = sum(Uresid.^2);
       uSStotal = (length(Uts)-1) * var(Uts);
       R2U(ii,jj) = 1 - uSSresid/uSStotal;      
        
    end
    
end

%V R2 map
VperR = reshape(Vperiod,size(Vperiod,1),size(Vperiod,2),size(Vperiod,3)*size(Vperiod,4));

R2V = zeros(size(Vperiod,1),size(Vperiod,2));
for ii = 1:size(Vperiod,1)
    for jj = 1:size(Vperiod,2)
       
        
       Vts = VperR(ii,jj,:);
       VTts = VTperR(ii,jj,:);
       pV = polyfit(VTts,Vts,1);
       pVfit = pV(1)*VTts+pV(2);
       Vresid = Vts - pVfit;
       VSSresid = sum(Vresid.^2);
       VSStotal = (length(Vts)-1) * var(Vts);
       R2V(ii,jj) = 1 - VSSresid/VSStotal;      
        
    end
    
end

%Q R2 map
QperR = reshape(Qperiod,size(Qperiod,1),size(Qperiod,2),size(Qperiod,3)*size(Qperiod,4));

R2Q = zeros(size(Qperiod,1),size(Qperiod,2));
for ii = 1:size(Qperiod,1)
    for jj = 1:size(Qperiod,2)
       
        
       Qts = QperR(ii,jj,:);
       VTts = VTperR(ii,jj,:);
       pQ = polyfit(VTts,Qts,1);
       pQfit = pQ(1)*VTts+pQ(2);
       Qresid = Qts - pQfit;
       QSSresid = sum(Qresid.^2);
       QSStotal = (length(Qts)-1) * var(Qts);
       R2Q(ii,jj) = 1 - QSSresid/QSStotal;      
        
    end
    
end


VTts = reshape(VTperiod,size(VTperiod,1)*size(VTperiod,2)*size(VTperiod,3)*size(VTperiod,4),1);
Qts = reshape(Qperiod,size(Qperiod,1)*size(Qperiod,2)*size(Qperiod,3)*size(Qperiod,4),1);
Uts = reshape(Uperiod,size(Uperiod,1)*size(Uperiod,2)*size(Uperiod,3)*size(Uperiod,4),1);
Vts = reshape(Vperiod,size(Vperiod,1)*size(Vperiod,2)*size(Vperiod,3)*size(Vperiod,4),1);
WSts = (Uts.^2+Vts.^2).^0.5;

pU = polyfit(VTts,Uts,1);
pV = polyfit(VTts,Vts,1);
pQ = polyfit(VTts,Qts,1);
pWS = polyfit(VTts,WSts,1);

pWSfit = pWS(1)*VTts+pWS(2);
WSresid = WSts - pWSfit;
wsSSresid = sum(WSresid.^2);
wsSStotal = (length(WSts)-1) * var(WSts);
WSrsq = 1 - wsSSresid/wsSStotal

pQfit = pQ(1)*VTts+pQ(2);
Qresid = Qts - pQfit;
qSSresid = sum(Qresid.^2);
qSStotal = (length(Qts)-1) * var(Qts);
Qrsq = 1 - qSSresid/qSStotal

pUfit = pU(1)*VTts+pU(2);
Uresid = Uts - pUfit;
uSSresid = sum(Uresid.^2);
uSStotal = (length(Uts)-1) * var(Uts);
Ursq = 1 - uSSresid/uSStotal

pVfit = pV(1)*VTts+pV(2);
Vresid = Vts - pVfit;
vSSresid = sum(Vresid.^2);
vSStotal = (length(Vts)-1) * var(Vts);
Vrsq = 1 - vSSresid/vSStotal


%%
%to make shit look pretty.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Q_ARmask700____19530101_19531231.mat');
titMon = 'o';
mons =['J','F','M','A','M','J','J','A','S','O','N','D'];
for yy = 1:length(monthsall)
    
    mon = monthsall(yy);
    
    monC = mons(mon);
    
    titMon = strcat(titMon,monC);      
end

%Pad LAT/LON
padvar = 4;
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

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
if ~make_it_tight,  clear subplot;  end

ax11 = figure;
subplot(2,3,2)
pUfit = pU(1)*VTts+pU(2);
plot((VTts),(Uts),'*')
xlabel('Vapor Transport [Kg/Kg * m/s]','FontSize',18);
ylabel('U wind [m/s]','FontSize',18);
hold on 
plot(VTts,pUfit,'LineWidth',3)
tit = strcat('U',{'         '},'R^2:',{' '},num2str(Ursq));
title(tit,'FontSize',21)

subplot(2,3,3)
pVfit = pV(1)*VTts+pV(2);
plot(VTts,Vts,'*')
xlabel('Vapor Transport [Kg/Kg * m/s]','FontSize',18);
ylabel('V wind [m/s]','FontSize',18);
hold on 
plot(VTts,pVfit,'LineWidth',3)
tit = strcat('V',{'         '},'R^2:',{' '},num2str(Vrsq));
title(tit,'FontSize',21)

subplot(2,3,1)
pQfit = pQ(1)*VTts+pQ(2);
plot(VTts,Qts,'*')
xlabel('Vapor Transport [Kg/Kg * m/s]','FontSize',18);
ylabel('Specific Humidity [Kg/Kg]','FontSize',18);
hold on 
plot(VTts,pQfit,'LineWidth',3)
tit = strcat('Q',{'         '},'R^2:',{' '},num2str(Qrsq));
title(tit,'FontSize',21)

aa3 = subplot(2,3,4);
padvar =4;
Qplot = padarray(R2Q,[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',4,'Color','k');
m_contourf(lon,lat,Qplot','Color',[1,1,1],'LineWidth',1);
hold on
RuBU = cbrewer('seq','YlOrBr',24);
colormap(aa3,RuBU);
title(strcat('Q',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR R^2'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
ylabel(h1,'R^2')


aa3 = subplot(2,3,5);
padvar =4;
Uplot = padarray(R2U,[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',4,'Color','k');
m_contourf(lon,lat,Uplot','Color',[1,1,1],'LineWidth',1);
hold on
RuBU = cbrewer('seq','YlGnBu',24);
colormap(aa3,RuBU);
title(strcat('U',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR R^2'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([0 1])
ylabel(h1,'R^2')

aa3 = subplot(2,3,6);
padvar =4;
Vplot = padarray(R2V,[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',4,'Color','k');
m_contourf(lon,lat,Vplot','Color',[1,1,1],'LineWidth',1);
hold on
RuBU = cbrewer('seq','YlGnBu',24);
colormap(aa3,RuBU);
title(strcat('V',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR R^2'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([0 1])
ylabel(h1,'R^2')

set(ax11,'Position',[450         130        1972        1215])
% 
figsave=strcat('/Users/wchapman/Desktop/QUVfigs/',num2str(titMon(2:end)),'_',num2str(Level(RRa)),'_R2maps.jpg');
saveas(ax11,figsave);
close all 






%%
%%Try Two Nino 
VTperR = reshape(VTcatMemNino,size(VTcatMemNino,1),size(VTcatMemNino,2),size(VTcatMemNino,3)*size(VTcatMemNino,4));

%U R2 map
UperR = reshape(UcatMemNino,size(UcatMemNino,1),size(UcatMemNino,2),size(UcatMemNino,3)*size(UcatMemNino,4));

R2U = zeros(size(Uperiod,1),size(Uperiod,2));
for ii = 1:size(Uperiod,1)
    for jj = 1:size(Uperiod,2)
       
        
       Uts = UperR(ii,jj,:);
       VTts = VTperR(ii,jj,:);
       pU = polyfit(VTts,Uts,1);
       pUfit = pU(1)*VTts+pU(2);
       Uresid = Uts - pUfit;
       uSSresid = sum(Uresid.^2);
       uSStotal = (length(Uts)-1) * var(Uts);
       R2U(ii,jj) = 1 - uSSresid/uSStotal;      
        
    end
    
end

%V R2 map
VperR = reshape(VcatMemNino,size(VcatMemNino,1),size(VcatMemNino,2),size(VcatMemNino,3)*size(VcatMemNino,4));

R2V = zeros(size(Vperiod,1),size(Vperiod,2));
for ii = 1:size(Vperiod,1)
    for jj = 1:size(Vperiod,2)
       
        
       Vts = VperR(ii,jj,:);
       VTts = VTperR(ii,jj,:);
       pV = polyfit(VTts,Vts,1);
       pVfit = pV(1)*VTts+pV(2);
       Vresid = Vts - pVfit;
       VSSresid = sum(Vresid.^2);
       VSStotal = (length(Vts)-1) * var(Vts);
       R2V(ii,jj) = 1 - VSSresid/VSStotal;      
        
    end
    
end

%V R2 map
QperR = reshape(QcatMemNino,size(QcatMemNino,1),size(QcatMemNino,2),size(QcatMemNino,3)*size(QcatMemNino,4));

R2Q = zeros(size(Qperiod,1),size(Qperiod,2));
for ii = 1:size(Qperiod,1)
    for jj = 1:size(Qperiod,2)
       
        
       Qts = QperR(ii,jj,:);
       VTts = VTperR(ii,jj,:);
       pQ = polyfit(VTts,Qts,1);
       pQfit = pQ(1)*VTts+pQ(2);
       Qresid = Qts - pQfit;
       QSSresid = sum(Qresid.^2);
       QSStotal = (length(Qts)-1) * var(Qts);
       R2Q(ii,jj) = 1 - QSSresid/QSStotal;      
        
    end
    
end

VTts = reshape(VTcatMemNino,size(VTcatMemNino,1)*size(VTcatMemNino,2)*size(VTcatMemNino,3)*size(VTcatMemNino,4),1);
Qts = reshape(QcatMemNino,size(QcatMemNino,1)*size(QcatMemNino,2)*size(QcatMemNino,3)*size(QcatMemNino,4),1);
Uts = reshape(UcatMemNino,size(UcatMemNino,1)*size(UcatMemNino,2)*size(UcatMemNino,3)*size(UcatMemNino,4),1);
Vts = reshape(VcatMemNino,size(VcatMemNino,1)*size(VcatMemNino,2)*size(VcatMemNino,3)*size(VcatMemNino,4),1);
WSts = (Uts.^2+Vts.^2).^0.5;

pU = polyfit(VTts,Uts,1);
pV = polyfit(VTts,Vts,1);
pQ = polyfit(VTts,Qts,1);
pWS = polyfit(VTts,WSts,1);

pWSfit = pWS(1)*VTts+pWS(2);
WSresid = WSts - pWSfit;
wsSSresid = sum(WSresid.^2);
wsSStotal = (length(WSts)-1) * var(WSts);
WSrsq = 1 - wsSSresid/wsSStotal

pQfit = pQ(1)*VTts+pQ(2);
Qresid = Qts - pQfit;
qSSresid = sum(Qresid.^2);
qSStotal = (length(Qts)-1) * var(Qts);
Qrsq = 1 - qSSresid/qSStotal

pUfit = pU(1)*VTts+pU(2);
Uresid = Uts - pUfit;
uSSresid = sum(Uresid.^2);
uSStotal = (length(Uts)-1) * var(Uts);
Ursq = 1 - uSSresid/uSStotal

pVfit = pV(1)*VTts+pV(2);
Vresid = Vts - pVfit;
vSSresid = sum(Vresid.^2);
vSStotal = (length(Vts)-1) * var(Vts);
Vrsq = 1 - vSSresid/vSStotal


%%
%to make shit look pretty.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Q_ARmask700____19530101_19531231.mat');

%Pad LAT/LON
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

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
if ~make_it_tight,  clear subplot;  end

ax11 = figure;
subplot(2,3,2)
pUfit = pU(1)*VTts+pU(2);
plot((VTts),(Uts),'*')
xlabel('Vapor Transport [Kg/Kg * m/s]','FontSize',18);
ylabel('U wind [m/s]','FontSize',18);
hold on 
plot(VTts,pUfit,'LineWidth',3)
tit = strcat('Nino U',{'         '},'R^2:',{' '},num2str(Ursq));
title(tit,'FontSize',21)

subplot(2,3,3)
pVfit = pV(1)*VTts+pV(2);
plot(VTts,Vts,'*')
xlabel('Vapor Transport [Kg/Kg * m/s]','FontSize',18);
ylabel('V wind [m/s]','FontSize',18);
hold on 
plot(VTts,pVfit,'LineWidth',3)
tit = strcat('Nino V',{'         '},'R^2:',{' '},num2str(Vrsq));
title(tit,'FontSize',21)

subplot(2,3,1)
pQfit = pQ(1)*VTts+pQ(2);
plot(VTts,Qts,'*')
xlabel('Vapor Transport [Kg/Kg * m/s]','FontSize',18);
ylabel('Specific Humidity [Kg/Kg]','FontSize',18);
hold on 
plot(VTts,pQfit,'LineWidth',3)
tit = strcat('Nino Q',{'         '},'R^2:',{' '},num2str(Qrsq));
title(tit,'FontSize',21)

aa3 = subplot(2,3,4);
padvar =4;
Qplot = padarray(R2Q,[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',4,'Color','k');
m_contourf(lon,lat,Qplot','Color',[1,1,1],'LineWidth',1);
hold on
RuBU = cbrewer('seq','YlOrBr',24);
colormap(aa3,RuBU);
title(strcat('Nino Q',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR R^2'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
ylabel(h1,'R^2')


aa3 = subplot(2,3,5);
padvar =4;
Uplot = padarray(R2U,[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',4,'Color','k');
m_contourf(lon,lat,Uplot','Color',[1,1,1],'LineWidth',1);
hold on
RuBU = cbrewer('seq','YlGnBu',24);
colormap(aa3,RuBU);
title(strcat('Nino U',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR R^2'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([0 1])
ylabel(h1,'R^2')

aa3 = subplot(2,3,6);
padvar =4;
Vplot = padarray(R2V,[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',4,'Color','k');
m_contourf(lon,lat,Vplot','Color',[1,1,1],'LineWidth',1);
hold on
RuBU = cbrewer('seq','YlGnBu',24);
colormap(aa3,RuBU);
title(strcat('Nino V',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR R^2'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([0 1])
ylabel(h1,'R^2')

set(ax11,'Position',[450         130        1972        1215])

figsave=strcat('/Users/wchapman/Desktop/QUVfigs/',num2str(titMon(2:end)),'_',num2str(Level(RRa)),'_Nino_R2maps.jpg');
saveas(ax11,figsave);

close all 



%%
%%Try Two Nino 
VTperR = reshape(VTcatMemNina,size(VTcatMemNina,1),size(VTcatMemNina,2),size(VTcatMemNina,3)*size(VTcatMemNina,4));

%U R2 map
UperR = reshape(UcatMemNina,size(UcatMemNina,1),size(UcatMemNina,2),size(UcatMemNina,3)*size(UcatMemNina,4));

R2U = zeros(size(Uperiod,1),size(Uperiod,2));
for ii = 1:size(Uperiod,1)
    for jj = 1:size(Uperiod,2)
       
        
       Uts = UperR(ii,jj,:);
       VTts = VTperR(ii,jj,:);
       pU = polyfit(VTts,Uts,1);
       pUfit = pU(1)*VTts+pU(2);
       Uresid = Uts - pUfit;
       uSSresid = sum(Uresid.^2);
       uSStotal = (length(Uts)-1) * var(Uts);
       R2U(ii,jj) = 1 - uSSresid/uSStotal;      
        
    end
    
end

%V R2 map
VperR = reshape(VcatMemNina,size(VcatMemNina,1),size(VcatMemNina,2),size(VcatMemNina,3)*size(VcatMemNina,4));

R2V = zeros(size(Vperiod,1),size(Vperiod,2));
for ii = 1:size(Vperiod,1)
    for jj = 1:size(Vperiod,2)
       
        
       Vts = VperR(ii,jj,:);
       VTts = VTperR(ii,jj,:);
       pV = polyfit(VTts,Vts,1);
       pVfit = pV(1)*VTts+pV(2);
       Vresid = Vts - pVfit;
       VSSresid = sum(Vresid.^2);
       VSStotal = (length(Vts)-1) * var(Vts);
       R2V(ii,jj) = 1 - VSSresid/VSStotal;      
        
    end
    
end

%V R2 map
QperR = reshape(QcatMemNina,size(QcatMemNina,1),size(QcatMemNina,2),size(QcatMemNina,3)*size(QcatMemNina,4));

R2Q = zeros(size(Qperiod,1),size(Qperiod,2));
for ii = 1:size(Qperiod,1)
    for jj = 1:size(Qperiod,2)
       
        
       Qts = QperR(ii,jj,:);
       VTts = VTperR(ii,jj,:);
       pQ = polyfit(VTts,Qts,1);
       pQfit = pQ(1)*VTts+pQ(2);
       Qresid = Qts - pQfit;
       QSSresid = sum(Qresid.^2);
       QSStotal = (length(Qts)-1) * var(Qts);
       R2Q(ii,jj) = 1 - QSSresid/QSStotal;      
        
    end
    
end

VTts = reshape(VTcatMemNina,size(VTcatMemNina,1)*size(VTcatMemNina,2)*size(VTcatMemNina,3)*size(VTcatMemNina,4),1);
Qts = reshape(QcatMemNina,size(QcatMemNina,1)*size(QcatMemNina,2)*size(QcatMemNina,3)*size(QcatMemNina,4),1);
Uts = reshape(UcatMemNina,size(UcatMemNina,1)*size(UcatMemNina,2)*size(UcatMemNina,3)*size(UcatMemNina,4),1);
Vts = reshape(VcatMemNina,size(VcatMemNina,1)*size(VcatMemNina,2)*size(VcatMemNina,3)*size(VcatMemNina,4),1);
WSts = (Uts.^2+Vts.^2).^0.5;

pU = polyfit(VTts,Uts,1);
pV = polyfit(VTts,Vts,1);
pQ = polyfit(VTts,Qts,1);
pWS = polyfit(VTts,WSts,1);

pWSfit = pWS(1)*VTts+pWS(2);
WSresid = WSts - pWSfit;
wsSSresid = sum(WSresid.^2);
wsSStotal = (length(WSts)-1) * var(WSts);
WSrsq = 1 - wsSSresid/wsSStotal

pQfit = pQ(1)*VTts+pQ(2);
Qresid = Qts - pQfit;
qSSresid = sum(Qresid.^2);
qSStotal = (length(Qts)-1) * var(Qts);
Qrsq = 1 - qSSresid/qSStotal

pUfit = pU(1)*VTts+pU(2);
Uresid = Uts - pUfit;
uSSresid = sum(Uresid.^2);
uSStotal = (length(Uts)-1) * var(Uts);
Ursq = 1 - uSSresid/uSStotal

pVfit = pV(1)*VTts+pV(2);
Vresid = Vts - pVfit;
vSSresid = sum(Vresid.^2);
vSStotal = (length(Vts)-1) * var(Vts);
Vrsq = 1 - vSSresid/vSStotal


%%
%to make shit look pretty.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Q_ARmask700____19530101_19531231.mat');

%Pad LAT/LON
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

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
if ~make_it_tight,  clear subplot;  end

ax11 = figure;
subplot(2,3,2)
pUfit = pU(1)*VTts+pU(2);
plot((VTts),(Uts),'*')
xlabel('Vapor Transport [Kg/Kg * m/s]','FontSize',18);
ylabel('U wind [m/s]','FontSize',18);
hold on 
plot(VTts,pUfit,'LineWidth',3)
tit = strcat('Nina U',{'         '},'R^2:',{' '},num2str(Ursq));
title(tit,'FontSize',21)

subplot(2,3,3)
pVfit = pV(1)*VTts+pV(2);
plot(VTts,Vts,'*')
xlabel('Vapor Transport [Kg/Kg * m/s]','FontSize',18);
ylabel('V wind [m/s]','FontSize',18);
hold on 
plot(VTts,pVfit,'LineWidth',3)
tit = strcat('Nina V',{'         '},'R^2:',{' '},num2str(Vrsq));
title(tit,'FontSize',21)

subplot(2,3,1)
pQfit = pQ(1)*VTts+pQ(2);
plot(VTts,Qts,'*')
xlabel('Vapor Transport [Kg/Kg * m/s]','FontSize',18);
ylabel('Specific Humidity [Kg/Kg]','FontSize',18);
hold on 
plot(VTts,pQfit,'LineWidth',3)
tit = strcat('Nina Q',{'         '},'R^2:',{' '},num2str(Qrsq));
title(tit,'FontSize',21)

aa3 = subplot(2,3,4);
padvar =4;
Qplot = padarray(R2Q,[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',4,'Color','k');
m_contourf(lon,lat,Qplot','Color',[1,1,1],'LineWidth',1);
hold on
RuBU = cbrewer('seq','YlOrBr',24);
colormap(aa3,RuBU);
title(strcat('Nina Q',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR R^2'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
ylabel(h1,'R^2')


aa3 = subplot(2,3,5);
padvar =4;
Uplot = padarray(R2U,[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',4,'Color','k');
m_contourf(lon,lat,Uplot','Color',[1,1,1],'LineWidth',1);
hold on
RuBU = cbrewer('seq','YlGnBu',24);
colormap(aa3,RuBU);
title(strcat('Nina U',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR R^2'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([0 1])
ylabel(h1,'R^2')

aa3 = subplot(2,3,6);
padvar =4;
Vplot = padarray(R2V,[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',4,'Color','k');
m_contourf(lon,lat,Vplot','Color',[1,1,1],'LineWidth',1);
hold on
RuBU = cbrewer('seq','YlGnBu',24);
colormap(aa3,RuBU);
title(strcat('Nina V',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR R^2'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([0 1])
ylabel(h1,'R^2')

set(ax11,'Position',[450         130        1972        1215])

figsave=strcat('/Users/wchapman/Desktop/QUVfigs/',num2str(titMon(2:end)),'_',num2str(Level(RRa)),'_Nina_R2maps.jpg');
saveas(ax11,figsave);

close all 


%%
%%Try four noNino 
VTperR = reshape(VTcatMemNoNino,size(VTcatMemNino,1),size(VTcatMemNoNino,2),size(VTcatMemNoNino,3)*size(VTcatMemNoNino,4));

%U R2 map
UperR = reshape(UcatMemNoNino,size(UcatMemNoNino,1),size(UcatMemNoNino,2),size(UcatMemNoNino,3)*size(UcatMemNoNino,4));

R2U = zeros(size(Uperiod,1),size(Uperiod,2));
for ii = 1:size(Uperiod,1)
    for jj = 1:size(Uperiod,2)
       
        
       Uts = UperR(ii,jj,:);
       VTts = VTperR(ii,jj,:);
       pU = polyfit(VTts,Uts,1);
       pUfit = pU(1)*VTts+pU(2);
       Uresid = Uts - pUfit;
       uSSresid = sum(Uresid.^2);
       uSStotal = (length(Uts)-1) * var(Uts);
       R2U(ii,jj) = 1 - uSSresid/uSStotal;      
        
    end
    
end

%V R2 map
VperR = reshape(VcatMemNoNino,size(VcatMemNoNino,1),size(VcatMemNoNino,2),size(VcatMemNoNino,3)*size(VcatMemNoNino,4));

R2V = zeros(size(Vperiod,1),size(Vperiod,2));
for ii = 1:size(Vperiod,1)
    for jj = 1:size(Vperiod,2)
       
        
       Vts = VperR(ii,jj,:);
       VTts = VTperR(ii,jj,:);
       pV = polyfit(VTts,Vts,1);
       pVfit = pV(1)*VTts+pV(2);
       Vresid = Vts - pVfit;
       VSSresid = sum(Vresid.^2);
       VSStotal = (length(Vts)-1) * var(Vts);
       R2V(ii,jj) = 1 - VSSresid/VSStotal;      
        
    end
    
end

%V R2 map
QperR = reshape(QcatMemNoNino,size(QcatMemNoNino,1),size(QcatMemNoNino,2),size(QcatMemNoNino,3)*size(QcatMemNoNino,4));

R2Q = zeros(size(Qperiod,1),size(Qperiod,2));
for ii = 1:size(Qperiod,1)
    for jj = 1:size(Qperiod,2)
       
        
       Qts = QperR(ii,jj,:);
       VTts = VTperR(ii,jj,:);
       pQ = polyfit(VTts,Qts,1);
       pQfit = pQ(1)*VTts+pQ(2);
       Qresid = Qts - pQfit;
       QSSresid = sum(Qresid.^2);
       QSStotal = (length(Qts)-1) * var(Qts);
       R2Q(ii,jj) = 1 - QSSresid/QSStotal;      
        
    end
    
end

VTts = reshape(VTcatMemNoNino,size(VTcatMemNoNino,1)*size(VTcatMemNoNino,2)*size(VTcatMemNoNino,3)*size(VTcatMemNoNino,4),1);
Qts = reshape(QcatMemNoNino,size(QcatMemNoNino,1)*size(QcatMemNoNino,2)*size(QcatMemNoNino,3)*size(QcatMemNoNino,4),1);
Uts = reshape(UcatMemNoNino,size(UcatMemNoNino,1)*size(UcatMemNoNino,2)*size(UcatMemNoNino,3)*size(UcatMemNoNino,4),1);
Vts = reshape(VcatMemNoNino,size(VcatMemNoNino,1)*size(VcatMemNoNino,2)*size(VcatMemNoNino,3)*size(VcatMemNoNino,4),1);
WSts = (Uts.^2+Vts.^2).^0.5;

pU = polyfit((VTts),(Uts),1);
pV = polyfit(VTts,Vts,1);
pQ = polyfit(VTts,(Qts),1);
pWS = polyfit(VTts,WSts,1);

pWSfit = pWS(1)*VTts+pWS(2);
WSresid = WSts - pWSfit;
wsSSresid = sum(WSresid.^2);
wsSStotal = (length(WSts)-1) * var(WSts);
WSrsq = 1 - wsSSresid/wsSStotal

pQfit = pQ(1)*VTts+pQ(2);
Qresid = (Qts) - pQfit;
qSSresid = sum(Qresid.^2);
qSStotal = (length(Qts)-1) * var(Qts);
Qrsq = 1 - qSSresid/qSStotal

pUfit = pU(1)*(VTts)+pU(2);
Uresid = Uts - pUfit;
uSSresid = sum(Uresid.^2);
uSStotal = (length(Uts)-1) * var(Uts);
Ursq = 1 - uSSresid/uSStotal

pVfit = pV(1)*VTts+pV(2);
Vresid = Vts - pVfit;
vSSresid = sum(Vresid.^2);
vSStotal = (length(Vts)-1) * var(Vts);
Vrsq = 1 - vSSresid/vSStotal


%%
%to make shit look pretty.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Q_ARmask700____19530101_19531231.mat');

%Pad LAT/LON
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

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
if ~make_it_tight,  clear subplot;  end

ax11 = figure;
subplot(2,3,2)
pUfit = pU(1)*VTts+pU(2);
plot((VTts),(Uts),'*')
xlabel('Vapor Transport [Kg/Kg * m/s]','FontSize',18);
ylabel('U wind [m/s]','FontSize',18);
hold on 
plot(VTts,pUfit,'LineWidth',3)
tit = strcat('Nuetral U',{'         '},'R^2:',{' '},num2str(Ursq));
title(tit,'FontSize',21)

subplot(2,3,3)
pVfit = pV(1)*VTts+pV(2);
plot(VTts,Vts,'*')
xlabel('Vapor Transport [Kg/Kg * m/s]','FontSize',18);
ylabel('V wind [m/s]','FontSize',18);
hold on 
plot(VTts,pVfit,'LineWidth',3)
tit = strcat('Nuetral V',{'         '},'R^2:',{' '},num2str(Vrsq));
title(tit,'FontSize',21)

subplot(2,3,1)
pQfit = pQ(1)*VTts+pQ(2);
plot(VTts,Qts,'*')
xlabel('Vapor Transport [Kg/Kg * m/s]','FontSize',18);
ylabel('Specific Humidity [Kg/Kg]','FontSize',18);
hold on 
plot(VTts,pQfit,'LineWidth',3)
tit = strcat('Nuetral Q',{'         '},'R^2:',{' '},num2str(Qrsq));
title(tit,'FontSize',21)

aa3 = subplot(2,3,4);
padvar =4;
Qplot = padarray(R2Q,[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',4,'Color','k');
m_contourf(lon,lat,Qplot','Color',[1,1,1],'LineWidth',1);
hold on
RuBU = cbrewer('seq','YlOrBr',24);
colormap(aa3,RuBU);
title(strcat('Nuetral Q',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR R^2'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
ylabel(h1,'R^2')


aa3 = subplot(2,3,5);
padvar =4;
Uplot = padarray(R2U,[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',4,'Color','k');
m_contourf(lon,lat,Uplot','Color',[1,1,1],'LineWidth',1);
hold on
RuBU = cbrewer('seq','YlGnBu',24);
colormap(aa3,RuBU);
title(strcat('Nuetral U',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR R^2'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([0 1])
ylabel(h1,'R^2')

aa3 = subplot(2,3,6);
padvar =4;
Vplot = padarray(R2V,[padvar,padvar],nan,'both');
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',4,'Color','k');
m_contourf(lon,lat,Vplot','Color',[1,1,1],'LineWidth',1);
hold on
RuBU = cbrewer('seq','YlGnBu',24);
colormap(aa3,RuBU);
title(strcat('Nuetral V',{' '},titMon(2:end),{' '},num2str(Level(RRa)),'hPa AR R^2'))
ax = gca;
ax.FontSize = 21;
h1  = colorbar;
caxis([0 1])
ylabel(h1,'R^2')

set(ax11,'Position',[450         130        1972        1215])

figsave=strcat('/Users/wchapman/Desktop/QUVfigs/',num2str(titMon(2:end)),'_',num2str(Level(RRa)),'_Nuetral_R2maps.jpg');
saveas(ax11,figsave);

close all 

end






