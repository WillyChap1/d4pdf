%Vert Cross section. 
clear all
close all

%based on 925 DJF Variance Plot. 
%Mean lon: 150W, 25-40N
%Nina lon: 150W  25-40N
%Nino lon: 148W 25- 40N
%Neutral lon:150W 25-40N 

yearsinrange = 1952:1:2010;

desLat = [22,42];
desLon = 210;

% =========================================================================
% ============================ Switches ===================================
% =========================================================================

Level = [1000,925,850,700,600,500,400,300]; %hPa levels
Ninomonths1 = [12];  %months in first year el nino [Dec]
Ninomonths2 = [1 2];  %months in second year el nino [JAN, FEB] +1 -

%months to test for map.
months1 = [12];  %months in first year ex. [NOV, DEC]  0
months2 = [1 2];  %months in second year ex. [JAN, FEB] +1 -

anomaNINO = [0.5,4];    %min and max anomaly nino
anomaNINA = [0.5,4];    %min and max anomaly nina 

daysinMonth = [31 28 31 30 31 30 31 31 30 31 30 31];
% =========================================================================

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
plot(NoNinoyears,anomnino(anomnino<=anomaNINO(1) & anomnino>=-anomaNINO(1)),'r*')

cd /Users/wchapman/SIO/research/100Amip/d4pdf/Q/ARmasked 

load('Q_ARmask1000____19520101_19521231.mat')

[~,latST] = (min(abs(desLat(1)-lat)));
[~,latEN] = (min(abs(desLat(2)-lat)));
[~,lonDO] = (min(abs(desLon-lon)));

QperiodLev=[];
UperiodLev=[];
VperiodLev=[];
% 
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

%% Let's get our data. 
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
        
        Qtemp = Qmask(lonDO,latST:latEN,:,indexdo);
        Utemp = Umask(lonDO,latST:latEN,:,indexdo);
        Vtemp = Vmask(lonDO,latST:latEN,:,indexdo);
        
        %QkeepMonth = cat(4,QkeepMonth,sum(Qtemp,4));
        %QkeepAll = cat(4,QkeepAll,Qtemp);
        %timekeepMonths = cat(1,timekeepMonths,timetot(indexdo(1)));
        %timekeep = cat(1,timekeep,timetot(indexdo));
        
        Qperiodtemp = cat(4,Qperiodtemp,Qtemp);
        Uperiodtemp = cat(4,Uperiodtemp,Utemp);
        Vperiodtemp = cat(4,Vperiodtemp,Vtemp);
        
        
        timePeriod = cat(1,timePeriod,timetot(indexdo(1)));
        
        
    end

    end
    
    if isempty(months2)==0
  
    fileQ = pp(hh+1).name;
    disp(fileQ)
    tic
    load(fileQ);
    fileU = strcat('/Users/wchapman/SIO/research/100Amip/d4pdf/UV/ARmasked/',qq(hh+1).name);
    disp(fileU)
    load(fileU);
    fileV = strcat('/Users/wchapman/SIO/research/100Amip/d4pdf/UV/ARmasked/',rr(hh+1).name);
    load(fileV);
    
    
    timetot = date1+minutes(time);
        
    for jr = min(months2):max(months2)    %Dec Nov. 
        [~,month,~]=datevec(timetot);
        indexdo = find(month==jr);
        Qtemp = Qmask(lonDO,latST:latEN,:,indexdo);
        Utemp = Umask(lonDO,latST:latEN,:,indexdo);
        Vtemp = Vmask(lonDO,latST:latEN,:,indexdo);
        
        %QkeepMonth = cat(4,QkeepMonth,sum(Qtemp,4)); 
        %QkeepAll = cat(4,QkeepAll,Qtemp);
        %timekeepMonths = cat(1,timekeepMonths,timetot(indexdo(1)));
        %timekeep = cat(1,timekeep,timetot(indexdo));
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
   
    QperiodLev = cat(5,QperiodLev,Qperiod);
    UperiodLev = cat(5,UperiodLev,Uperiod);
    VperiodLev = cat(5,VperiodLev,Vperiod);


end

%load('/Users/wchapman/SIO/research/100Amip/datasets/DJF_QUVperiodLev.mat')
save('/Users/wchapman/SIO/research/100Amip/datasets/DJF_QUVperiodLev.mat','QperiodLev','timekeepPeriod','UperiodLev','Vperiodlev');
VTperiodLev = ((QperiodLev.*UperiodLev).^2+(QperiodLev.*VperiodLev).^2).^0.5;

%%  Divde into Nino, Nina, and Nuetral 
clear UNino UNina UNoNino VNino VNina VNoNino QNino QNina QNoNino VTNino VTNina VTNoNino
%==========================================================================
%========================       Nino .       ==============================
%==========================================================================

[yearsIn,~] = datevec(timekeepPeriod);

for ff = 1:length(Ninoyears)
    
    NinoIndex = find(yearsIn==Ninoyears(ff));
    
    QNino(:,:,:,ff) = QperiodLev(:,:,:,NinoIndex,:); 
    UNino(:,:,:,ff) = UperiodLev(:,:,:,NinoIndex,:);
    VNino(:,:,:,ff) = VperiodLev(:,:,:,NinoIndex,:);
    VTNino(:,:,:,ff) = VTperiodLev(:,:,:,NinoIndex,:);

    timeNino(ff) = timekeepPeriod(NinoIndex);
    
end
QNino = permute(QNino,[5,1,2,4,3]);
UNino = permute(UNino,[5,1,2,4,3]);
VNino = permute(VNino,[5,1,2,4,3]);
VTNino = permute(VTNino,[5,1,2,4,3]);

%==========================================================================
%========================       Nina .       ==============================
%==========================================================================

[yearsIn,~] = datevec(timekeepPeriod);

for ff = 1:length(Ninayears)
    
    NinaIndex = find(yearsIn==Ninayears(ff));
    
    QNina(:,:,:,ff) = QperiodLev(:,:,:,NinaIndex,:); 
    UNina(:,:,:,ff) = UperiodLev(:,:,:,NinaIndex,:);
    VNina(:,:,:,ff) = VperiodLev(:,:,:,NinaIndex,:);
    VTNina(:,:,:,ff) = VTperiodLev(:,:,:,NinaIndex,:);

    timeNina(ff) = timekeepPeriod(NinaIndex);
    
end
QNina = permute(QNina,[5,1,2,4,3]);
UNina = permute(UNina,[5,1,2,4,3]);
VNina = permute(VNina,[5,1,2,4,3]);
VTNina = permute(VTNina,[5,1,2,4,3]);

%==========================================================================
%========================       Neutral .       ==============================
%==========================================================================

[yearsIn,~] = datevec(timekeepPeriod);

for ff = 1:length(NoNinoyears)
    
    NoNinoIndex = find(yearsIn==NoNinoyears(ff));
    
    QNoNino(:,:,:,ff,:) = QperiodLev(:,:,:,NoNinoIndex,:); 
    UNoNino(:,:,:,ff,:) = UperiodLev(:,:,:,NoNinoIndex,:);
    VNoNino(:,:,:,ff,:) = VperiodLev(:,:,:,NoNinoIndex,:);
    VTNoNino(:,:,:,ff,:) = VTperiodLev(:,:,:,NoNinoIndex,:);

    timeNoNino(ff) = timekeepPeriod(NoNinoIndex);
    
end
QNoNino = permute(QNoNino,[5,1,2,4,3]);
UNoNino = permute(UNoNino,[5,1,2,4,3]);
VNoNino = permute(VNoNino,[5,1,2,4,3]);
VTNoNino = permute(VTNoNino,[5,1,2,4,3]);

%%
%mean it all.
QmeanLev = mean(QperiodLev,4);
QmeanLev = mean(QmeanLev,3);

VmeanLev = mean(VperiodLev,4);
VmeanLev = mean(VmeanLev,3);

VTmeanLev = mean(VTperiodLev,4);
VTmeanLev = mean(VTmeanLev,3);

UmeanLev = mean(UperiodLev,4);
UmeanLev = mean(UmeanLev,3);

QmeanLevNino = mean(QNino,4);
QmeanLevNino = mean(QmeanLevNino,3);

VmeanLevNino = mean(VNino,4);
VmeanLevNino = mean(VmeanLevNino,3);

VTmeanLevNino = mean(VTNino,4);
VTmeanLevNino = mean(VTmeanLevNino,3);

UmeanLevNino = mean(UNino,4);
UmeanLevNino = mean(UmeanLevNino,3);

QmeanLevNina = mean(QNina,4);
QmeanLevNina = mean(QmeanLevNina,3);

VmeanLevNina = mean(VNina,4);
VmeanLevNina = mean(VmeanLevNina,3);

VTmeanLevNina = mean(VTNina,4);
VTmeanLevNina = mean(VTmeanLevNina,3);

UmeanLevNina = mean(UNina,4);
UmeanLevNina = mean(UmeanLevNina,3);

QmeanLevNoNino = mean(QNoNino,4);
QmeanLevNoNino = mean(QmeanLevNoNino,3);

VmeanLevNoNino = mean(VNoNino,4);
VmeanLevNoNino = mean(VmeanLevNoNino,3);

VTmeanLevNoNino = mean(VTNoNino,4);
VTmeanLevNoNino = mean(VTmeanLevNoNino,3);

UmeanLevNoNino = mean(UNoNino,4);
UmeanLevNoNino = mean(UmeanLevNoNino,3);


%variance it all
Qanom = QperiodLev - QmeanLev;
Qint = std(Qanom,0,3);
Qint = mean(Qint,4);
Qint = squeeze(Qint);
Qint=Qint';
Qint = flipud(Qint);

Uanom = UperiodLev - UmeanLev;
Uint = std(Uanom,0,3);
Uint = mean(Uint,4);
Uint = squeeze(Uint);
Uint=Uint';
Uint = flipud(Uint);

Vanom = VperiodLev - VmeanLev;
Vint = std(Vanom,0,3);
Vint = mean(Vint,4);
Vint = squeeze(Vint);
Vint=Vint';
Vint = flipud(Vint);

VTanom = VTperiodLev - VTmeanLev;
VTint = std(VTanom,0,3);
VTint = mean(VTint,4);
VTint = squeeze(VTint);
VTint=VTint';
VTint = flipud(VTint);

QanomNino = QNino - QmeanLevNino;
QintNino = std(QanomNino,0,3);
QintNino = mean(QintNino,4);
QintNino = squeeze(QintNino);
QintNino=QintNino';
QintNino = flipud(QintNino);

UanomNino = UNino - UmeanLevNino;
UintNino = std(UanomNino,0,3);
UintNino = mean(UintNino,4);
UintNino = squeeze(UintNino);
UintNino=UintNino';
UintNino = flipud(UintNino);

VanomNino = VNino - VmeanLevNino;
VintNino = std(VanomNino,0,3);
VintNino = mean(VintNino,4);
VintNino = squeeze(VintNino);
VintNino=VintNino';
VintNino = flipud(VintNino);

VTanomNino = VTNino - VTmeanLevNino;
VTintNino = std(VTanomNino,0,3);
VTintNino = mean(VTintNino,4);
VTintNino = squeeze(VTintNino);
VTintNino=VTintNino';
VTintNino = flipud(VTintNino);

QanomNina = QNina - QmeanLevNina;
QintNina = std(QanomNina,0,3);
QintNina = mean(QintNina,4);
QintNina = squeeze(QintNina);
QintNina=QintNina';
QintNina = flipud(QintNina);

UanomNina = UNina - UmeanLevNina;
UintNina = std(UanomNina,0,3);
UintNina = mean(UintNina,4);
UintNina = squeeze(UintNina);
UintNina=UintNina';
UintNina = flipud(UintNina);

VanomNina = VNina - VmeanLevNina;
VintNina = std(VanomNina,0,3);
VintNina = mean(VintNina,4);
VintNina = squeeze(VintNina);
VintNina=VintNina';
VintNina = flipud(VintNina);

VTanomNina = VTNina - VTmeanLevNina;
VTintNina = std(VTanomNina,0,3);
VTintNina = mean(VTintNina,4);
VTintNina = squeeze(VTintNina);
VTintNina=VTintNina';
VTintNina = flipud(VTintNina);

QanomNoNino = QNoNino - QmeanLevNoNino;
QintNoNino = std(QanomNoNino,0,3);
QintNoNino = mean(QintNoNino,4);
QintNoNino = squeeze(QintNoNino);
QintNoNino=QintNoNino';
QintNoNino = flipud(QintNoNino);

UanomNoNino = UNoNino - UmeanLevNoNino;
UintNoNino = std(UanomNoNino,0,3);
UintNoNino = mean(UintNoNino,4);
UintNoNino = squeeze(UintNoNino);
UintNoNino=UintNoNino';
UintNoNino = flipud(UintNoNino);

VanomNoNino = VNoNino - VmeanLevNoNino;
VintNoNino = std(VanomNoNino,0,3);
VintNoNino = mean(VintNoNino,4);
VintNoNino = squeeze(VintNoNino);
VintNoNino=VintNoNino';
VintNoNino = flipud(VintNoNino);

VTanomNoNino = VTNoNino - VTmeanLevNoNino;
VTintNoNino = std(VTanomNoNino,0,3);
VTintNoNino = mean(VTintNoNino,4);
VTintNoNino = squeeze(VTintNoNino);
VTintNoNino=VTintNoNino';
VTintNoNino = flipud(VTintNoNino);

QmeanLev = squeeze(QmeanLev);
QmeanLev = flipud(QmeanLev');
VmeanLev = squeeze(VmeanLev);
VmeanLev = flipud(VmeanLev');
VTmeanLev = squeeze(VTmeanLev);
VTmeanLev = flipud(VTmeanLev');
UmeanLev = squeeze(UmeanLev);
UmeanLev = flipud(UmeanLev');

QmeanLevNino = squeeze(QmeanLevNino);
QmeanLevNino = flipud(QmeanLevNino');
VmeanLevNino = squeeze(VmeanLevNino);
VmeanLevNino = flipud(VmeanLevNino');
VTmeanLevNino = squeeze(VTmeanLevNino);
VTmeanLevNino = flipud(VTmeanLevNino');
UmeanLevNino = squeeze(UmeanLevNino);
UmeanLevNino = flipud(UmeanLevNino');


QmeanLevNina = squeeze(QmeanLevNina);
QmeanLevNina = flipud(QmeanLevNina');
VmeanLevNina = squeeze(VmeanLevNina);
VmeanLevNina = flipud(VmeanLevNina');
VTmeanLevNina = squeeze(VTmeanLevNina);
VTmeanLevNina = flipud(VTmeanLevNina');
UmeanLevNina = squeeze(UmeanLevNina);
UmeanLevNina = flipud(UmeanLevNina');


QmeanLevNoNino = squeeze(QmeanLevNoNino);
QmeanLevNoNino = flipud(QmeanLevNoNino');
VmeanLevNoNino = squeeze(VmeanLevNoNino);
VmeanLevNoNino = flipud(VmeanLevNoNino');
VTmeanLevNoNino = squeeze(VTmeanLevNoNino);
VTmeanLevNoNino = flipud(VTmeanLevNoNino');
UmeanLevNoNino = squeeze(UmeanLevNoNino);
UmeanLevNoNino = flipud(UmeanLevNoNino');


%% plotting!  

levs=[300,400,500,600,700,850,925,1000];
lats = lat(latST:latEN);

%options: [gap]   %hght [Top Bot]  %Wid  [L R]                          
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
if ~make_it_tight,  clear subplot;  end

load('/Users/wchapman/SIO/research/100Amip/d4pdf/Q/ARmasked/Q_ARmask300____19520101_19521231.mat');
titMon = [];
mons =['J','F','M','A','M','J','J','A','S','O','N','D'];
for yy = 1:length(monthsall)
    
    mon = monthsall(yy);
    
    monC = mons(mon);
    
    titMon = strcat(titMon,monC);      
end 
ax11=figure;

aa2=subplot(4,5,6);
padvar =4;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
m_line([desLon desLon],[desLat(1) desLat(2)],'LineWidth',5,'Color','r')
hold on
title(strcat('R.O.I.'))
ax = gca;
ax.FontSize = 25;

aa2=subplot(4,5,2);
contourf(lats,levs,Qint);
ggg = cbrewer('seq','YlGnBu',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'Q Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
ylabel(yy,'Kg/Kg')

caxis([0 0.25])

aa2=subplot(4,5,3);
contourf(lats,levs,QintNino);
ggg = cbrewer('seq','YlGnBu',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'Q- Nino Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
caxis([0 0.25])
ylabel(yy,'Kg/Kg')



aa2=subplot(4,5,4);
contourf(lats,levs,QintNina);
ggg = cbrewer('seq','YlGnBu',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'Q- Nina Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
caxis([0 0.25])
ylabel(yy,'Kg/Kg')


aa2=subplot(4,5,5);
contourf(lats,levs,QintNoNino);
ggg = cbrewer('seq','YlGnBu',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'Q- Neutral Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
caxis([0 0.25])
ylabel(yy,'Kg/Kg')


aa2=subplot(4,5,7);
contourf(lats,levs,Uint);
ggg = cbrewer('seq','YlOrRd',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'U- Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
caxis([100 1200])
ylabel(yy,'m/s')



aa2=subplot(4,5,8);
contourf(lats,levs,UintNino);
ggg = cbrewer('seq','YlOrRd',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'U- Nino - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar;
caxis([100 1200])
ylabel(yy,'m/s')


aa2=subplot(4,5,9);
contourf(lats,levs,UintNina);
ggg = cbrewer('seq','YlOrRd',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'U- Nina - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar;
caxis([100 1200])
ylabel(yy,'m/s')


aa2=subplot(4,5,10);
contourf(lats,levs,UintNoNino);
ggg = cbrewer('seq','YlOrRd',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'U- Neutral - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar;
caxis([100 1200])
ylabel(yy,'m/s')


aa2=subplot(4,5,12);
contourf(lats,levs,Vint);
ggg = cbrewer('seq','PuBuGn',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'V - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
caxis([50 550])
ylabel(yy,'m/s')


aa2=subplot(4,5,13);
contourf(lats,levs,VintNino);
ggg = cbrewer('seq','PuBuGn',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'V- Nino - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
caxis([50 550])
ylabel(yy,'m/s')


aa2=subplot(4,5,14);
contourf(lats,levs,VintNina);
ggg = cbrewer('seq','PuBuGn',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'V- Nina - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
caxis([50 550])
ylabel(yy,'m/s')


aa2=subplot(4,5,15);
contourf(lats,levs,VintNoNino);
ggg = cbrewer('seq','PuBuGn',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'V- Neutral - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar;
caxis([50 550])
ylabel(yy,'m/s')


aa2=subplot(4,5,17);
contourf(lats,levs,VTint);
ggg = cbrewer('seq','GnBu',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'VT - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
caxis([0 600])
ylabel(yy,'Kg/(ms)')


aa2=subplot(4,5,18);
contourf(lats,levs,VTintNino);
ggg = cbrewer('seq','GnBu',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'VT- Nino - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
caxis([0 600])
ylabel(yy,'Kg/(ms)')


aa2=subplot(4,5,19);
contourf(lats,levs,VTintNina);
ggg = cbrewer('seq','GnBu',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'VT- Nina - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
caxis([0 600])
ylabel(yy,'Kg/(ms)')


aa2=subplot(4,5,20);
contourf(lats,levs,VTintNoNino);
ggg = cbrewer('seq','GnBu',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'VT- Neutral - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar;
caxis([0 600])
ylabel(yy,'Kg/(ms)')


set(ax11,'Position',[  16          82        2515        1263]);

%%
ax12=figure;

aa2=subplot(4,5,6);
padvar =4;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
m_line([desLon desLon],[desLat(1) desLat(2)],'LineWidth',5,'Color','r')
hold on
title(strcat('R.O.I.'))
ax = gca;
ax.FontSize = 25;

aa2=subplot(4,5,2);
contourf(lats,levs,Qint);
ggg = cbrewer('seq','YlGnBu',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'Q - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
ylabel(yy,'Kg/Kg')


aa2=subplot(4,5,3);
contourf(lats,levs,QintNino-Qint);
ggg = cbrewer('div','RdYlBu',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'Q- Nino - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
ylabel(yy,'Kg/Kg')
caxis([-.025,.025])


aa2=subplot(4,5,4);
contourf(lats,levs,QintNina-Qint);
ggg = cbrewer('div','RdYlBu',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'Q- Nina - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
ylabel(yy,'Kg/Kg')
caxis([-.025,.025])


aa2=subplot(4,5,5);
contourf(lats,levs,QintNoNino-Qint);
ggg = cbrewer('div','RdYlBu',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'Q- Neutral - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
ylabel(yy,'Kg/Kg')
caxis([-.025,.025])


aa2=subplot(4,5,7);
contourf(lats,levs,Uint);
ggg = cbrewer('seq','YlOrRd',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'U - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
caxis([100 1200])
ylabel(yy,'m/s')



aa2=subplot(4,5,8);
contourf(lats,levs,UintNino-Uint);
ggg = cbrewer('div','RdGy',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'U- Nino - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar;
ylabel(yy,'m/s')
caxis([-300 300])


aa2=subplot(4,5,9);
contourf(lats,levs,UintNina-Uint);
ggg = cbrewer('div','RdGy',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'U- Nina - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar;
ylabel(yy,'m/s')
caxis([-300 300])


aa2=subplot(4,5,10);
contourf(lats,levs,UintNoNino-Uint);
ggg = cbrewer('div','RdGy',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'U - Neutral - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar;
ylabel(yy,'m/s')
caxis([-300 300])


aa2=subplot(4,5,12);
contourf(lats,levs,Vint);
ggg = cbrewer('seq','PuBuGn',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'V - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
ylabel(yy,'m/s')
caxis([50 550])


aa2=subplot(4,5,13);
contourf(lats,levs,VintNino-Vint);
ggg = cbrewer('div','PRGn',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'V- Nino - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
ylabel(yy,'m/s')
caxis([-50 50])


aa2=subplot(4,5,14);
contourf(lats,levs,VintNina-Vint);
ggg = cbrewer('div','PRGn',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'V- Nina - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
ylabel(yy,'m/s')
caxis([-50 50])


aa2=subplot(4,5,15);
contourf(lats,levs,VintNoNino-Vint);
ggg = cbrewer('div','PRGn',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'V- Neutral - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar;
ylabel(yy,'m/s')
caxis([-50 50])


aa2=subplot(4,5,17);
contourf(lats,levs,VTint);
ggg = cbrewer('seq','GnBu',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'VT - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
caxis([0 600])
ylabel(yy,'Kg/(ms)')


aa2=subplot(4,5,18);
contourf(lats,levs,VTintNino-VTint);
ggg = cbrewer('div','BrBG',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'VT- Nino - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
ylabel(yy,'Kg/(ms)')
caxis([-90 90])


aa2=subplot(4,5,19);
contourf(lats,levs,VTintNina-VTint);
ggg = cbrewer('div','BrBG',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'VT- Nina - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
ylabel(yy,'Kg/(ms)')
caxis([-90 90])


aa2=subplot(4,5,20);
contourf(lats,levs,VTintNoNino-VTint);
ggg = cbrewer('div','BrBG',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'VT- Neutral - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar;
ylabel(yy,'Kg/(ms)')
caxis([-90 90])


set(ax12,'Position',[  16          82        2515        1263]);




annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Anomaly', ...
    'EdgeColor', 'none', ...
    'Position', [0.665    0.9000    1.0000    0.1000],'FontSize',25,'FontWeight','bold');

annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Mean', ...
    'EdgeColor', 'none', ...
    'Position', [0.285   0.9000    1.0000    0.1000],'FontSize',25,'FontWeight','bold');
%vert line 
annotation('line',[0.40 0.40],[0.978 0.02],'color','k');
annotation('line',[0.395 0.395],[0.978 0.02],'color','k');
annotation('line',[0.22 0.395],[0.978 0.978],'color','k');
annotation('line',[0.4 0.97],[0.978 0.978],'color','k');

%%
ax13=figure;

aa2=subplot(4,5,6);
padvar =4;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k')
m_line([desLon desLon],[desLat(1) desLat(2)],'LineWidth',5,'Color','r')
hold on
title(strcat('R.O.I.'))
ax = gca;
ax.FontSize = 25;

aa2=subplot(4,5,2);
contourf(lats,levs,Qint);
ggg = cbrewer('seq','YlGnBu',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'Q - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
ylabel(yy,'Kg/Kg')


aa2=subplot(4,5,3);
contourf(lats,levs,QintNino-Qint);
ggg = cbrewer('div','RdYlBu',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'Q- Nino - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
ylabel(yy,'Kg/Kg')
caxis([-.025,.025])


aa2=subplot(4,5,4);
contourf(lats,levs,QintNina-Qint);
ggg = cbrewer('div','RdYlBu',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'Q- Nina - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
ylabel(yy,'Kg/Kg')
caxis([-.025,.025])


aa2=subplot(4,5,5);
contourf(lats,levs,QintNoNino-Qint);
ggg = cbrewer('div','RdYlBu',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'Q- Neutral - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
ylabel(yy,'Kg/Kg')
caxis([-.025,.025])


aa2=subplot(4,5,7);
contourf(lats,levs,Uint);
ggg = cbrewer('seq','YlOrRd',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'U - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
caxis([100 1200])
ylabel(yy,'m/s')



aa2=subplot(4,5,8);
contourf(lats,levs,UintNino-Uint);
ggg = cbrewer('div','RdGy',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'U- Nino - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar;
ylabel(yy,'m/s')
caxis([-300 300])


aa2=subplot(4,5,9);
contourf(lats,levs,UintNina-Uint);
ggg = cbrewer('div','RdGy',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'U- Nina - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar;
ylabel(yy,'m/s')
caxis([-300 300])


aa2=subplot(4,5,10);
contourf(lats,levs,UintNoNino-Uint);
ggg = cbrewer('div','RdGy',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'U - Neutral - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar;
ylabel(yy,'m/s')
caxis([-300 300])


aa2=subplot(4,5,12);
contourf(lats,levs,Vint);
ggg = cbrewer('seq','PuBuGn',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'V - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
ylabel(yy,'m/s')
caxis([50 550])


aa2=subplot(4,5,13);
contourf(lats,levs,VintNino-Vint);
ggg = cbrewer('div','PRGn',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'V- Nino - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
ylabel(yy,'m/s')
caxis([-50 50])


aa2=subplot(4,5,14);
contourf(lats,levs,VintNina-Vint);
ggg = cbrewer('div','PRGn',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'V- Nina - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
ylabel(yy,'m/s')
caxis([-50 50])


aa2=subplot(4,5,15);
contourf(lats,levs,VintNoNino-Vint);
ggg = cbrewer('div','PRGn',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'V- Neutral - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar;
ylabel(yy,'m/s')
caxis([-50 50])


aa2=subplot(4,5,17);
contourf(lats,levs,VTint);
ggg = cbrewer('seq','GnBu',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'VT - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
caxis([0 600])
ylabel(yy,'Kg/(ms)')


aa2=subplot(4,5,18);
contourf(lats,levs,(VTintNino-VTint)./VTint);
ggg = cbrewer('div','BrBG',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'VT- Nino - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
ylabel(yy,'Kg/(ms)')
caxis([-90 90])


aa2=subplot(4,5,19);
contourf(lats,levs,VTintNina-VTint);
ggg = cbrewer('div','BrBG',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'VT- Nina - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar; 
ylabel(yy,'Kg/(ms)')
caxis([-90 90])


aa2=subplot(4,5,20);
contourf(lats,levs,VTintNoNino-VTint);
ggg = cbrewer('div','BrBG',24);
colormap(aa2,ggg);
set(gca,'YDir','reverse');
title(strcat(titMon,{' '},'VT- Neutral - Internal Variability'));
ylabel('Pressure [hPa]');
xlabel('Latitude [N^{o}]');
yy=colorbar;
ylabel(yy,'Kg/(ms)')
caxis([-90 90])


set(ax13,'Position',[  16          82        2515        1263]);


annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Standard Anomaly', ...
    'EdgeColor', 'none', ...
    'Position', [0.643    0.9000    1.0000    0.1000],'FontSize',25,'FontWeight','bold');

annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Mean', ...
    'EdgeColor', 'none', ...
    'Position', [0.285   0.9000    1.0000    0.1000],'FontSize',25,'FontWeight','bold');
%vert line 
annotation('line',[0.40 0.40],[0.978 0.02],'color','k');
annotation('line',[0.395 0.395],[0.978 0.02],'color','k');
annotation('line',[0.22 0.395],[0.978 0.978],'color','k');
annotation('line',[0.4 0.97],[0.978 0.978],'color','k');


