%make weird spaggetti plot
close all
load('/Users/wchapman/SIO/research/100Amip/datasets/DJF_QUVperiodLev.mat')
levs=[300,400,500,600,700,850,925,1000];
VTperiodLev = ((QperiodLev.*UperiodLev).^2+(QperiodLev.*VperiodLev).^2).^0.5;


%% erase this later TEMP!!!!

datestart = '1952-12-31 00:00';
formatin = 'yyyy-MM-dd HH:mm';
date1 = datetime(datestart,'InputFormat',formatin);
times = 0:56;
timekeepPeriod = date1+years(times);



%%  Divde into Nino, Nina, and Nuetral 
clear UNino UNina UNoNino VNino VNina VNoNino QNino QNina QNoNino VTNino VTNina VTNoNino
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

Qmean = mean(QperiodLev,3);
Qmean = mean(Qmean,2);
Qmean =mean(Qmean,4);
Qmean = squeeze(Qmean);

Umean = mean(UperiodLev,3);
Umean = mean(Umean,2);
Umean =mean(Umean,4);
Umean = squeeze(Umean);

Vmean = mean(VperiodLev,3);
Vmean = mean(Vmean,2);
Vmean =mean(Vmean,4);
Vmean = squeeze(Vmean);

VTmean = mean(VTperiodLev,3);
VTmean = mean(VTmean,2);
VTmean =mean(VTmean,4);
VTmean = squeeze(VTmean);

QmeanNino = mean(QNino,3);
QmeanNino = mean(QmeanNino,2);
QmeanNino =mean(QmeanNino,4);
QmeanNino = squeeze(QmeanNino);

UmeanNino = mean(UNino,3);
UmeanNino = mean(UmeanNino,2);
UmeanNino =mean(UmeanNino,4);
UmeanNino = squeeze(UmeanNino);

VmeanNino = mean(VNino,3);
VmeanNino = mean(VmeanNino,2);
VmeanNino =mean(VmeanNino,4);
VmeanNino = squeeze(VmeanNino);

VTmeanNino = mean(VTNino,3);
VTmeanNino = mean(VTmeanNino,2);
VTmeanNino =mean(VTmeanNino,4);
VTmeanNino = squeeze(VTmeanNino);

QmeanNina = mean(QNina,3);
QmeanNina = mean(QmeanNina,2);
QmeanNina =mean(QmeanNina,4);
QmeanNina = squeeze(QmeanNina);

UmeanNina = mean(UNina,3);
UmeanNina = mean(UmeanNina,2);
UmeanNina =mean(UmeanNina,4);
UmeanNina = squeeze(UmeanNina);

VmeanNina = mean(VNina,3);
VmeanNina = mean(VmeanNina,2);
VmeanNina =mean(VmeanNina,4);
VmeanNina = squeeze(VmeanNina);

VTmeanNina = mean(VTNina,3);
VTmeanNina = mean(VTmeanNina,2);
VTmeanNina =mean(VTmeanNina,4);
VTmeanNina = squeeze(VTmeanNina);

QmeanNoNino = mean(QNoNino,3);
QmeanNoNino = mean(QmeanNoNino,2);
QmeanNoNino =mean(QmeanNoNino,4);
QmeanNoNino = squeeze(QmeanNoNino);

UmeanNoNino = mean(UNoNino,3);
UmeanNoNino = mean(UmeanNoNino,2);
UmeanNoNino =mean(UmeanNoNino,4);
UmeanNoNino = squeeze(UmeanNoNino);

VmeanNoNino = mean(VNoNino,3);
VmeanNoNino = mean(VmeanNoNino,2);
VmeanNoNino =mean(VmeanNoNino,4);
VmeanNoNino = squeeze(VmeanNoNino);

VTmeanNoNino = mean(VTNoNino,3);
VTmeanNoNino = mean(VTmeanNoNino,2);
VTmeanNoNino =mean(VTmeanNoNino,4);
VTmeanNoNino = squeeze(VTmeanNoNino);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%stddev. 

VTmeanLev = mean(VTperiodLev,4);
VTmeanLev = mean(VTmeanLev,3);
VTanom = VTperiodLev - VTmeanLev;
VTint = std(VTanom,0,3);
VTint = mean(VTint,4);
VTint = squeeze(VTint);
VTint = mean(VTint,1);

QmeanLev = mean(QperiodLev,4);
QmeanLev = mean(QmeanLev,3);
Qanom = QperiodLev - QmeanLev;
Qint = std(Qanom,0,3);
Qint = mean(Qint,4);
Qint = squeeze(Qint);
Qint = mean(Qint,1);

VmeanLev = mean(VperiodLev,4);
VmeanLev = mean(VmeanLev,3);
Vanom = VperiodLev - VmeanLev;
Vint = std(Vanom,0,3);
Vint = mean(Vint,4);
Vint = squeeze(Vint);
Vint = mean(Vint,1);


UmeanLev = mean(UperiodLev,4);
UmeanLev = mean(UmeanLev,3);
Uanom = UperiodLev - UmeanLev;
Uint = std(Uanom,0,3);
Uint = mean(Uint,4);
Uint = squeeze(Uint);
Uint = mean(Uint,1);

VTmeanLevNino = mean(VTNino,4);
VTmeanLevNino = mean(VTmeanLevNino,3);
VTanomNino = VTNino - VTmeanLevNino;
VTintNino = std(VTanomNino,0,3);
VTintNino = mean(VTintNino,4);
VTintNino = squeeze(VTintNino);
VTintNino = mean(VTintNino,1);

QmeanLevNino = mean(QNino,4);
QmeanLevNino = mean(QmeanLevNino,3);
QanomNino = QNino - QmeanLevNino;
QintNino = std(QanomNino,0,3);
QintNino = mean(QintNino,4);
QintNino = squeeze(QintNino);
QintNino = mean(QintNino,1);

VmeanLevNino = mean(VNino,4);
VmeanLevNino = mean(VmeanLevNino,3);
VanomNino = VNino - VmeanLevNino;
VintNino = std(VanomNino,0,3);
VintNino = mean(VintNino,4);
VintNino = squeeze(VintNino);
VintNino = mean(VintNino,1);


UmeanLevNino = mean(UNino,4);
UmeanLevNino = mean(UmeanLevNino,3);
UanomNino = UNino - UmeanLevNino;
UintNino = std(UanomNino,0,3);
UintNino = mean(UintNino,4);
UintNino = squeeze(UintNino);
UintNino = mean(UintNino,1);

VTmeanLevNina = mean(VTNina,4);
VTmeanLevNina = mean(VTmeanLevNina,3);
VTanomNina = VTNina - VTmeanLevNina;
VTintNina = std(VTanomNina,0,3);
VTintNina = mean(VTintNina,4);
VTintNina = squeeze(VTintNina);
VTintNina = mean(VTintNina,1);

QmeanLevNina = mean(QNina,4);
QmeanLevNina = mean(QmeanLevNina,3);
QanomNina = QNina - QmeanLevNina;
QintNina = std(QanomNina,0,3);
QintNina = mean(QintNina,4);
QintNina = squeeze(QintNina);
QintNina = mean(QintNina,1);

VmeanLevNina = mean(VNina,4);
VmeanLevNina = mean(VmeanLevNina,3);
VanomNina = VNina - VmeanLevNina;
VintNina = std(VanomNina,0,3);
VintNina = mean(VintNina,4);
VintNina = squeeze(VintNina);
VintNina = mean(VintNina,1);


UmeanLevNina = mean(UNina,4);
UmeanLevNina = mean(UmeanLevNina,3);
UanomNina = UNina - UmeanLevNina;
UintNina = std(UanomNina,0,3);
UintNina = mean(UintNina,4);
UintNina = squeeze(UintNina);
UintNina = mean(UintNina,1);

VTmeanLevNoNino = mean(VTNoNino,4);
VTmeanLevNoNino = mean(VTmeanLevNoNino,3);
VTanomNoNino = VTNoNino - VTmeanLevNoNino;
VTintNoNino = std(VTanomNoNino,0,3);
VTintNoNino = mean(VTintNoNino,4);
VTintNoNino = squeeze(VTintNoNino);
VTintNoNino = mean(VTintNoNino,1);

QmeanLevNoNino = mean(QNoNino,4);
QmeanLevNoNino = mean(QmeanLevNoNino,3);
QanomNoNino = QNoNino - QmeanLevNoNino;
QintNoNino = std(QanomNoNino,0,3);
QintNoNino = mean(QintNoNino,4);
QintNoNino = squeeze(QintNoNino);
QintNoNino = mean(QintNoNino,1);

VmeanLevNoNino = mean(VNoNino,4);
VmeanLevNoNino = mean(VmeanLevNoNino,3);
VanomNoNino = VNoNino - VmeanLevNoNino;
VintNoNino = std(VanomNoNino,0,3);
VintNoNino = mean(VintNoNino,4);
VintNoNino = squeeze(VintNoNino);
VintNoNino = mean(VintNoNino,1);


UmeanLevNoNino = mean(UNoNino,4);
UmeanLevNoNino = mean(UmeanLevNoNino,3);
UanomNoNino = UNoNino - UmeanLevNoNino;
UintNoNino = std(UanomNoNino,0,3);
UintNoNino = mean(UintNoNino,4);
UintNoNino = squeeze(UintNoNino);
UintNoNino = mean(UintNoNino,1);


%% plot. 
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.005], [0.05 0.09], [0.07 0.05]);
if ~make_it_tight,  clear subplot;  end


gg1=figure;
subplot(1,4,2:4);
aa=plot(levs,flipud(VTmean),'LineWidth',4,'Color',[0.4627    0.4039    0.3098]);
jbfill(levs',flipud(VTmean-VTint'),flipud(VTmean+VTint'),[0.4627    0.4039    0.3098],'k',1,0.15);

hold on
bb=plot(levs,flipud(Vmean),'LineWidth',4,'Color',[ 0    0.2667    0.1059]);
jbfill(levs',flipud(Vmean-Vint'),flipud(Vmean+Vint'),[ 0    0.2667    0.1059],'k',1,0.18);

hold on
cc=plot(levs,flipud(Umean),'LineWidth',4,'Color',[0.5961         0    0.2627]);
jbfill(levs',flipud(Umean-Uint'),flipud(Umean+Uint'),[0.5961         0    0.2627],'k',1,0.15);
dood = legend([aa,bb,cc],'Vapor Transport kg/(ms)','V-wind (m/s)','U-wind (m/s)'...
    ,'Location','SouthEast');

yticks([50 500 1000 1500])
dood.FontSize=17;
ylim([0,1500])
set(gca,'xtick',[])
set(gca,'XDir','reverse');
view([90 -90])
title('DJF AR VT, U, & V')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 20)

subplot(1,4,1)
hold on
bb=plot(levs,flipud(Qmean),'LineWidth',4,'Color',[0.0902    0.4039    0.6824]);
xlabel('[hPa]')
jbfill(levs',flipud(Qmean-Qint'),flipud(Qmean+Qint'),[0.0902    0.4039    0.6824],'k',1,0.15);
dood = legend('Q (Kg/Kg)');
dood.FontSize=17;
set(gca,'XDir','reverse');
view([90 -90])
title('DJF AR Q')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 20)

set(gg1,'Position',[183         472        1094         871])


annotation('textbox', [0 0.9 1 0.1], ...
    'String','[Lon: 210E Lat: 22N-42N]     DJF: Q,U,V,VT & Internal Variability', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center','FontSize',25,'FontWeight','bold')
%%
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.09 0.09], [0.05 0.09], [0.05 0.05]);
if ~make_it_tight,  clear subplot;  end

bb = figure;
set(bb,'Position',[183         196        1606        1147])
subplot(2,2,1);
aaa = plot(levs,flipud(VTmean),'LineWidth',4,'Color',[0.0902    0.4039    0.6824]);
xlabel('[hPa]')
jbfill(levs',flipud(VTmean-VTint'),flipud(VTmean+VTint'),[0.0902    0.4039    0.6824],'k',1,0.15);
hold on
bbb = plot(levs,flipud(VTmeanNino),'LineWidth',4,'Color',[0.5961         0    0.2627]);
xlabel('[hPa]')
ylabel('Vapor Transport Kg/(m/s)')
jbfill(levs',flipud(VTmeanNino-VTintNino'),flipud(VTmeanNino+VTintNino'),[0.5961         0    0.2627],'k',1,0.15);

hold on
ccc= plot(levs,flipud(VTmeanNina),'LineWidth',4,'Color',[0.1020    0.1020    0.1020]);
xlabel('[hPa]')
jbfill(levs',flipud(VTmeanNina-VTintNina'),flipud(VTmeanNina+VTintNina'),[0.1020    0.1020    0.1020],'k',1,0.15);

dood = legend([aaa,bbb,ccc],'VT (Kg/Kg)','VTnino (Kg/Kg)','VTnina (Kg/Kg)');
dood.FontSize=17;

set(gca,'XDir','reverse');
view([90 -90])
title('AR VT')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 20)
hold on


subplot(2,2,2);
aaa = plot(levs,flipud(Qmean),'LineWidth',4,'Color',[0.0902    0.4039    0.6824]);
xlabel('[hPa]')
jbfill(levs',flipud(Qmean-Qint'),flipud(Qmean+Qint'),[0.0902    0.4039    0.6824],'k',1,0.15);
hold on
bbb = plot(levs,flipud(QmeanNino),'LineWidth',4,'Color',[0.5961         0    0.2627]);
xlabel('[hPa]')
ylabel('Specific Humidity [Kg/Kg]')
jbfill(levs',flipud(QmeanNino-QintNino'),flipud(QmeanNino+QintNino'),[0.5961         0    0.2627],'k',1,0.15);

hold on
ccc= plot(levs,flipud(QmeanNina),'LineWidth',4,'Color',[0.1020    0.1020    0.1020]);
xlabel('[hPa]')
jbfill(levs',flipud(QmeanNina-QintNina'),flipud(QmeanNina+QintNina'),[0.1020    0.1020    0.1020],'k',1,0.15);

dood = legend([aaa,bbb,ccc],'Q (Kg/Kg)','Qnino (Kg/Kg)','Qnina (Kg/Kg)');
dood.FontSize=17;

set(gca,'XDir','reverse');
view([90 -90])
title('AR Q')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 20)
hold on



subplot(2,2,3);
aaa = plot(levs,flipud(Vmean),'LineWidth',4,'Color',[0.0902    0.4039    0.6824]);
xlabel('[hPa]')
jbfill(levs',flipud(Vmean-Vint'),flipud(Vmean+Vint'),[0.0902    0.4039    0.6824],'k',1,0.15);
hold on
bbb = plot(levs,flipud(VmeanNino),'LineWidth',4,'Color',[0.5961         0    0.2627]);
xlabel('[hPa]')
ylabel('Meridional Wind Speed [m/s]')
jbfill(levs',flipud(VmeanNino-VintNino'),flipud(VmeanNino+VintNino'),[0.5961         0    0.2627],'k',1,0.15);

hold on
ccc= plot(levs,flipud(VmeanNina),'LineWidth',4,'Color',[0.1020    0.1020    0.1020]);
xlabel('[hPa]')
jbfill(levs',flipud(VmeanNina-VintNina'),flipud(VmeanNina+VintNina'),[0.1020    0.1020    0.1020],'k',1,0.15);

dood = legend([aaa,bbb,ccc],'V (Kg/Kg)','Vnino (Kg/Kg)','Vnina (Kg/Kg)');
dood.FontSize=17;

set(gca,'XDir','reverse');
view([90 -90])
title('AR V')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 20)
hold on

subplot(2,2,4);

aaa = plot(levs,flipud(Umean),'LineWidth',4,'Color',[0.0902    0.4039    0.6824]);
xlabel('[hPa]')
jbfill(levs',flipud(Umean-Uint'),flipud(Umean+Uint'),[0.0902    0.4039    0.6824],'k',1,0.15);
hold on
bbb = plot(levs,flipud(UmeanNino),'LineWidth',4,'Color',[0.5961         0    0.2627]);
xlabel('[hPa]')
ylabel('Zonal Wind Speed [m/s]')

jbfill(levs',flipud(UmeanNino-UintNino'),flipud(UmeanNino+UintNino'),[0.5961         0    0.2627],'k',1,0.15);

hold on
ccc= plot(levs,flipud(UmeanNina),'LineWidth',4,'Color',[0.1020    0.1020    0.1020]);
xlabel('[hPa]')
jbfill(levs',flipud(UmeanNina-UintNina'),flipud(UmeanNina+UintNina'),[0.1020    0.1020    0.1020],'k',1,0.15);

dood = legend([aaa,bbb,ccc],'U (Kg/Kg)','Unino (Kg/Kg)','Unina (Kg/Kg)');
dood.FontSize=17;

set(gca,'XDir','reverse');
view([90 -90])
title('AR U')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 20)
hold on


annotation('textbox', [0 0.9 1 0.1], ...
    'String','[Lon: 210E Lat: 22N-42N]     DJF: Q,U,V,VT & Internal Variability + Nino & Nina', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center','FontSize',25,'FontWeight','bold')

annotation('line',[0.49 0.49],[0.97 0.02],'color','k','LineWidth',2);
annotation('line',[0.01 0.978],[0.47 0.47],'color','k','LineWidth',2);
annotation('line',[0.02 0.978],[0.97 0.97],'color','k','LineWidth',2);

