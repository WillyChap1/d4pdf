%make mean SST global --- remove the seasonal cycle. 
clear all
close all
count =1;
firstlag=0;
jr=0;
for mm =firstlag:firstlag+2

cd /Volumes/HotSlop/d4pdf/SST/
%load my variables --->
load('SSTglobeAve.mat')

%change this for frequency vs IVT ==>
load('/Users/wchapman/SIO/research/100Amip/datasets/AR_IVT_monthlyCount.mat')
%<==

sstAVE = sstAVE(:,:,1:length(timemonth));
timeTot = timeTot(1:length(timemonth));
%----------------------- Switches -----------------------------------------
switcher = 'ARcount';   %either ARcount or cluster regresions
clustnum = jr;
pcNum = 1;
latROI = [20 35];
lonROI = [200 235];
%--------------------------------------------------------------------------
%months to test on ARs 
trymonths =[1 2 12];

%sea surface lag to try 
lag = mm;  %months. 

%% cluster regressions 
if strcmp(switcher,'cluster')
    % TS created using script in Users/wchapman/Git/d4pdf --> D4PDFClustersTS.m
    load('/Users/wchapman/SIO/research/100Amip/datasets/6clustTS.mat')
    ARts = MonthsTs(:,clustnum);
    
    %MOD!!!!!!
    %ARts = MonthsTs(:,4)+MonthsTs(:,5)+MonthsTs(:,8);
    
    
    % remove seasonal cycle and mean from SST 
    timemon = (1:length(timeTot))';
    for ii = 1:length(lon)
        for jj = 1:length(lat)
            TSr = squeeze(sstAVE(ii,jj,:)); 
            An = [ones(length(timeTot),1) timemon(:) sin(2*pi*timemon/12) cos(2*pi*timemon/12)];
            x2 = (An'*An)\An'*TSr;
            Tfit = An*x2; 
            sstAVE(ii,jj,:) = squeeze(sstAVE(ii,jj,:))-Tfit;
        end 
    end
    
    %remove seasonal cycle... by removing mean state per month (anomaly).
    
    
    [~,monther,~]=datevec(timeTot);

    % Demean by month 

    for k = 1:12
        indy = monther ==k;
        monthlymeans(k) = mean(ARts(indy));
        ARts(indy) = bsxfun(@minus,ARts(indy),monthlymeans(k));
    end 
    
    
end
%% AR count regressions 
if strcmp('ARcount',switcher)
    
    if mm == firstlag
    cd /Users/wchapman/SIO/research/100Amip/d4pdf/
    Im = imread('RegionSelect2.png');
    BW = roipoly(Im);
    BB = imresize(BW,[34,45]);

    load('/Users/wchapman/SIO/research/100Amip/d4pdf/allyear/bigsquare/EXPtot_1_13_19520101_19521231.mat')
    
    close all
    %repositions BBs. 
    BBs = double(BB);
    BBs=BBs';
    BBs = fliplr(BBs);
    BBs(BBs==0)=nan;
    cd /Volumes/HotSlop/d4pdf/SST/ 
    end 
    
    
    ARts =BBs.*ARmonth;
    ARts = squeeze(nanmean(nanmean(nanmean(ARts,1),2),3));
    load SSTglobeAve.mat lat lon
    
% remove seasonal cycle and mean. 
    timemon = (1:length(timeTot))';
    for ii = 1:length(lon)
         for jj = 1:length(lat)
            TSr = squeeze(sstAVE(ii,jj,:)); 
            An = [ones(length(timeTot),1) timemon(:) sin(2*pi*timemon/12) cos(2*pi*timemon/12)];
            x2 = (An'*An)\An'*TSr;
            Tfit = An*x2; 
            sstAVE(ii,jj,:) = squeeze(sstAVE(ii,jj,:))-Tfit;
         end 
    end

    An = [ones(length(timeTot),1) timemon(:) sin(2*pi*timemon/12) cos(2*pi*timemon/12) sin(4*pi*timemon/12) cos(4*pi*timemon/12) sin(6*pi*timemon/12) cos(6*pi*timemon/12) sin(8*pi*timemon/12) cos(8*pi*timemon/12)];
    x2 =(An'*An)\An'*ARts;
    ARfit = An*x2;
    ARts=ARts-ARfit;
end

%%
if strcmp('EOF',switcher)
    
    [~,latinS] = (min(abs(latmonth-latROI(1))));
    [~,latinE] = (min(abs(latmonth-latROI(2))));
    [~,loninS] = (min(abs(lonmonth-lonROI(1))));
    [~,loninE] = (min(abs(lonmonth-lonROI(2))));
    
    ARts = ARmonth(loninS:loninE,latinS:latinE,:,:);
    ARts = squeeze(mean(ARts,3));
    
    [eof_map,pcs,expVAR]=eof3d(ARts,timeTot,10);
   
    ARts = pcs(:,pcNum);
    
    timemon = (1:length(timeTot))';
    for ii = 1:length(lon)
         for jj = 1:length(lat)
            TSr = squeeze(sstAVE(ii,jj,:)); 
            An = [ones(length(timeTot),1) timemon(:) sin(2*pi*timemon/12) cos(2*pi*timemon/12)];
            x2 = (An'*An)\An'*TSr;
            Tfit = An*x2; 
            sstAVE(ii,jj,:) = squeeze(sstAVE(ii,jj,:))-Tfit;
         end 
    end
    
end
%% exclude months and include the lag 
[~,monthin] = datevec(timemonth);
Lia = ismember(trymonths,monthin);
idx=find(Lia);
map=bsxfun(@eq, trymonths(idx),monthin(:));  %this is tricky. look it up. help bsxfun
map = logical(sum(map,2));

timedoAR = timemonth(map);
ARts=ARts(map);

if lag>0
    ARts =ARts((1+lag):end);
    map = circshift(map,-lag);
    sstAVE = sstAVE(:,:,map);
    sstAVE = sstAVE(:,:,1:(end-lag));
    timedoAR = timedoAR((1+lag):end);
    timedoSST = timemonth(map);
    timedoSST = timedoSST(1:end-lag);

else
    sstAVE =sstAVE(:,:,map);
    timedoSST = timemonth(map);
end

%%
for ii = 1:length(lon)
    
    if ~(mod(ii,25))
        disp(strcat('Im at',{' '},num2str(ii)))
    end
    
    for jj = 1:length(lat)
        
        
        TS11 = squeeze(sstAVE(ii,jj,:));       
        TS11 = detrend(TS11);
        A = [ARts];
        r2=corr(A,TS11);
        [~,pvv] = corr(A,TS11);
        if pvv<=0.05
            r2Mat(ii,jj)=r2;
        else 
            r2Mat(ii,jj)=0;
        end
        
            
    end
    
end

%%

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.04], [0.02 0.09], [0.03 0.03]);
if ~make_it_tight,  clear subplot;  end

if count==1
    ax13 = figure;
end
subplot(2,2,count)
plotSST = r2Mat;
r2Mat(r2Mat==0)=NaN;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
% if strcmp('ARcount',switcher) || strcmp('EOF',switcher)
% m_line([lonROI(1) lonROI(1) lonROI(2) lonROI(2) lonROI(1)],[latROI(1)...
%     latROI(2) latROI(2) latROI(1) latROI(1)],'LineWidth',3,'Color','b')
% axr = legend('R.O.I.');
% axr.FontSize = 15;
% end
hold on
m_grid('box','fancy','tickdir','in');
m_coast('patch',[.9 .9 .9],'edgecolor','k');
m_contourf(lon,lat,plotSST',18);
title(strcat('lag:',num2str(lag),{' '},'months - Cluster:',{' '},num2str(clustnum)),'FontSize',18) 
RuBU = cbrewer('div','RdBu',15);
colormap(ax13,RuBU);
caxis([-0.7,0.7]);
colorbar

set(ax13,'Position',[427         195        1636        1149])
count = count+1;
end


aaa = subplot(2,2,4);
BBs = double(BB);
m_proj('miller','lon',[double(lonmonth(1)),double(lonmonth(end))],'lat',[double(latmonth(1)),double(latmonth(end))]);
hold on
m_contour(double(lonmonth),double(latmonth),flipud((BBs)),1,'LineWidth',4);
shading flat
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',1,'Color','k');
RuBU = cbrewer('seq','Blues',15);
colormap(aaa,RuBU);
title('Selected AR Region');

annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Global SST Correlation DJF AR Averge vs SST', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center','FontSize',30)



%'R^2 of SST regressed on AR Count in R.O.I.'