clear all
%close all
count =1;
firstlag=0;

cd /Volumes/HotSlop/d4pdf/SST/
%load my variables --->
load('SSTglobeAve.mat')
load('/Users/wchapman/SIO/research/100Amip/datasets/ARmonthlyCount.mat')

sstAVE = sstAVE(:,:,1:length(timemonth));
timeTot = timeTot(1:length(timemonth));
%----------------------- Switches -----------------------------------------
switcher = 'ARcount';   %either ARcount or cluster regresions
clustnum = 1;
pcNum = 1;
mm =0;
latROI = [20 35];
lonROI = [200 235];
%--------------------------------------------------------------------------
%months to test on ARs 
trymonths =[1 2 3 4 5 6 7 8 9 11 12];

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
    figure
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
            An = [ones(length(timeTot),1) timemon(:) sin(2*pi*timemon/12) cos(2*pi*timemon/12) sin(4*pi*timemon/12) cos(4*pi*timemon/12)];
            x2 = (An'*An)\An'*TSr;
            Tfit = An*x2; 
            sstAVE(ii,jj,:) = squeeze(sstAVE(ii,jj,:))-Tfit;
         end 
    end
    plot(timeTot,ARts)
    hold on
    An = [ones(length(timeTot),1) timemon(:) sin(2*pi*timemon/12) cos(2*pi*timemon/12) sin(4*pi*timemon/12) cos(4*pi*timemon/12) sin(6*pi*timemon/12) cos(6*pi*timemon/12) sin(8*pi*timemon/12) cos(8*pi*timemon/12)];
    x2 =(An'*An)\An'*ARts;
    
    ARfit = An*x2;
    ARts=ARts-ARfit;
    plot(timeTot,ARts)
    hold on
    plot(timeTot,ARfit)
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


% % 
% lower_freq = 0.0; %what should we keep- low . keep everything 1 per 5 years
% upper_freq = 0.6; %what should we keep- high .              & 2 per year 

dt = 1/12;     %Time between sample (units of months);
Fs = 1/dt;     %sampling period; 
N = length(ARts); %length plot of signal
dF = (Fs/N); 
f = (-Fs/2:dF:Fs/2-dF)';
%create my band pass filter: 
% BPF = ((lower_freq < abs(f)) & (abs(f) < upper_freq));

t = dt*(0:N-1); 
spek = ifftshift(fft(ARts));

figure 
semilogx((f),abs(spek));


figure
speka = ifftshift(fft(ARfit));
semilogx((f),abs(speka));
title('FFT of AR time series','FontSize',25)
ylabel('Intensity','FontSize',20)
xlabel('frequency [yr^{-1}]','FontSize',20)

% %filter it. 
% BPF = double(BPF);
% filtspek = zeros(length(spek),1);
% filtspek(BPF==1) = spek(BPF==1);

% subplot(2,1,2);
% plot(f,abs(filtspek));
% MonsTsFilt=ifft(ifftshift(filtspek));

% figure
% plot(MonsTsFilt)
% hold on
% plot(ARts);
% title('High Pass vs signal')
% legend('High Pass', 'Signal')