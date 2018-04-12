%find Lagged ROI correlations. 
%first. Create Time Series. 
clear all
close all
cd /Users/wchapman/SIO/research/100Amip/d4pdf/allyear/AllEns/
load('DJF_ARcatMean.mat')                   %lat lon total_DJF_ARs Time(0).
%load('NDJFM_ARcatMean.mat')
%==========================================================================
%define Region of Interest. 
%==========================================================================
monthsall = [12,1,2];     %set DJF. 

desLon =[220 235];
desLat =[30 48];

[~,latind1] = min(abs(lat - desLat(1)));
[~,latind2] = min(abs(lat - desLat(2)));
[~,lonind1] = min(abs(lon - desLon(1)));
[~,lonind2] = min(abs(lon - desLon(2)));
%==========================================================================

%create DJF time series. 

MonsTs = ARcatMemMean(lonind1:lonind2,latind1:latind2,:,:);
MonsTs = mean(MonsTs,3);
MonsTs = reshape(MonsTs,size(MonsTs,1)*size(MonsTs,2),size(MonsTs,4));
MonsTs = mean(MonsTs,1);


titMon = [];
mons =['J','F','M','A','M','J','J','A','S','O','N','D'];
for yy = 1:length(monthsall)
    
    mon = monthsall(yy);
    
    monC = mons(mon);
    
    titMon = strcat(titMon,monC);      
end 
figure
MonsTs = detrend(MonsTs);
plot(timeCat,MonsTs)
datetick('x','yy')
tit =strcat(titMon,{' '},'TS of AR landFall events in ROI');
title(tit)
xlabel('Time')
ylabel('Number of Events')
hold on 

% load('/Users/wchapman/SIO/research/100Amip/datasets/DJF_pc34noNino.mat')
% MonsTs = pc34_2';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FFT hi-pass filter %%%%%%%%%%%%%%%%%%%%%%%%%
%5 year filter 
lower_freq = 0.0; %what should we keep- low . keep everything 1 per 5 years
upper_freq = 0.6; %what should we keep- high .              & 2 per year 

dt = 1;     %Time between sample (1per year);
Fs = 1/dt;  %sampling period; 
N = length(MonsTs); %length of signal
dF = Fs/N; 
f = (-Fs/2:dF:Fs/2-dF)';
%create my band pass filter: 
BPF = ((lower_freq < abs(f)) & (abs(f) < upper_freq));

t = dt*(0:N-1); 
MonsTs = detrend(MonsTs);
spek = fftshift(fft(MonsTs));

figure 
subplot(2,1,1)
plot(f,abs(spek));

%filter it. 
BPF = double(BPF);
filtspek = zeros(length(spek),1);
filtspek(BPF==1) = spek(BPF==1);

subplot(2,1,2);
plot(f,abs(filtspek));

MonsTsFilt=ifft(ifftshift(filtspek));

figure
plot(MonsTsFilt)
hold on
plot(MonsTs);
title('High Pass vs signal')
legend('High Pass', 'Signal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FFT hi-pass filter %%%%%%%%%%%%%%%%%%%%%%%%%

MonsTs = MonsTsFilt'; 

% load('/Users/wchapman/SIO/research/100Amip/datasets/DJF_WholeROI_pc.mat')
% MonsTs = pc_3';

%==========================================================================
%find lagged correlations. USING MonsTs. 
%==========================================================================
%%
Ninomonthsall = [12,1,2];

nino34 = nino_3_get(1952,2009,Ninomonthsall,1);  %function to get our nino index shit..
nino34 = detrend(ninso34);

[corraa,pval] = corr(nino34,MonsTs');

grabperiod = [8 9 10 11 12 1 2 3 4 5 6 7 8 9 10 11 12 1 2 3 4 5 6 7 8 9 10 11 12 1 2 3 4 5]; %DJF
yearstart =1953;
yearend = 2010;
count =1;
Xmonstr =[];
for ii = length(grabperiod):-1:3
        
    Ninomonthsall = grabperiod(ii-2:ii); 
    
    
    if Ninomonthsall(1) == 12
        yearstart = yearstart-1;
        yearend = yearend -1;
    end
    
    %get index
    nino34 = nino_3_get(yearstart,yearend,Ninomonthsall,1);  %function to get our nino index shit..
    PMM = PMM_get(yearstart,yearend,Ninomonthsall,1);
    
    %detrend demean
    nino34 = detrend(nino34);
    PMM = detrend(PMM);
    
    
    %hi-pass filter them 
    lower_freq = 0.0; %what should we keep- low . keep everything 1 per 5 years
    upper_freq = 0.6; %what should we keep- high .              & 2 per year 

    dt = 1;     %Time between sample (1per year);
    Fs = 1/dt;  %sampling period; 
    N = length(nino34); %length of signal
    dF = Fs/N; 
    f = (-Fs/2:dF:Fs/2-dF)';
    %create my band pass filter: 
    BPF = ((lower_freq < abs(f)) & (abs(f) < upper_freq));
    
    %nino 
    spekN = fftshift(fft(nino34'));
    BPF = double(BPF);
    filtspekN = zeros(length(spekN),1);
    filtspekN(BPF==1) = spekN(BPF==1);
    nino34=ifft(ifftshift(filtspekN));
    
    %PMM
    spekP = fftshift(fft(PMM'));
    BPF = double(BPF);
    filtspekP = zeros(length(spekP),1);
    filtspekP(BPF==1) = spekP(BPF==1);
    PMM=ifft(ifftshift(filtspekP));

    
    %Correlate 
    corraN(ii)= corr(nino34,MonsTs');
    [~,pvalsN(ii)] = corr(nino34,MonsTs');
    corraP(ii)= corr(PMM,MonsTs');
    [~,pvalsP(ii)] = corr(PMM,MonsTs');
    
   count = count+1;
   disp(ii)
   Xmonstr = vertcat(Xmonstr, mons(Ninomonthsall)); 
    
end 


%just for the NIO.
yearstart =1955;
yearend = 2010;
count =1;
Xmonstr2 =[];
for ii = length(grabperiod):-1:3
        
    Ninomonthsall = grabperiod(ii-2:ii);   
    if Ninomonthsall(1) == 12
        yearstart = yearstart-1;
        yearend = yearend -1;
    end
    
    %get index
    NIO = NIO_d4pdf_get(yearstart,yearend,Ninomonthsall,1);
    
    TWNP = TWNP_d4pdf_get(yearstart,yearend,Ninomonthsall,1);
    %detrend demean
    
    SST = SST_d4pdf_get(yearstart,yearend,Ninomonthsall,1);
    
    NIO = detrend(NIO);
    TWNP = detrend(TWNP);
    SST = detrend(SST);
    
    
    %hi-pass filter them 
    lower_freq = 0.0; %what should we keep- low . keep everything 1 per 10 years
    upper_freq = 0.6; %what should we keep- high .              & 2 per year 

    dt = 1;     %Time between sample (1per year);
    Fs = 1/dt;  %sampling period; 
    N = length(NIO); %length of signal
    dF = Fs/N; 
    f = (-Fs/2:dF:Fs/2-dF)';
    %create my band pass filter: 
    BPF = ((lower_freq < abs(f)) & (abs(f) < upper_freq));
    
    %NIO 
    spekNI = fftshift(fft(NIO'));
    BPF = double(BPF);
    filtspekNI = zeros(length(spekNI),1);
    filtspekNI(BPF==1) = spekNI(BPF==1);
    NIO=ifft(ifftshift(filtspekNI));
    
    %TWNP
    spekT = fftshift(fft(TWNP'));
    BPF = double(BPF);
    filtspekT = zeros(length(spekT),1);
    filtspekT(BPF==1) = spekT(BPF==1);
    TWNP=ifft(ifftshift(filtspekT));
    
    
    %Correlate 
    corraNI(ii)= corr(NIO,MonsTs(3:end)');
    [~,pvalsNI(ii)] = corr(NIO,MonsTs(3:end)');    
    
    corraT(ii)= corr(TWNP,MonsTs(3:end)');
    [~,pvalsT(ii)] = corr(TWNP,MonsTs(3:end)');
    
    
    corraS(ii)= corr(SST,MonsTs(3:end)');
    [~,pvalsS(ii)] = corr(SST,MonsTs(3:end)');
    
    
    
   count = count+1;
   disp(ii)
   Xmonstr2 = vertcat(Xmonstr2, mons(Ninomonthsall)); 
    
end 


corraN(corraN==0)=[];
corraP(corraP==0)=[];
corraNI(corraNI==0)=[];
corraT(corraT==0)=[];
corraS(corraS==0)=[];


pvalsN(pvalsN==0)=[];
pvalsP(pvalsP==0)=[];
pvalsNI(pvalsNI==0)=[];
pvalsT(pvalsT==0)=[];
pvalsS(pvalsS==0)=[];



Xmonstr = flipud(Xmonstr);

ax1=figure;
xdummy = 1:1:32;

xdummynanN = xdummy;
xdummynanN((pvalsN>0.05))=nan;

xdummynanP = xdummy;
xdummynanP((pvalsP>0.05))=nan;

xdummynanNI = xdummy;
xdummynanNI((pvalsNI>0.05))=nan;

xdummynanT = xdummy;
xdummynanT((pvalsT>0.05))=nan;

xdummynanS = xdummy;
xdummynanS((pvalsS>0.05))=nan;

col1 = [rand rand rand];
col2 = [rand rand rand];
col3 = [rand rand rand];
col4 = [rand rand rand];
col5 = [rand rand rand];

%plot(xdummy,corraN,'Color',[0.4000 0.7608 0.6471])
plot(xdummy,corraN,'Color',col1)

hold on
%plot(xdummynanN,corraN,'o','Color',[0.4000 0.7608 0.6471])
%h1=plot(xdummynanN,corraN,'Color',[0.4000 0.7608 0.6471],'LineWidth',4);
h1=plot(xdummynanN,corraN,'Color',col1,'LineWidth',4);

hold on
%plot(xdummy,corraP,'Color',[0.9882 0.5529 0.3843])
plot(xdummy,corraP,'Color',col2)

%h2=plot(xdummynanP,corraP,'Color',[0.9882 0.5529 0.3843],'LineWidth',4);
h2=plot(xdummynanP,corraP,'Color',col2,'LineWidth',4);

hold on
plot(xdummy,corraNI,'Color',col3)
h3=plot(xdummynanNI,corraNI,'Color',col3,'LineWidth',4);

% hold on
% plot(xdummy,corraT,'Color',col4)
% h4=plot(xdummynanT,corraT,'Color',col4,'LineWidth',4);

hold on
plot(xdummy,corraT,'Color',col5)
h5=plot(xdummynanT,corraT,'Color',col5,'LineWidth',4);


indexestotick = [2 5 8 11 14 17 20 23 26 29];
hold on 
plot([length(Xmonstr)-3,length(Xmonstr)-3],[-1,1],'r--')
plot([0,length(Xmonstr)+7],[0,0],'k')



doggy= Xmonstr(indexestotick,:);


for uu = 1:length(doggy)
    doggycell{uu,1} = doggy(uu,:);
end 

doggycell{10,1} = 'DJF_{0}';
doggycell{6,1} = 'DJF_{1}';
doggycell{2,1} = 'DJF_{-2}';


xlim([0,34]);
xticks(indexestotick);
ylabel('Correlation','FontSize',25)
xlabel('Months','FontSize',25)
tit =strcat('Lagged Corr AR count',{' '},'(',titMon,')',{' '},'vs. Index');
title(tit,'FontSize',30)
xticklabels(doggycell)
lgd=legend([h1,h2,h3,h5],'Nino34','PMM','NIO','SST Global','Location','southwest');
lgd.FontSize = 25;
xAX = get(gca,'XAxis');
xAX.FontSize = 22;
yAX = get(gca,'YAxis');
yAX.FontSize = 22;
set(ax1,'Position',[1000         615         909         723])




figure
lon1d = desLon(1);
lon2d = desLon(2);
lat1d = desLat(1);
lat2d = desLat(2);
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
boxlon=[lon1d lon1d lon2d lon2d lon1d];
boxlat=[lat1d lat2d lat2d lat1d lat1d];
m_line(boxlon,boxlat,'linewi',2,'color','r')









