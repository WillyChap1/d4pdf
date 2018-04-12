%D4PDFClustersTS.m
%Kmeans for "find Events" 
close all
clear all
cd /Users/wchapman/SIO/research/100Amip/datasets
numCLUST = 8;
switcher = 'month';


if strcmp(switcher,'event')
    filename = 'AllEnseOneMat_LF_230_240_30_55.mat';
    load(filename);
    lon1d = str2double(filename(18:20)); 
    lon2d = str2double(filename(22:24));
    lat1d = str2double(filename(26:27));
    lat2d = str2double(filename(29:30));
    
else 
    filename = 'ARmonthlyCount.mat';
    load(filename);
    lon = load('AllEnseOneMat_LF_230_240_30_55.mat','lon');
    lon =lon.lon;
    lat = load('AllEnseOneMat_LF_230_240_30_55.mat','lat');
    lat= lat.lat;
    compARcat = zeros(size(ARmonth,1),size(ARmonth,2),size(ARmonth,3)*size(ARmonth,4));
    index1 = 1;
    for kk = 1:size(ARmonth,3)
       compARcat(:,:,index1:index1+707) =  ARmonth(:,:,kk,:);
       comptimeMat(index1:index1+707,:) = [kk*ones(708,1) datenum(timemonth)];
       index1 = index1+708;
    end
    
    lat1d = 20;
    lat2d = 22;
    lon1d =200;
    lon2d = 201;
end
%%%%%%%%%%%%%%%%%%%%% .  select months   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Months1 = [1 2 12];
Months2 = [];


if isempty(Months2)
    MonthSelect = Months1;
else 
    MonthSelect = cat(2,Months1,Months2);   %combined months for plotting.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


crazyDateVec = datevec(comptimeMat(:,2));
c = ismember(crazyDateVec(:,2),MonthSelect);
indMonth = find(c);

compARcat = compARcat(:,:,indMonth);
comptimeMat = comptimeMat(indMonth,:);

%set the ARs for just shape files: 
settoone =1; 
BigCat = reshape(compARcat,[size(compARcat,1)*size(compARcat,2),size(compARcat,3)]);
rotBig = BigCat';

% if settoone == 1
%    BigCat(BigCat>1)= 1;  
% end 


[Tryingit,CC,sumd,D] = kmeans(rotBig,numCLUST,'dist', 'sqEuclidean', 'start', 'cluster');
MonthSelectDatevec = crazyDateVec(indMonth,:);
MonthSelectDatevec(:,4) = Tryingit;
%%
%years loop
disp('Creating TSs')
MonthsTs=[];
timeTS =[];
for ii = 1952:2010
   
   %months loop
   for jj = 1:12 
       yearind = find(MonthSelectDatevec(:,1) == ii);
       subsetMSDV = MonthSelectDatevec(yearind,:);
       monthind = find(subsetMSDV(:,2)==jj);
       subsetMSDV = subsetMSDV(monthind,:);
       
       MonthsTS(jj,1) = length(find(subsetMSDV(:,4)==1));
       MonthsTS(jj,2) = length(find(subsetMSDV(:,4)==2));
       MonthsTS(jj,3) = length(find(subsetMSDV(:,4)==3));
       MonthsTS(jj,4) = length(find(subsetMSDV(:,4)==4));
       MonthsTS(jj,5) = length(find(subsetMSDV(:,4)==5));
       MonthsTS(jj,6) = length(find(subsetMSDV(:,4)==6));
       MonthsTS(jj,7) = length(find(subsetMSDV(:,4)==7));
       MonthsTS(jj,8) = length(find(subsetMSDV(:,4)==8));
       timevec(jj,:) = [ii,jj];
       
   end 
    
    MonthsTs = cat(1,MonthsTs,MonthsTS);
    timeTS = cat(1,timeTS,timevec);
end
%MonthsTs = MonthsTs(any(MonthsTs,2),:);
save('/Users/wchapman/SIO/research/100Amip/datasets/6clustTS.mat','MonthsTs','timeTS');
plot(MonthsTs)

%% plotting

oneK = find(Tryingit == 1);
twoK = find(Tryingit == 2);
threeK = find(Tryingit == 3);
fourK = find(Tryingit == 4);
fiveK = find(Tryingit == 5);
sixK = find(Tryingit == 6);
sevK = find(Tryingit == 7);
eigK = find(Tryingit == 8);

oneKmat = compARcat(:,:,oneK);
twoKmat = compARcat(:,:,twoK);
threeKmat = compARcat(:,:,threeK);
fourKmat = compARcat(:,:,fourK);
fiveKmat = compARcat(:,:,fiveK);
sixKmat = compARcat(:,:,sixK);
sevKmat = compARcat(:,:,sevK);
eigKmat = compARcat(:,:,eigK);

oneKsum = sum(oneKmat,3)./length(oneK);
twoKsum = sum(twoKmat,3)./length(twoK);
threeKsum = sum(threeKmat,3)./length(threeK);
fourKsum = sum(fourKmat,3)./length(fourK);
fiveKsum = sum(fiveKmat,3)./length(fiveK);
sixKsum = sum(sixKmat,3)./length(sixK);
sevKsum = sum(sevKmat,3)./length(sevK);
eigKsum = sum(eigKmat,3)./length(eigK);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%      %Plotting%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create MonthTitle.
titMon = [];
mons =['J','F','M','A','M','J','J','A','S','O','N','D'];
for yy = 1:length(MonthSelect)
    
    mon = MonthSelect(yy);
    
    monC = mons(mon);
    
    titMon = strcat(titMon,monC);
       
end 

                                %options: [gap]   %hght [Top Bot]  %Wid  [L R]                          
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
if ~make_it_tight,  clear subplot;  end



heyhey = figure;
hold on
ax11=subplot(2,2,1);
padvar =4;
padCOMPone = padarray(oneKsum,[padvar,padvar],nan,'both');
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

ARplot = padCOMPone;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
boxlon=[lon1d lon1d lon2d lon2d lon1d];
boxlat=[lat1d lat2d lat2d lat1d lat1d];
m_line(boxlon,boxlat,'linewi',2,'color','r')
m_contourf(lon,lat,ARplot',15);
RuBU = cbrewer('seq','Greys',15);
colormap(ax11,RuBU);
tit = strcat(titMon,{' '},'D4pdf Kmeans 1: Mean AR composite',{' '},'ARs:',{' '},num2str(length(oneK)));
title(tit)
%caxis([0 100])
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
ylabel(h1,'Number of Event Periods');
%-=========================================================================

ax12=subplot(2,2,2);
padvar =4;
padCOMPtwo = padarray(twoKsum,[padvar,padvar],nan,'both');

ARplot = padCOMPtwo;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,ARplot',15);
RuBU = cbrewer('seq','Greys',15);
colormap(ax12,RuBU);
tit = strcat(titMon,{' '},'D4pdf Kmeans 2: Mean AR composite',{' '},'ARs:',{' '},num2str(length(twoK)));
title(tit)%caxis([0 100])
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
ylabel(h1,'Number of Event Periods');

%-=========================================================================

ax13 = subplot(2,2,3);
padvar =4;
padCOMPthree = padarray(threeKsum,[padvar,padvar],nan,'both');

ARplot = padCOMPthree;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,ARplot',15);
RuBU = cbrewer('seq','Greys',15);
colormap(ax13,RuBU);
tit = strcat(titMon,{' '},'D4pdf Kmeans 3: Mean AR composite',{' '},'ARs:',{' '},num2str(length(threeK)));
title(tit)%caxis([0 100])
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
ylabel(h1,'Number of Event Periods');


%=========================================================================

ax14 = subplot(2,2,4);
padvar =4;
padCOMPfour = padarray(fourKsum,[padvar,padvar],nan,'both');

ARplot = padCOMPfour;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,ARplot',15);
RuBU = cbrewer('seq','Greys',15);
colormap(ax14,RuBU);
tit = strcat(titMon,{' '},'D4pdf Kmeans 4: Mean AR composite',{' '},'ARs:',{' '},num2str(length(fourK)));
title(tit)
%caxis([0 100])
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
ylabel(h1,'Number of Event Periods');
pos = get(heyhey,'position');
set(heyhey,'position',[193          81        1831        1264])
%=========================================================================

yogabba = figure; 
ax14 = subplot(2,2,1);
padvar =4;
padCOMPfour = padarray(fiveKsum,[padvar,padvar],nan,'both');

ARplot = padCOMPfour;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,ARplot',15);
RuBU = cbrewer('seq','Greys',15);
colormap(ax14,RuBU);
tit = strcat(titMon,{' '},'D4pdf Kmeans 5: Mean AR composite',{' '},'ARs:',{' '},num2str(length(fiveK)));
title(tit)
%caxis([0 100])
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
ylabel(h1,'Number of Event Periods');

%=========================================================================


ax14 = subplot(2,2,2);
padvar =4;
padCOMPfour = padarray(sixKsum,[padvar,padvar],nan,'both');

ARplot = padCOMPfour;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,ARplot',15);
RuBU = cbrewer('seq','Greys',15);
colormap(ax14,RuBU);
tit = strcat(titMon,{' '},'D4pdf Kmeans 6: Mean AR composite',{' '},'ARs:',{' '},num2str(length(sixK)));
title(tit)
%caxis([0 100])
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
ylabel(h1,'Number of Event Periods');


%=========================================================================


ax14 = subplot(2,2,3);
padvar =4;
padCOMPfour = padarray(sevKsum,[padvar,padvar],nan,'both');

ARplot = padCOMPfour;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,ARplot',15);
RuBU = cbrewer('seq','Greys',15);
colormap(ax14,RuBU);
tit = strcat(titMon,{' '},'D4pdf Kmeans 7: Mean AR composite',{' '},'ARs:',{' '},num2str(length(sevK)));
title(tit)
%caxis([0 100])
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
ylabel(h1,'Number of Event Periods');


%=========================================================================


ax14 = subplot(2,2,4);
padvar =4;
padCOMPfour = padarray(eigKsum,[padvar,padvar],nan,'both');

ARplot = padCOMPfour;
m_proj('miller','lon',[double(lon(1)), double(lon(end))],'lat',[double(lat(1)), double(lat(end))]);
hold on
m_grid('box','fancy','tickdir','in');
m_coast('linewidth',2,'Color','k');
m_contourf(lon,lat,ARplot',15);
RuBU = cbrewer('seq','Greys',15);
colormap(ax14,RuBU);
tit = strcat(titMon,{' '},'D4pdf Kmeans 8: Mean AR composite',{' '},'ARs:',{' '},num2str(length(eigK)));
title(tit)
%caxis([0 100])
ax = gca;
ax.FontSize = 25;
h1  = colorbar;
ylabel(h1,'Number of Event Periods');
pos = get(yogabba,'position');
set(yogabba,'position',[193          81        1831        1264])


