clear all 
close all

cd /Users/wchapman/SIO/research/100Amip/datasets
filename = 'ARmonthlyCount.mat';
load(filename);
lon = load('AllEnseOneMat_LF_230_240_30_55.mat','lon');
lon =lon.lon;
lat = load('AllEnseOneMat_LF_230_240_30_55.mat','lat');
lat= lat.lat;

compARcat = squeeze(mean(ARmonth,3));
%get desired months
MonthSelect = [12,1,2];
crazyDateVec = datevec(timemonth);
c = ismember(crazyDateVec(:,2),MonthSelect);
indMonth = find(c);

compARcat = compARcat(:,:,indMonth);
timenew = timemonth(indMonth);

%perform EOF
[eof_map,pc,expvar] = eof3d(compARcat,timenew,10);

%get clusters 
load('/Users/wchapman/SIO/research/100Amip/datasets/6clustTS.mat')

ARts = MonthsTs(:,5);

monther = timeTS(:,2);
for k = 1:12
        indy = monther ==k;
        monthlymeans(k) = mean(ARts(indy));
        ARts(indy) = bsxfun(@minus,ARts(indy),monthlymeans(k));
end 

c = ismember(timeTS(:,2),MonthSelect);
indMonth = find(c);
timeTRY = timeTS(indMonth,:);
ARts = ARts(indMonth);

corr(pc(:,1),ARts)



