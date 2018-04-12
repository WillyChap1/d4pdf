%find events. 
clear all

cd /Volumes/HotSlop/d4pdf/ivtEns/IVTmasked
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   desired Lat/Lon    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%test between this lat and lon. 

desLon =[230 240];
desLat =[30 55];

%%%%%%%%%%%%%%%%%%%%%%%%%%   desired Lat/Lon    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = dir('IVTmask*.mat');

ARmonth = [];
timemonth =[];
monther =1:12;
timeinc = 0:6:8784;

for hh = 1:length(s)
    filename = s(hh).name;
    disp(filename)
    %set my timer increment.  
    year = filename(14:17);
    datestart = strcat(year,'-01-01 00:00');
    formatin = 'yyyy-MM-dd HH:mm';
    date1 = datetime(datestart,'InputFormat',formatin);
    timetot = date1+hours(timeinc);
    [yearsdo,~]=datevec(timetot);
    timetot = timetot(yearsdo==str2double(year));
    
    [yearsdo,monthsdo,~]=datevec(timetot);
    
    %load variables.
    load(filename);
    ARtot = IVTmask(:,:,:,1:length(timetot)); %shorten to just this year. 
    
    ARmonthTemp=[];
    timemonthTemp = [];
    for jj = 1:length(monther)
        ARtemp = ARtot(:,:,:,monthsdo ==jj);
        ARtemp = nansum(ARtemp,4);
        ARmonthTemp = cat(4,ARmonthTemp,ARtemp);
        timetemp = find(monthsdo==jj,1);
        timetemp = timetot(timetemp);
        timemonthTemp =cat(1,timemonthTemp,timetemp);
    end
    
    ARmonth=cat(4,ARmonth,ARmonthTemp);
    timemonth = cat(1,timemonth,timemonthTemp);

end
latmonth = lat;
lonmonth = lon; 
save('/Users/wchapman/SIO/research/100Amip/datasets/AR_IVT_monthlyCount.mat',...
    'ARmonth','timemonth','latmonth','lonmonth');

