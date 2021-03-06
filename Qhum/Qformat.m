%prep IVT ensembles 
clear all

lev = 850;

dirsearch = strcat('Q',num2str(lev),'*.nc');
cd /Users/wchapman/SIO/research/100Amip/d4pdf/Qlev/
s = dir(dirsearch);

filename = s(1).name;
lat =ncread(filename,'lat');
lon =ncread(filename,'lon');
nlat =length(lat);
nlon =length(lon);

for ii = 1952:2011 
    disp(ii)
    tic
    Qtot = zeros(nlon,nlat,length(s),1464);
    for jj = 1:length(s)
        disp(jj)
        tic
        filename = s(jj).name;
       
        time = ncread(filename,'time');
        q = ncread(filename,'qa');
        lat =ncread(filename,'lat');
        lon =ncread(filename,'lon');
        toc
        
        datestart = '01-Jan-1951 00:00:00';
        formatin = 'dd-MMM-yyyy HH:mm:ss';
        date1 = datetime(datestart,'InputFormat',formatin);
        
        timevector = date1 + minutes(time);
        dateMat = datevec(timevector);
       
        indexget = find(dateMat(:,1)==ii);
       
        Qtot(:,:,jj,1:length(indexget)) = q(:,:,indexget);
       
        timer = timevector(indexget);
        
    end
    toc
    Qtot = Qtot(:,:,:,length(timer));
    file2save = strcat(filename(1:4),'______',num2str(ii),'0101_',num2str(ii),'1231.mat');
    save(file2save,'lon','lat','timer','Qtot')
    
end 