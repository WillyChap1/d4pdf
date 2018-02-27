%prep IVT ensembles 
clear all
cd /Users/wchapman/SIO/research/100Amip/d4pdf/500Ens
s = dir('ivt*.nc');

for ii = 1952:2011 
    disp(ii)
    tic
    for jj = 1:length(s)
        
        filename = s(jj).name;
        time = ncread(filename,'time');
        ivt = ncread(filename,'ivtt');
        lat =ncread(filename,'lat');
        lon =ncread(filename,'lon');

        datestart = '01-Jan-1951 00:00:00';
        formatin = 'dd-MMM-yyyy HH:mm:ss';
        date1 = datetime(datestart,'InputFormat',formatin);

        timevector = date1 + minutes(time);
        dateMat = datevec(timevector);
    
    
        indexget = find(dateMat(:,1)==ii);
    
        IVTtot(:,:,jj,:) = ivt(:,:,indexget);
        timer = timevector(indexget);

    end
    toc
    file2save = strcat('IVTtot______',num2str(ii),'0101_',num2str(ii),'1231.mat');
    save(file2save,'lon','lat','timer','IVTtot')
    
end 