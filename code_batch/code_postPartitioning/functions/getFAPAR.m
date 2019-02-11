function siteFAPAR = getFAPAR(site, fluxTime)
% this function will...
% take site name, years and time series
% read in the fAPAR data
% and interpolate it to the native resolution

% returns 


% if ~isvector(fluxTime)
%     error('Input must be a vector')
% end

% read the fAPAR data for this site
filename=strcat('../../modis-land-product-subset-master/Rcode/formattedData/',...
            site,'_MOD15A2_Fpar_1km.csv');
datainFAPAR = csvread(filename,9,2);

% get the median fAPAR across pixels
% set zero and 255 to NaN
datainFAPARdata=datainFAPAR(:,3:end);
indX=(datainFAPARdata>100) | (datainFAPARdata <=5);
datainFAPARdata(indX)=NaN;
tmp=smooth(datainFAPARdata,3);
c = smooth(datainFAPARdata(:),3);
C1 = reshape(c,length(datainFAPARdata),9);
C1(indX)=NaN;
datainFAPAR(:,3:end)=C1;

% column 6 is the center pixel
% work with this for now
siteFAPAR=datainFAPAR(:,6);
areaMedian=nanmedian(datainFAPAR(:,2:end),2);



% for each time entry find the associated fAPAR
for ii=1:length(fluxTime)
    
    cYear=fluxTime(ii,1);
    cDOY=fluxTime(ii,2);
   
    indX=datainFAPAR(:,1)==cYear & datainFAPAR(:,2)==cDOY;
    if sum(indX)>0
    fAPAR(ii)=siteFAPAR(indX);
    else 
        fAPAR(ii)=NaN;
    end
    
end

% interpolate to get continuous fAPAR
indX=~isnan(fAPAR);
[~, b]=find(fAPAR>=0);
vq = interp1(b,fAPAR(~isnan(fAPAR)),1:length(fAPAR),'linear','extrap');
smoothedFAPAR=smooth(vq,48*14);
%[tmpDataSmooth,r1,vr]=ssa_tk2(vq,48*14,0); too slow to be worth implementing
% delete all 'smoothed' entries up to the first valid data point
smoothedFAPAR(1:b)= NaN;

siteFAPAR=smoothedFAPAR;

end

