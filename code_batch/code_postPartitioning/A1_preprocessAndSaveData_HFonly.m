% this script will...
% load FLUXNET2015 csv files
% the REddyProc partitioned data
% and Rick Wehr's data

% and save a file for later use containing all data together

% requires functions:
%   date2doy.m
%   getFAPAR.m

% for Harvard forest only

% T. Keenan, November 2018

clearvars

saveData4diurnalCycleAnalysis=1;

addpath('./functions')

home1='../data_FLUXNET2015release3/FLUXNET2015data/';

% get the list of sites
d = dir(home1);
isub = [d(:).isdir]; % returns logical vector
nameFolds = {d(isub).name}';
% remove the non-sites '.' and '..'
nameFolds(ismember(nameFolds,{'.','..'})) = [];
% extract only the FULLSET site data
IndexC = strfind(nameFolds, 'FULLSET');
Index = find(not(cellfun('isempty', IndexC)));

sites=nameFolds(Index);
clear Index IndexC isub d nameFolds


indxHarvard=159;	% index changes depending on how many sites were processed

ii =indxHarvard; 	% could loop through each site by changing this to a for loop
outFileData=[];
close all

cSite=sites(ii);
cSiteShort=cSite{1}(5:10);
disp(cSiteShort)

home=strcat(home1, cSite);
filename1=strcat(home,'/*FULLSET_HH*.csv');
filename=dir(filename1{:});

if isempty(filename)
    filename1=strcat(home,'/*FULLSET_HR*.csv');
    filename=dir(filename1{:});
end

% load the fluxnet processed data
tmp=strcat(home,'/',filename.name);
hourlyData=importdata(tmp{:});

% 1. create year, DOY vector from Timestamp of flux data
tmp2=floor(hourlyData.data(:,1)/10000);
hourz = rem(hourlyData.data(:,1),10000) ;
t = datetime(tmp2,'ConvertFrom','yyyymmdd');
doy=date2doy(datenum(t));
DateVector = datevec(t);
DateVector(:,4)=repmat((1:24)'.*ones(24,1),length(DateVector)/24,1);
fluxTime=[DateVector(:,1),doy];

% get fAPAR for whole time series
% if time series is pre MODIS, no fAPAR is returned
fAPAR=getFAPAR(cSiteShort,fluxTime);

% set the location of data of interest
indNEE=strcmp(hourlyData.colheaders,'NEE_VUT_MEAN');
indNEEqc=strcmp(hourlyData.colheaders,'NEE_VUT_USTAR50_QC');
indGPPnt=strcmp(hourlyData.colheaders,'GPP_NT_VUT_MEAN');
indGPPdt=strcmp(hourlyData.colheaders,'GPP_DT_VUT_MEAN');
indRecont=strcmp(hourlyData.colheaders,'RECO_NT_VUT_MEAN');
indRecodt=strcmp(hourlyData.colheaders,'RECO_DT_VUT_MEAN');
indLE=strcmp(hourlyData.colheaders,'LE_F_MDS');
indSWf=strcmp(hourlyData.colheaders,'SW_IN_F');
indSWpot=strcmp(hourlyData.colheaders,'SW_IN_POT');
indTaf=strcmp(hourlyData.colheaders,'TA_F');
indTafMDS=strcmp(hourlyData.colheaders,'TA_F_MDS');
indTsfMDS=strcmp(hourlyData.colheaders,'TS_F_MDS_1');
indVPD=strcmp(hourlyData.colheaders,'VPD_F');
indNight=strcmp(hourlyData.colheaders,'NIGHT');
indWD = strcmp(hourlyData.colheaders,'WD');

% some data are reported as -9999 (go figure!)
indX = hourlyData.data == -9999;
hourlyData.data(indX)=NaN;

% for Harvard - append REddyProc partitioned data to the file
if ii==indxHarvard
    % load the REddyProc NT
    REddyProcNT=csvread(strcat('../data_REddyProcOutput/',cSiteShort,'REddyProc_NT_VUT_USTAR50.csv'),1,1);
    REddyProcNT(REddyProcNT==-9999)=NaN;
    REddyProcNTReco=REddyProcNT(:,6);
    REddyProcNTGPP=REddyProcNT(:,7);
    REddyProcNT_wInhib=csvread(strcat('../data_REddyProcOutput/',cSiteShort,'REddyProc_NT_VUT_USTAR50_wInhib.csv'),1,1);
    REddyProcNT_wInhib(REddyProcNT_wInhib==-9999)=NaN;
    REddyProcNTReco_wInhib=REddyProcNT_wInhib(:,6);
    REddyProcNTGPP_wInhib=REddyProcNT_wInhib(:,7);
    
    % load the REddyProc DT
    % data has 3 columns
    % new GPP has the year info also
    REddyProcDTGPP=csvread(strcat('../data_REddyProcOutput/',cSiteShort,'_GPP_DT_VUT_USTAR50.csv'),1,1);
    REddyProcDTReco=csvread(strcat('../data_REddyProcOutput/',cSiteShort,'_Reco_DT_VUT_USTAR50.csv'),1,1);
    REddyProcDTGPP_wInhib=csvread(strcat('../data_REddyProcOutput/',cSiteShort,'_GPP_DT_VUT_USTAR50_wInhib.csv'),1,1);
    REddyProcDTReco_wInhib=csvread(strcat('../data_REddyProcOutput/',cSiteShort,'_Reco_DT_VUT_USTAR50_wInhib.csv'),1,1);
    
    % remove years that do not match up from the FLUXNET file
    hourlyYearz=floor(hourlyData.data(:,1)/100000000);
    indX=ismember(hourlyYearz,REddyProcDTGPP(:,1));
    hourlyData.data=hourlyData.data(indX,:);
    % get the unique years
    yearx=floor(hourlyData.data(:,1)./100000000);
    years=unique(yearx);
    tmp2=floor(hourlyData.data(:,1)/10000);
    hourz = rem(hourlyData.data(:,1),10000) ;
    t = datetime(tmp2,'ConvertFrom','yyyymmdd');
    doy=date2doy(datenum(t));
    
    % and trim the associated date vector
    DateVector=DateVector(indX,:);
    
    hourlyDataLength = size(hourlyData.data,2);
    
    hourlyData.data=horzcat(hourlyData.data,...
        REddyProcNTGPP,REddyProcNTReco,REddyProcDTGPP(:,2:end),REddyProcDTReco,...
        REddyProcNTGPP_wInhib,REddyProcNTReco_wInhib,REddyProcDTGPP_wInhib(:,2:end),REddyProcDTReco_wInhib);
    
    indGPPREddyProcNT=hourlyDataLength +1;
    indRecoREddyProcNT=hourlyDataLength +2;
    indGPPREddyProcDT=hourlyDataLength +3;
    indRecoREddyProcDT=hourlyDataLength +4;
    
    indGPPREddyProcNT_wInhib=hourlyDataLength +5;
    indRecoREddyProcNT_wInhib=hourlyDataLength +6;
    indGPPREddyProcDT_wInhib=hourlyDataLength +7;
    indRecoREddyProcDT_wInhib=hourlyDataLength +8;
    
end

% Add the Harvard Forest isotope data
if strcmp(cSiteShort,'US-Ha1')
    % Load Rick Wehr's data
    
    wehrData=readtable('../data_RickWehr/wehrDataHourly.csv','Delimiter',',');
    
    wehrTime=wehrData(:,1:4);
    wehrTimeArray=table2array(wehrTime);
    DateVector2=DateVector(:,1:4);
    indX=wehrTimeArray(:,4)==0;
    wehrTimeArray(indX,4)=24;
    [a,b]= ismember(DateVector2,wehrTimeArray,'rows');
    [a1,b1]= ismember(wehrTimeArray,DateVector2,'rows');
    
    wehrDataArray=table2array(wehrData(:,7:end));
    wehrDataArrayOverlap=wehrDataArray(a1,:);
    wehrDataAligned=nan(length(hourlyData.data),4+width(wehrData)-6);
    wehrDataAligned(1:length(DateVector2),1:4)=DateVector2;
    wehrDataAligned(a,5:end)=wehrDataArrayOverlap;
    
    % log where data will be located
    indxWehrNEE=size(hourlyData.data,2)+1;
    indxWehrGPP_IFP=size(hourlyData.data,2)+2;
    indxWehrReco_IFP=size(hourlyData.data,2)+3;
    indxWehrGPP_NT=size(hourlyData.data,2)+4;
    indxWehrReco_NT=size(hourlyData.data,2)+5;
    indxWehrNEE_DT=size(hourlyData.data,2)+6;
    indxWehrGPP_DT=size(hourlyData.data,2)+7;
    indxWehrReco_DT=size(hourlyData.data,2)+8;
    indxWehrLAI=size(hourlyData.data,2)+9;
    
    % append to hourlyData
    hourlyData.data=horzcat(hourlyData.data,wehrDataAligned(:,5:13));
end


%%
% loop through years to concatenate
outFileData=[];
for jj=1:length(years)
    
    cYear=years(jj);
    disp(cYear)
    
    indX=yearx==cYear;
    cYearData=hourlyData.data(indX,:);
    cYearFAPAR=fAPAR(indX);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Collect the data to write out
    outFileDataYear=[cYearData(:,indNEE),cYearData(:,indGPPdt),cYearData(:,indGPPnt),...
        cYearData(:,indRecodt),cYearData(:,indRecont),...
        cYearData(:,indLE),...
        cYearData(:,indSWf),cYearData(:,indTaf),cYearData(:,indTafMDS),cYearData(:,indTsfMDS),...
        cYearData(:,indVPD),cYearFAPAR,...
        cYearData(:,indGPPREddyProcDT),cYearData(:,indRecoREddyProcDT),... % Original GPP Reco from ReddyPro GL processing
        cYearData(:,indGPPREddyProcNT),cYearData(:,indRecoREddyProcNT),... % Original GPP Reco from ReddyPro GL processing
        cYearData(:,indGPPREddyProcDT_wInhib),cYearData(:,indRecoREddyProcDT_wInhib),... % Original GPP Reco from ReddyPro GL processing
        cYearData(:,indGPPREddyProcNT_wInhib),cYearData(:,indRecoREddyProcNT_wInhib),... % Original GPP Reco from ReddyPro GL processing
        cYearData(:,indxWehrNEE),cYearData(:,indxWehrGPP_IFP),cYearData(:,indxWehrReco_IFP),...  % new GPP
        cYearData(:,indxWehrGPP_NT),cYearData(:,indxWehrReco_NT),cYearData(:,indxWehrNEE_DT),...
        cYearData(:,indxWehrGPP_DT),cYearData(:,indxWehrReco_DT),cYearData(:,indxWehrLAI),...
        cYearData(:,indWD),...
        ];
    
    outFileData=vertcat(outFileData,outFileDataYear);
    
end

% save an output file for further analysis
outFileName2=strcat('./data_intermediateData/',sites(ii));

% remove no data's to reduce file size
outData=horzcat(DateVector(:,1),DateVector(:,2),doy,hourz,outFileData);
indX=isnan(outData(:,1:end-11));
outData(any(indX,2),:) = [];

T = array2table(outData,...
    'VariableNames',{'Year','Month','DOY','HOD','NEE','GPPdt','GPPnt',...
    'Recodt','Recont',...
    'LE','SWin','AirT','AirTMDS','SoilTMDS','VPD','fAPAR',...
    'GPPREddyProcDT','RecoREddyProcDT',...
    'GPPREddyProcNT','RecoREddyProcNT',...
    'GPPREddyProcDT_wInhib','RecoREddyProcDT_wInhib',...
    'GPPREddyProcNT_wInhib','RecoREddyProcNT_wInhib',...
    'wehrNEE','wehrGPP_IFP','wehrReco_IFP',...
    'wehrGPP_NT','wehrReco_NT','wehrNEE_DT',...
    'wehrGPP_DT','wehrReco_DT','LAI',...
    'WD',...
    });
if saveData4diurnalCycleAnalysis == 1
    save(outFileName2{:},'T')
end

