% this script will...

% load the parameters estimated by REddyProc
% compare RrefDay with RrefNight
% and extract monthly values for plotting later

% now using only growing season data for calculating % inhibition

% the results here are not dependent on whether the input is with or
% without inhibition, as it just compares R_night to R_day

% T. Keenan, November 2018

clearvars

noFAPAR=0;      % set to 1 if not using FAPAR
saveOutputData=1;

airTsoilT='airT'; % options, 'airT','soilT' define the temperature used for partitioning

addpath('./functions')

% get the list of sites for which data has been processed
homeREddyProcOut = '../data_REddyProcOutput/';
d = dir(homeREddyProcOut);
nameFolds = {d.name}';
% remove the non-sites '.' and '..'
nameFolds(ismember(nameFolds,{'.','..'})) = [];
% extract only individual sites
IndexC = strfind(nameFolds, '_parameters.');
Index = find(not(cellfun('isempty', IndexC)));
processedSites=nameFolds(Index);

sites=processedSites;
clear Index IndexC isub d nameFolds

% load the tier 1 tier 2 designation (Tier 1 open, Tier 2 not)
tmp=readtable('../data_FLUXNET2015release3/data-availability-openOrNot_v3.csv','Delimiter',',');
openOrNot=table2array(tmp);
indX=ismember(openOrNot,'Tier 1');
yearzOpenOrNot=1991:2014;

% loop through each site
for ii =1:length(sites) % harvard forest is site # 128
    close all
    cSite=sites(ii);
    cSiteShort=cSite{1}(1:6);
    disp(cSiteShort)
    disp(airTsoilT)
    
    if strcmp(airTsoilT,'airT')
        param=readtable(strcat(homeREddyProcOut,cSiteShort,'REddyProc_parameters.csv'),'TreatAsEmpty','NA');
    else
        param= readtable(strcat(homeREddyProcOut,cSiteShort,'REddyProc_parameters4_dev2.csv'),'TreatAsEmpty','NA');
    end
    DT_GPP=csvread(strcat(homeREddyProcOut,cSiteShort,'_GPP_DT_VUT_USTAR50.csv'),1,1);
    
    numHourz=length(unique(param.Hour));
    GPP95 = running_percentile(DT_GPP(:,2),numHourz*2,95); % running 5-day 95th percentile
    
    %%
    if noFAPAR==0
        %% load the fAPAR data
        % for ii =1:length(sites)
        homeFAPAR='../modis-land-product-subset-master/Rcode/formattedData/';
        filename=strcat(homeFAPAR, cSiteShort,'_MOD15A2_Fpar_1km.csv');
        try
            fAPARData=readtable(filename);
            fAPAR=fAPARData.R1C1;

            % remove the first year fAPAR
            indX=fAPARData.Year==2000;
            fAPARData.R1C1(indX)=NaN;
        catch
            fAPARData.R1C1=fAPARData.R1C1+NaN;
        end
    end
    
    %%
    % remove any data that is Tier 2
    % find the openOrNot vector for the current site
    indX=ismember(openOrNot,cSiteShort);
    cSiteOpenOrNot=openOrNot(indX(:,1),2:end);
    % find tier 2 years to delete
    yearzNotOpen = yearzOpenOrNot(ismember(cSiteOpenOrNot,'Tier 2'));
    % if there are years that are not open, delete them
    if ~isempty(yearzNotOpen)
        indXTier2=ismember(param.Year,yearzNotOpen);
    end
    
    for ii = 7:8
        switch ii
            case 7
                current='R_night';
                cTitle='R Night';
            case 8
                current='R_ref';
                cTitle='R Day';
        end
        
        % LOW GPP FILTER!!!!
        % remove R_night and R_ref where there is no GPP        
        Y95 = prctile( GPP95 , 95 );
        Y05 = prctile( GPP95 , 05 );
        threshold = Y05+0.1*(Y95-Y05);
        y=param.(current);
    
        y(GPP95<threshold)=NaN;
        param.(current)=y;
        
        y=param.(current);
        y(y==-9999)=NaN;
        param.(current)=y;
        
        % filter for Tier 2 data
        if ~isempty(yearzNotOpen)
            y=param.(current);
            y(indXTier2)=NaN;
            param.(current)=y;
        end
    
    end
    
    % get seasonal patterns of Rref and Rref_Night
    %% mean monthly values for each
    
    paramNo20112012 = param;
    if strcmp(cSiteShort,'US-Ha1')
        indX=param.Year==2010 | param.Year==2011;
        paramNo20112012(indX,:)=[];
        GPP95(indX)=[];
    end
    
    % clean the param file by removing large outliers
    paramNo20112012.R_night(paramNo20112012.R_night>40)=NaN;
    paramNo20112012.R_ref(paramNo20112012.R_ref>40)=NaN;
    paramNo20112012.R_night(paramNo20112012.R_night<=0)=NaN;
    paramNo20112012.R_ref(paramNo20112012.R_ref<=0)=NaN;
    
    paramNo20112012.percentInhibition=100*(paramNo20112012.R_night-paramNo20112012.R_ref)./paramNo20112012.R_night;
    
    % clean the percent inhibition by removing outliers
    paramNo20112012.percentInhibition(paramNo20112012.percentInhibition<-200)=NaN;
    paramNo20112012.percentInhibition(paramNo20112012.percentInhibition>200)=NaN;
    paramNo20112012.percentInhibition(~isfinite(paramNo20112012.percentInhibition))=NaN;
    
    % get the mean monthly parameters for each year
    test= accumarray([paramNo20112012.Year,paramNo20112012.Month],paramNo20112012.R_ref,[],@nanmedian,NaN);
    monthlyParam.Rref2=nanmean(test);
    monthlyParam.RrefStd2=nanstd(test);
    test= accumarray([paramNo20112012.Year,paramNo20112012.Month],paramNo20112012.R_night,[],@nanmedian,NaN);
    monthlyParam.RrefNight2=nanmean(test);
    monthlyParam.RrefNightStd2=nanstd(test);
    
    test= accumarray([paramNo20112012.Year,paramNo20112012.Month],paramNo20112012.percentInhibition,[],@nanmedian,NaN);
    monthlyParam.percentInh2=nanmean(test);
    monthlyParam.percentInhStd2=nanstd(test);
    
    % get the mean & stdev monthly value over all years
    % monthlyParamz.Rref=nanmean(monthlyParam.Rref2);
    monthlyParam.Rref= accumarray(paramNo20112012.Month,paramNo20112012.R_ref,[],@nanmedian,NaN);
    monthlyParam.RrefStd= accumarray(paramNo20112012.Month,paramNo20112012.R_ref,[],@nanstd,NaN);
    monthlyParam.RrefNight= accumarray(paramNo20112012.Month,paramNo20112012.R_night,[],@nanmedian,NaN);
    monthlyParam.RrefNightStd= accumarray(paramNo20112012.Month,paramNo20112012.R_night,[],@nanstd,NaN);
    monthlyParam.percentInh= accumarray(paramNo20112012.Month,paramNo20112012.percentInhibition,[],@nanmedian,NaN);
    monthlyParam.percentInhStd= accumarray(paramNo20112012.Month,paramNo20112012.percentInhibition,[],@nanstd,NaN);
    
    % get the monthly 95th percentile GPP
    monthlyParam.GPP95= accumarray(paramNo20112012.Month,GPP95,[],@nanmean,NaN);
    
    monthlyParam.GPP95_growingSeason=monthlyParam.GPP95;
    monthlyParam.GPP95_growingSeason(monthlyParam.GPP95_growingSeason<1)=NaN;
    
    yearz = fAPARData.Year;
    doy = fAPARData.Doy;
    [~, mm, ~, ~, ~] = datevec(datenum(yearz,1,doy));
    
    % get the mean fAPAR from the 3x3 window
    % and replace 255s with NaN
    tmp=fAPARData{:,5:end};
    tmp(tmp==250)=NaN;
    medianFAPARdata=nanmedian(tmp,2);

    if noFAPAR==0
        monthlyParam.fAPARmean= accumarray(mm,medianFAPARdata,[],@nanmedian,NaN);
        monthlyParam.fAPARstd= accumarray(mm,medianFAPARdata,[],@nanstd,NaN);
    end
    
    numYears=length(unique(paramNo20112012.Year))-length(yearzNotOpen);
    % to avoid divide by zero issues, set zero numyears to 1.1
    % num years is only used for plotting
    % and removing years with less than 5 years
    % so setting to 1.1 is ok.
    if numYears == 0 
        numYears = 1.1;
    end
    
    % save a datafile for each site with the final figure data
    outFileName=strcat('./data_meanMonthly/Tier1only_gs95/',cSiteShort,'_',airTsoilT);
    
    % remove no data's to reduce file size
    monthly.RefNight=monthlyParam.RrefNight2;
    monthly.RefNightSte =monthlyParam.RrefNightStd2/(numYears^0.5);
    monthly.RefDay=monthlyParam.Rref2;
    monthly.RefDaySte =monthlyParam.RrefStd2/(numYears^0.5);
    monthly.inhibition=monthlyParam.percentInh2;
    monthly.inhibitionSte=monthlyParam.percentInhStd2/(numYears^0.5);
    outData=horzcat([1:12]',monthly.RefNight',monthly.RefNightSte',...
        monthly.RefDay',monthly.RefDaySte',monthly.inhibition',...
        monthly.inhibitionSte',numYears*ones(1,12)',...
        monthlyParam.fAPARmean, monthlyParam.fAPARstd/(numYears^0.5));
   
    T = array2table(outData,...
        'VariableNames',{'Month','RrefNight','RrefNightSte',...
        'RrefDay','RrefDaySte',...
        'Inhibition','InhibitionSte','numYears'...
        'fAPAR','fAPARste',...
        });
    if saveOutputData == 1
        save(outFileName,'T')
    end
    
end
