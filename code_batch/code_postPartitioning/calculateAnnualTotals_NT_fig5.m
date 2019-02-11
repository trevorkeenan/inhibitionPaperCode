% this script will...
% load the post-processed FLUXNET2015-REddyProc file

% get annual sums and compare

% T. Keenan, November 2018

clear all

addpath('./functions')

saveFigures=0;
loadData=1; %1 loads original, 0 loads preprocessed (by this script with loadData=0)

if loadData==1
    clearvars -except loadData saveFigures
end

convertUmolCO2togC=(12/44)*1000000;
convertPerStoPerH= 3600;
conversion=3.2757*convertPerStoPerH/convertUmolCO2togC;

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

if loadData==1
    %%
    % loop through each processed site, find F15 site data and merge
    annualTotalsGPP=nan(2016,length(processedSites));
    annualTotalsGPP_wInhib=nan(2016,length(processedSites));
    annualTotalsReco=nan(2016,length(processedSites));
    annualTotalsReco_wInhib=nan(2016,length(processedSites));
    
    annualTotalsGPP_gs=nan(2016,length(processedSites));
    annualTotalsGPP_gs_wInhib=nan(2016,length(processedSites));
    annualTotalsReco_gs=nan(2016,length(processedSites));
    annualTotalsReco_gs_wInhib=nan(2016,length(processedSites));
    
    annualTotalsReco_gs_dt=nan(2016,length(processedSites));
    annualTotalsReco_gs_dt_wInhib=nan(2016,length(processedSites));
    annualTotalsGPP_gs_dt=nan(2016,length(processedSites));
    annualTotalsGPP_gs_dt_wInhib=nan(2016,length(processedSites));
    
    for ii =1:length(sites)
        
        close all
        cSite=sites(ii);
        cSiteShort=cSite{1}(1:6);
        disp(cSiteShort)
        
        % set indices of the data of interest
        indReco=6;
        indGPP=7;
        indReco_wInhib=7;
        indGPP_wInhib=8;
        
        % load the NT data for with and without inhibition
        NT=csvread(strcat(homeREddyProcOut,cSiteShort,'REddyProc_NT_VUT_USTAR50.csv'),1,1);
        % check size is correct
        sizex=size(NT);
        if sizex(2)<7
            NT=zeros(sizex(1),7);
        end
        NT_wInhib=csvread(strcat(homeREddyProcOut,cSiteShort,'REddyProc_NT_VUT_USTAR50_wInhib.csv'),1,1);
        % check size is correct
        sizex=size(NT_wInhib);
        if sizex(2)<8
            NT_wInhib=zeros(sizex(1),8);
        end
        % the param file has the year etc.
        param= readtable(strcat(homeREddyProcOut,cSiteShort,'REddyProc_parameters.csv'),'TreatAsEmpty','NA');
        
             numHourz=length(unique(param.Hour));
     if numHourz==48
         conversionT=conversion/2; % values scaled to per half hour before summing
     else
          conversionT=conversion;
     end
     
        % get annual GPP for both with and without inhibition
        tmp= accumarray(param.Year,NT(:,indGPP),[],@nansum,NaN);
        annualTotalsGPP(1:length(tmp),ii)=tmp*conversionT;
        tmp= accumarray(param.Year,NT_wInhib(:,indGPP_wInhib),[],@nansum,NaN);
        annualTotalsGPP_wInhib(1:length(tmp),ii)=tmp*conversionT;
        
        % Reco
        tmp= accumarray(param.Year,NT(:,indReco),[],@nansum,NaN);
        annualTotalsReco(1:length(tmp),ii)=tmp*conversionT;
        tmp= accumarray(param.Year,NT_wInhib(:,indReco_wInhib),[],@nansum,NaN);
        annualTotalsReco_wInhib(1:length(tmp),ii)=tmp*conversionT;
        
        numHourz=length(unique(param.Hour));
        GPP95 = running_percentile(NT(:,indGPP),numHourz*2,95); % running 5-day 95th percentile
        
        % find off season times to NAN
        indX=GPP95<2;
        
        % find day-times to NAN (here set night-time data to NAN)
        indXnight=NT(:,indGPP)<0.5;
        
        % delete all the non - growing season data
        GPP_gs = NT(:,indGPP);
        GPP_gs_wInhib = NT_wInhib(:,indGPP_wInhib);
        GPP_gs(indX)=NaN;
        GPP_gs_wInhib(indX)=NaN;
        
        % delete all the non - growing season night-time data
        GPP_gs_dt = GPP_gs;
        GPP_gs_dt_wInhib = GPP_gs_wInhib;
        GPP_gs_dt(indXnight)=NaN;
        GPP_gs_dt_wInhib(indXnight)=NaN;
        
        Reco_gs = NT(:,indReco);
        Reco_gs_wInhib = NT_wInhib(:,indReco_wInhib);
        Reco_gs(indX)=NaN;
        Reco_gs_wInhib(indX)=NaN;
        
        % delete all the non - growing season night-time data
        Reco_gs_dt = Reco_gs;
        Reco_gs_dt_wInhib = Reco_gs_wInhib;
        Reco_gs_dt(indXnight)=NaN;
        Reco_gs_dt_wInhib(indXnight)=NaN;
        
        % get annual growing season GPP for both with and without inhibition
        tmp= accumarray(param.Year,GPP_gs(:,1),[],@nansum,NaN)*conversionT;
        annualTotalsGPP_gs(1:length(tmp),ii)=tmp;
        tmp= accumarray(param.Year,GPP_gs_wInhib(:,1),[],@nansum,NaN)*conversionT;
        annualTotalsGPP_gs_wInhib(1:length(tmp),ii)=tmp;
        
        % get annual growing season DAYTIME GPP for both with and without inhibition
        tmp= accumarray(param.Year,GPP_gs_dt(:,1),[],@nansum,NaN)*conversionT;
        annualTotalsGPP_gs_dt(1:length(tmp),ii)=tmp;
        tmp= accumarray(param.Year,GPP_gs_dt_wInhib(:,1),[],@nansum,NaN)*conversionT;
        annualTotalsGPP_gs_dt_wInhib(1:length(tmp),ii)=tmp;
        
        
        % get annual growing season Reco
        tmp= accumarray(param.Year,Reco_gs(:,1),[],@nansum,NaN)*conversionT;
        annualTotalsReco_gs(1:length(tmp),ii)=tmp;
        tmp= accumarray(param.Year,Reco_gs_wInhib(:,1),[],@nansum,NaN)*conversionT;
        annualTotalsReco_gs_wInhib(1:length(tmp),ii)=tmp;
        
        % get annual growing season night time Reco
        tmp= accumarray(param.Year,Reco_gs_dt(:,1),[],@nansum,NaN)*conversionT;
        annualTotalsReco_gs_dt(1:length(tmp),ii)=tmp;
        tmp= accumarray(param.Year,Reco_gs_dt_wInhib(:,1),[],@nansum,NaN)*conversionT;
        annualTotalsReco_gs_dt_wInhib(1:length(tmp),ii)=tmp;
        
    end
    
    
    % set infinite values to NaN
    annualTotalsGPP(isinf(annualTotalsGPP))=NaN;
    annualTotalsGPP_wInhib(isinf(annualTotalsGPP_wInhib))=NaN;
    annualTotalsGPP_gs(isinf(annualTotalsGPP_gs))=NaN;
    annualTotalsGPP_gs_wInhib(isinf(annualTotalsGPP_gs_wInhib))=NaN;
    annualTotalsGPP_gs_dt(isinf(annualTotalsGPP_gs_dt))=NaN;
    annualTotalsGPP_gs_dt_wInhib(isinf(annualTotalsGPP_gs_wInhib))=NaN;
    
    annualTotalsGPP(annualTotalsGPP>4000)=NaN;
    annualTotalsGPP_wInhib(annualTotalsGPP_wInhib>4000)=NaN;
    annualTotalsGPP_gs(annualTotalsGPP_gs>4000)=NaN;
    annualTotalsGPP_gs_wInhib(annualTotalsGPP_gs_wInhib>4000)=NaN;
    annualTotalsGPP_gs_dt(annualTotalsGPP_gs_dt>4000)=NaN;
    annualTotalsGPP_gs_dt_wInhib(annualTotalsGPP_gs_dt_wInhib>4000)=NaN;
    
    annualTotalsReco(isinf(annualTotalsReco))=NaN;
    annualTotalsReco_wInhib(isinf(annualTotalsReco_wInhib))=NaN;
    annualTotalsReco_gs(isinf(annualTotalsReco_gs))=NaN;
    annualTotalsReco_gs_wInhib(isinf(annualTotalsReco_gs_wInhib))=NaN;
    annualTotalsReco_gs_dt(isinf(annualTotalsReco_gs_dt))=NaN;
    annualTotalsReco_gs_dt_wInhib(isinf(annualTotalsReco_gs_dt_wInhib))=NaN;
    
    annualTotalsReco(annualTotalsReco>4000)=NaN;
    annualTotalsReco(annualTotalsReco>4000)=NaN;
    annualTotalsReco_wInhib(annualTotalsReco_wInhib>4000)=NaN;
    annualTotalsReco_gs(annualTotalsReco_gs>4000)=NaN;
    annualTotalsReco_gs_wInhib(annualTotalsReco_gs_wInhib>4000)=NaN;
    annualTotalsReco_gs_dt(annualTotalsReco_gs_dt>4000)=NaN;
    annualTotalsReco_gs_dt_wInhib(annualTotalsReco_gs_dt_wInhib>4000)=NaN;
    
    % set negative values to NaN
    annualTotalsGPP(annualTotalsGPP<=0)=NaN;
    annualTotalsGPP_wInhib(annualTotalsGPP_wInhib<=0)=NaN;
    annualTotalsGPP_gs(annualTotalsGPP_gs<=0)=NaN;
    annualTotalsGPP_gs_wInhib(annualTotalsGPP_gs_wInhib<=0)=NaN;
    annualTotalsGPP_gs_dt(annualTotalsGPP_gs_dt<=0)=NaN;
    annualTotalsGPP_gs_dt_wInhib(annualTotalsGPP_gs_dt_wInhib<=0)=NaN;
    
    annualTotalsReco(annualTotalsReco<=0)=NaN;
    annualTotalsReco_wInhib(annualTotalsReco_wInhib<=0)=NaN;
    annualTotalsReco_gs(annualTotalsReco_gs<=0)=NaN;
    annualTotalsReco_gs_wInhib(annualTotalsReco_gs_wInhib<=0)=NaN;
    annualTotalsReco_gs_dt(annualTotalsReco_gs_dt<=0)=NaN;
    annualTotalsReco_gs_dt_wInhib(annualTotalsReco_gs_dt_wInhib<=0)=NaN;
    
    
end

if loadData==1
    save('./data_inter/dataFromCompAnnTotsNT_GPP.mat','annualTotalsGPP','annualTotalsGPP_gs','annualTotalsGPP_gs_dt',...
        'annualTotalsGPP_gs_wInhib','annualTotalsGPP_wInhib','annualTotalsGPP_gs_dt_wInhib')
    save('./data_inter/dataFromCompAnnTotsNT_Reco.mat','annualTotalsReco','annualTotalsReco_gs','annualTotalsReco_gs_dt_wInhib',...
        'annualTotalsReco_gs_wInhib','annualTotalsReco_wInhib','annualTotalsReco_gs_dt')
end
if ~exist('annualTotalsGPP','var')
    load('./data_inter/dataFromCompAnnTotsNT_GPP.mat')
    load('./data_inter/dataFromCompAnnTotsNT_Reco.mat')
end


%%
%   Calculate the bias associated with each period

% for GPP
xdata=annualTotalsGPP(:);
ydata=annualTotalsGPP_wInhib(:);

tmp=(nanmean(xdata)-nanmean(ydata))./nanmean(xdata);
bias.GPP1mean(1)=nanmean(tmp);
tmp=(xdata-ydata)./xdata;
tmp(tmp<-1)=NaN;
bias.GPP2_mean(1)=nanmean(tmp);
bias.GPP2_median(1)=nanmedian(tmp);
bias.GPP2_std=nanstd(tmp);

% for GPP during growing season
xdata=annualTotalsGPP_gs(:);
ydata=annualTotalsGPP_gs_wInhib(:);

tmp=(nanmean(xdata)-nanmean(ydata))./nanmean(xdata);
bias.GPPgs(1)=nanmean(tmp);
tmp=(xdata-ydata)./xdata;
tmp(tmp<-1)=NaN;
bias.GPP2gs_mean(1)=nanmean(tmp);
bias.GPP2gs_median(1)=nanmedian(tmp);
bias.GPP2gs_std=nanstd(tmp);


% for GPP during growing season and daytime
xdata=annualTotalsGPP_gs_dt(:);
ydata=annualTotalsGPP_gs_dt_wInhib(:);

tmp=(nanmean(xdata)-nanmean(ydata))./nanmean(xdata);
bias.GPPgsdt(1)=nanmean(tmp);
tmp=(xdata-ydata)./xdata;
tmp(tmp<-1)=NaN;
bias.GPP2gsdt_mean(1)=nanmean(tmp);
bias.GPP2gsdt_median(1)=nanmedian(tmp);
bias.GPP2gsdt_std=nanstd(tmp);


% for Reco (all)
xdata=annualTotalsReco(:);
ydata=annualTotalsReco_wInhib(:);

tmp=(xdata-ydata)./xdata;
tmp2=(nanmean(xdata)-nanmean(ydata))/nanmean(xdata);
bias.Reco(1)=tmp2;
tmp=(xdata-ydata)./xdata;
tmp(tmp<-1)=NaN;
bias.Reco2_mean(1)=nanmean(tmp);
bias.Reco2_median(1)=nanmedian(tmp);
bias.Reco2_std=nanstd(tmp);

% for Reco (gs)
xdata=annualTotalsReco_gs(:);
ydata=annualTotalsReco_gs_wInhib(:);

tmp2=(nanmean(xdata)-nanmean(ydata))/nanmean(xdata);
bias.RecoGs(1)=tmp2;
tmp=(xdata-ydata)./xdata;
tmp(tmp<-1)=NaN;
bias.Reco2gs_mean(1)=nanmean(tmp);
bias.Reco2gs_median(1)=nanmedian(tmp);
bias.Reco2gs_std=nanstd(tmp);


% for Reco (gs dt)
xdata=annualTotalsReco_gs_dt(:);
ydata=annualTotalsReco_gs_dt_wInhib(:);

tmp2=(nanmean(xdata)-nanmean(ydata))/nanmean(xdata);
bias.RecoGsDt(1)=tmp2;
tmp=(xdata-ydata)./xdata;
tmp(tmp<-1)=NaN;
bias.Reco2gsdt_mean(1)=nanmean(tmp);
bias.Reco2gsdt_median(1)=nanmedian(tmp);
bias.Reco2gsdt_std=nanstd(tmp);

bias.count=sum(~isnan(tmp));

save('./data_inter/bias_NT.mat','bias')


