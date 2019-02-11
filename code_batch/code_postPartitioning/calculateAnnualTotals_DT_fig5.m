% this script will...
% load the post-processed FLUXNET2015-REddyProc file

% get annual sums and compare

% T. Keenan, November 2018

clear all

addpath('./functions')

saveFigures=1;
loadData=0; %1 loads original, 0 loads preprocessed (by this script with loadData=0, for speed)

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

annualTotalsReco_gs_nt=nan(2016,length(processedSites));
annualTotalsReco_gs_nt_wInhib=nan(2016,length(processedSites));

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
    DT_GPP=csvread(strcat(homeREddyProcOut,cSiteShort,'_GPP_DT_VUT_USTAR50.csv'),1,1);
    % check size is correct
    sizex=size(DT_GPP);
    if sizex(2)<2
        DT_GPP=zeros(sizex(1),2);
    end
    DT_GPP_wInhib=csvread(strcat(homeREddyProcOut,cSiteShort,'_GPP_DT_VUT_USTAR50_wInhib.csv'),1,1);
    % check size is correct
    sizex=size(DT_GPP_wInhib);
    if sizex(2)<2
        DT_GPP_wInhib=zeros(sizex(1),2);
    end
    % load the NT data for with and without inhibition
    DT_Reco=csvread(strcat(homeREddyProcOut,cSiteShort,'_Reco_DT_VUT_USTAR50.csv'),1,1);
    % check size is correct
    sizex=size(DT_Reco);
    if sizex(2)<1
        DT_Reco=zeros(sizex(1),1);
    end
    DT_Reco_wInhib=csvread(strcat(homeREddyProcOut,cSiteShort,'_Reco_DT_VUT_USTAR50_wInhib.csv'),1,1);
    % check size is correct
    sizex=size(DT_Reco_wInhib);
    if sizex(2)<1
        DT_Reco_wInhib=zeros(sizex(1),1);
    end
    
    % get annual totals for growing season only
    % need to load the param file to get the hours for this site
    param=readtable(strcat(homeREddyProcOut,cSiteShort,'REddyProc_parameters.csv'),'TreatAsEmpty','NA');
   
    numHourz=length(unique(param.Hour));
    if numHourz==48
         conversionT=conversion/2; % values scaled to per half hour before summing
    else
          conversionT=conversion;
    end

    % get annual GPP for both with and without inhibition
    tmp= accumarray(DT_GPP(:,1),DT_GPP(:,2),[],@nansum,NaN);
    annualTotalsGPP(1:length(tmp),ii)=tmp*conversionT;
    tmp= accumarray(DT_GPP(:,1),DT_GPP_wInhib(:,2),[],@nansum,NaN);
    annualTotalsGPP_wInhib(1:length(tmp),ii)=tmp*conversionT;
    
    % Reco
    tmp= accumarray(DT_GPP(:,1),DT_Reco(:,1),[],@nansum,NaN);
    annualTotalsReco(1:length(tmp),ii)=tmp*conversionT;
    tmp= accumarray(DT_GPP(:,1),DT_Reco_wInhib(:,1),[],@nansum,NaN);
    annualTotalsReco_wInhib(1:length(tmp),ii)=tmp*conversionT;
    
    GPP95 = running_percentile(DT_GPP(:,2),numHourz*2,95); % running 5-day 95th percentile
    
    % set off season times to NAN
    indX=GPP95<2;
    
    % set day-times to NAN
    indXnight=DT_GPP(:,2)>0.5;
    
    % delete all the non - growing season data
    GPP_gs = DT_GPP(:,2);
    GPP_gs_wInhib = DT_GPP_wInhib(:,2);
    GPP_gs(indX)=NaN;
    GPP_gs_wInhib(indX)=NaN;
        
    Reco_gs = DT_Reco(:,1);
    Reco_gs_wInhib = DT_Reco_wInhib(:,1);
    Reco_gs(indX)=NaN;
    Reco_gs_wInhib(indX)=NaN;
    
    % delete all the non - growing season day-time data
    Reco_gs_nt = Reco_gs;
    Reco_gs_nt_wInhib = Reco_gs_wInhib;
    Reco_gs_nt(indXnight)=NaN;
    Reco_gs_nt_wInhib(indXnight)=NaN;
    
     % get annual growing season GPP for both with and without inhibition
    tmp= accumarray(DT_GPP(:,1),GPP_gs(:,1),[],@nansum,NaN)*conversionT;
    annualTotalsGPP_gs(1:length(tmp),ii)=tmp;
    tmp= accumarray(DT_GPP(:,1),GPP_gs_wInhib(:,1),[],@nansum,NaN)*conversionT;
    annualTotalsGPP_gs_wInhib(1:length(tmp),ii)=tmp;
    
    % get annual growing season Reco
    tmp= accumarray(DT_GPP(:,1),Reco_gs(:,1),[],@nansum,NaN)*conversionT;
    annualTotalsReco_gs(1:length(tmp),ii)=tmp;
    tmp= accumarray(DT_GPP(:,1),Reco_gs_wInhib(:,1),[],@nansum,NaN)*conversionT;
    annualTotalsReco_gs_wInhib(1:length(tmp),ii)=tmp;
 
    % get annual growing season night time Reco
    tmp= accumarray(DT_GPP(:,1),Reco_gs_nt(:,1),[],@nansum,NaN)*conversionT;
    annualTotalsReco_gs_nt(1:length(tmp),ii)=tmp;
    tmp= accumarray(DT_GPP(:,1),Reco_gs_nt_wInhib(:,1),[],@nansum,NaN)*conversionT;
    annualTotalsReco_gs_nt_wInhib(1:length(tmp),ii)=tmp;
    
end

    % set infinite values to NaN
    annualTotalsGPP(isinf(annualTotalsGPP))=NaN;
    annualTotalsGPP_wInhib(isinf(annualTotalsGPP_wInhib))=NaN;
    annualTotalsGPP_gs(isinf(annualTotalsGPP_gs))=NaN;
    annualTotalsGPP_gs_wInhib(isinf(annualTotalsGPP_gs_wInhib))=NaN;

    annualTotalsReco(isinf(annualTotalsReco))=NaN;
    annualTotalsReco_wInhib(isinf(annualTotalsReco_wInhib))=NaN;
    annualTotalsReco_gs(isinf(annualTotalsReco_gs))=NaN;
    annualTotalsReco_gs_wInhib(isinf(annualTotalsReco_gs_wInhib))=NaN;
    annualTotalsReco_gs_nt(isinf(annualTotalsReco_gs_nt))=NaN;
    annualTotalsReco_gs_nt_wInhib(isinf(annualTotalsReco_gs_nt_wInhib))=NaN;

    % set negative values to NaN
    annualTotalsGPP(annualTotalsGPP<0)=NaN;
    annualTotalsGPP_wInhib(annualTotalsGPP_wInhib<0)=NaN;
    annualTotalsGPP_gs(annualTotalsGPP_gs<0)=NaN;
    annualTotalsGPP_gs_wInhib(annualTotalsGPP_gs_wInhib<0)=NaN;
    
    annualTotalsReco(annualTotalsReco>4000)=NaN;
    annualTotalsReco(annualTotalsReco>4000)=NaN;
    annualTotalsReco_wInhib(annualTotalsReco_wInhib>4000)=NaN;
    annualTotalsReco_gs(annualTotalsReco_gs>4000)=NaN;
    annualTotalsReco_gs_wInhib(annualTotalsReco_gs_wInhib>4000)=NaN;
    annualTotalsReco_gs_nt(annualTotalsReco_gs_nt>4000)=NaN;
    annualTotalsReco_gs_nt_wInhib(annualTotalsReco_gs_nt_wInhib>4000)=NaN;

    annualTotalsReco(annualTotalsReco<=0)=NaN;
    annualTotalsReco_wInhib(annualTotalsReco_wInhib<=0)=NaN;
    annualTotalsReco_gs(annualTotalsReco_gs<=0)=NaN;
    annualTotalsReco_gs_wInhib(annualTotalsReco_gs_wInhib<=0)=NaN;
    annualTotalsReco_gs_nt(annualTotalsReco_gs_nt<=0)=NaN;
    annualTotalsReco_gs_nt_wInhib(annualTotalsReco_gs_nt_wInhib<=0)=NaN;
  
end

if loadData==1
    save('./data_inter/dataFromCompAnnTotsDT_GPP.mat','annualTotalsGPP','annualTotalsGPP_gs',...
        'annualTotalsGPP_gs_wInhib','annualTotalsGPP_wInhib')
    save('./data_inter/dataFromCompAnnTotsDT_Reco.mat','annualTotalsReco','annualTotalsReco_gs','annualTotalsReco_gs_nt_wInhib',...
        'annualTotalsReco_gs_wInhib','annualTotalsReco_wInhib','annualTotalsReco_gs_nt')
end
if ~exist('annualTotalsGPP','var')
    load('./data_inter/dataFromCompAnnTotsDT_GPP.mat')
    load('./data_inter/dataFromCompAnnTotsDT_Reco.mat')
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Calculate the bias associated with each period

% for GPP
xdata=annualTotalsGPP(:);
ydata=annualTotalsGPP_wInhib(:);

tmp=(nanmean(xdata)-nanmean(ydata))./nanmean(xdata);
bias.GPP(1)=nanmean(tmp);
tmp=(xdata-ydata)./xdata;
tmp(tmp<-1)=NaN;
bias.GPP2_mean(1)=nanmean(tmp);
bias.GPP2_median(1)=nanmedian(tmp);
bias.GPP2_std=nanstd(tmp);


% for Reco (all)
xdata=annualTotalsReco(:);
ydata=annualTotalsReco_wInhib(:);


%%

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


% for Reco (gs nt)
xdata=annualTotalsReco_gs_nt(:);
ydata=annualTotalsReco_gs_nt_wInhib(:);

tmp2=(nanmean(xdata)-nanmean(ydata))/nanmean(xdata);
bias.RecoGsNt(1)=tmp2;
tmp=(xdata-ydata)./xdata;
tmp(tmp<-1)=NaN;
bias.Reco2gsnt_mean(1)=nanmean(tmp);
bias.Reco2gsnt_median(1)=nanmedian(tmp);
bias.Reco2gsnt_std=nanstd(tmp);

bias.count=sum(~isnan(tmp));

save('./data_inter/bias_DT.mat','bias')



