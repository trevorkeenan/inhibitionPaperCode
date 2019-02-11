% this script will...

% read in pre-processed data and calculate the mean diurnal cycle
% in June-July
% and August-September

% currently for Harvard forest only.

% T. Keenan, November 2018

clearvars
% close all

saveFigures=1;
addpath('./functions')

home='./data_intermediateData/';

% get the list of sites
d = dir(home);
sites = {d(:).name}';
% remove the non-sites '.' and '..'
sites(ismember(sites,{'.','..','.DS_Store'})) = [];
clear Index IndexC isub d nameFolds


% define some colores
color.gray=[0.2 0.2 0.2];
color.darkGray=[0.2 0.2 0.2];
color.lightGray=[0.7 0.7 0.7];
color.barkBrown=[0.396 0.263 0.129]; % 39.6% red, 26.3% green and 12.9% blue
color.lightBrown=0.5*[0.396 0.263 0.129];
color.saddleBrown=[139 69 19]/255;
color.darkBrown=[200 69 19]/255;
color.blue=[0 0 1];
color.darkBlue=[0 0 0.5];
color.deepBlue=[0 0 0.2];
color.green=[0 1 0];
color.darkGreen=[0 0.5 0];
color.deepGreen=[0 0.2 0];

% 1. extract june july data - get mean diurnal cycle
% 2. extract august sept. data - get mean diurnal cycle
ii =1; %% could loop through each site


cSite=sites(ii);
cSiteShort=cSite{1}(5:10);
disp(cSite)
disp(cSiteShort)

filename=strcat(home, cSite);

tmp=filename{:};
tmp=tmp(1:end-3);
tmp=strcat(tmp,'mat');
hourlyData2=load(tmp);
hourlyData=hourlyData2.T;
indX=hourlyData.Year<=1992;
hourlyData(indX,:)=[];

hourlyData.HOD=hourlyData.HOD+100;

% Limit to a particular wind direction
%     indX=hourlyData.WD<180 | hourlyData.WD > 270;
%     hourlyData.RecoREddyProcNT(indX)=NaN;
%     hourlyData.RecoREddyProcNT_wInhib(indX)=NaN;
%     hourlyData.RecoREddyProcDT(indX)=NaN;
%     hourlyData.RecoREddyProcDT_wInhib(indX)=NaN;

% get data for June July and August September
indX= hourlyData.DOY> 151 & hourlyData.DOY< 212; % 151 = June 1st; 212 = July 31st
data.juneJuly = hourlyData(indX,:);
indX= hourlyData.DOY> 213 & hourlyData.DOY< 274; % 213 = Aug 1st; 274 = Sept 30st
data.augSept = hourlyData(indX,:);

% get mean diurnal cycles for each
variables=data.juneJuly.Properties.VariableNames; %get variable names
for jj=4:length(variables)
    tmp1(:,jj)= accumarray(data.juneJuly.HOD/100,data.juneJuly.(variables{jj}),[],@nanmedian,NaN);
    tmp2(:,jj)= accumarray(data.augSept.HOD/100,data.augSept.(variables{jj}),[],@nanmedian,NaN);
    tmp1std(:,jj)= accumarray(data.juneJuly.HOD/100,data.juneJuly.(variables{jj}),[],@nanstd,NaN);
    tmp2std(:,jj)= accumarray(data.augSept.HOD/100,data.augSept.(variables{jj}),[],@nanstd,NaN);
end

countJuneJulyIso=sum(~isnan(data.juneJuly.wehrReco_IFP));
countAugSeptIso=sum(~isnan(data.augSept.wehrReco_IFP));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load Rick Wehr's data independently
wehrData=readtable('../data_RickWehr/HF-PartitionedFluxes2.csv','Delimiter',',');
% process Wehr data to hourly
for jj=1:height(wehrData)
    wehrTime(jj,:) = strsplit(wehrData.DD_MM_YYYY_HH_MM_SS_EST_{jj},{'/','_',':'});
    
end % format of wehrTime dd mm yy hr:mm:ss

wehr.Reco_iso=wehrData.Reco_IFP_umol_m2_s_;
wehr.Reco_iso(wehr.Reco_iso==9999)=NaN;

% extract all non- june/july data OR august/september
monthz=cell2table(wehrTime)        ;
indX=strcmp(monthz.wehrTime2,'06') | strcmp(monthz.wehrTime2,'07');
wehr.JuneJuly=wehr;
wehr.JuneJuly.Reco_iso(~indX,:)=NaN;
indX=strcmp(monthz.wehrTime2,'08') | strcmp(monthz.wehrTime2,'09');
wehr.AugSept=wehr;
wehr.AugSept.Reco_iso(~indX,:)=NaN;

% delete all data not from SW direction (WD <270, WD> 180)
wehr.JuneJuly.Reco_iso_SW=wehr.JuneJuly.Reco_iso;
wehr.AugSept.Reco_iso_SW=wehr.AugSept.Reco_iso;
indX=wehrData.WindDirection_degFromNorth_>270 | wehrData.WindDirection_degFromNorth_ <180;
wehr.JuneJuly.Reco_iso_SW(indX,:)=NaN;
wehr.AugSept.Reco_iso_SW(indX,:)=NaN;

wehr.JuneJuly.Reco_iso_SE=wehr.JuneJuly.Reco_iso;
wehr.AugSept.Reco_iso_SE=wehr.AugSept.Reco_iso;
indX=wehrData.WindDirection_degFromNorth_>=90 & wehrData.WindDirection_degFromNorth_ < 180;
wehr.JuneJuly.Reco_iso_SE(indX,:)=NaN;
wehr.AugSept.Reco_iso_SE(indX,:)=NaN;

wehr.JuneJuly.Reco_iso_NE=wehr.JuneJuly.Reco_iso;
wehr.AugSept.Reco_iso_NE=wehr.AugSept.Reco_iso;
indX=wehrData.WindDirection_degFromNorth_<=90;
wehr.JuneJuly.Reco_iso_NE(indX,:)=NaN;
wehr.AugSept.Reco_iso_NE(indX,:)=NaN;

wehr.JuneJuly.Reco_iso_NW=wehr.JuneJuly.Reco_iso;
wehr.AugSept.Reco_iso_NW=wehr.AugSept.Reco_iso;
indX=wehrData.WindDirection_degFromNorth_>=270;
wehr.JuneJuly.Reco_iso_NW(indX,:)=NaN;
wehr.AugSept.Reco_iso_NW(indX,:)=NaN;


%%
% calculate wehrReco in two hour block means, to be consistent with
% Wehr et al. 2016
% wehrReco_IFP = 27
% wehrReco_IFP_sw = 35
variables{35}='wehrReco_IFP_SW';
variables{36}='wehrReco_IFP_SE';
variables{37}='wehrReco_IFP_NW';
variables{38}='wehrReco_IFP_NE';

hourz=str2num(cell2mat(monthz.wehrTime4));
X = hourz;
[B,idx] = histc(X,0:2:24);

% For June/July
V = accumarray(idx(:),wehr.JuneJuly.Reco_iso,[],@nanmean);
Vstd = accumarray(idx(:),wehr.JuneJuly.Reco_iso,[],@nanstd);
V_SW = accumarray(idx(:),wehr.JuneJuly.Reco_iso_SW,[],@nanmean);
V_SWstd = accumarray(idx(:),wehr.JuneJuly.Reco_iso_SW,[],@nanstd);
V_SE = accumarray(idx(:),wehr.JuneJuly.Reco_iso_SE,[],@nanmean);
V_SEstd = accumarray(idx(:),wehr.JuneJuly.Reco_iso_SE,[],@nanstd);
V_NE = accumarray(idx(:),wehr.JuneJuly.Reco_iso_NE,[],@nanmean);
V_NEstd = accumarray(idx(:),wehr.JuneJuly.Reco_iso_NE,[],@nanstd);
V_NW = accumarray(idx(:),wehr.JuneJuly.Reco_iso_NW,[],@nanmean);
V_NWstd = accumarray(idx(:),wehr.JuneJuly.Reco_iso_NW,[],@nanstd);

output = nan(1,2*numel(V));
output(2:2:end) = V;
output(1:2:end-1) = V;
tmp1(:,27)=output;
output = nan(1,2*numel(Vstd));
output(2:2:end) = Vstd;
output(1:2:end-1) = Vstd;
tmp1std(:,27)=output;

output = nan(1,2*numel(V_SW));
output(2:2:end) = V_SW;
output(1:2:end-1) = V_SW;
tmp1(:,35)=output;
output = nan(1,2*numel(V_SWstd));
output(2:2:end) = V_SWstd;
output(1:2:end-1) = V_SWstd;
tmp1std(:,35)=output;

output = nan(1,2*numel(V_SE));
output(2:2:end) = V_SE;
output(1:2:end-1) = V_SE;
tmp1(:,36)=output;
output = nan(1,2*numel(V_SEstd));
output(2:2:end) = V_SEstd;
output(1:2:end-1) = V_SEstd;
tmp1std(:,36)=output;

output = nan(1,2*numel(V_NW));
output(2:2:end) = V_NW;
output(1:2:end-1) = V_NW;
tmp1(:,37)=output;
output = nan(1,2*numel(V_SWstd));
output(2:2:end) = V_SWstd;
output(1:2:end-1) = V_SWstd;
tmp1std(:,37)=output;

output = nan(1,2*numel(V_NE));
output(2:2:end) = V_NE;
output(1:2:end-1) = V_NE;
tmp1(:,38)=output;
output = nan(1,2*numel(V_SWstd));
output(2:2:end) = V_SWstd;
output(1:2:end-1) = V_SWstd;
tmp1std(:,38)=output;

% For Aug/Sept
V = accumarray(idx(:),wehr.AugSept.Reco_iso,[],@nanmean);
Vstd = accumarray(idx(:),wehr.AugSept.Reco_iso,[],@nanstd);
V_SW = accumarray(idx(:),wehr.AugSept.Reco_iso_SW,[],@nanmean);
V_SWstd = accumarray(idx(:),wehr.AugSept.Reco_iso_SW,[],@nanstd);
V_SE = accumarray(idx(:),wehr.AugSept.Reco_iso_SE,[],@nanmean);
V_NE = accumarray(idx(:),wehr.AugSept.Reco_iso_NE,[],@nanmean);
V_NW = accumarray(idx(:),wehr.AugSept.Reco_iso_NW,[],@nanmean);

output = nan(1,2*numel(V));
output(2:2:end) = V;
output(1:2:end-1) = V;
tmp2(:,27)=output;
output = nan(1,2*numel(Vstd));
output(2:2:end) = Vstd;
output(1:2:end-1) = Vstd;
tmp2std(:,27)=output;

output = nan(1,2*numel(V_SW));
output(2:2:end) = V_SW;
output(1:2:end-1) = V_SW;
tmp2(:,35)=output;
output = nan(1,2*numel(V_SWstd));
output(2:2:end) = V_SWstd;
output(1:2:end-1) = V_SWstd;
tmp2std(:,35)=output;

output = nan(1,2*numel(V_SE));
output(2:2:end) = V_SE;
output(1:2:end-1) = V_SE;
tmp2(:,36)=output;
output = nan(1,2*numel(V_SEstd));
output(2:2:end) = V_SEstd;
output(1:2:end-1) = V_SEstd;
tmp2std(:,36)=output;

output = nan(1,2*numel(V_NW));
output(2:2:end) = V_NW;
output(1:2:end-1) = V_NW;
tmp2(:,37)=output;
output = nan(1,2*numel(V_SWstd));
output(2:2:end) = V_SWstd;
output(1:2:end-1) = V_SWstd;
tmp2std(:,37)=output;

output = nan(1,2*numel(V_NE));
output(2:2:end) = V_NE;
output(1:2:end-1) = V_NE;
tmp2(:,38)=output;
output = nan(1,2*numel(V_SWstd));
output(2:2:end) = V_SWstd;
output(1:2:end-1) = V_SWstd;
tmp2std(:,38)=output;



data.diurnal.juneJuly= array2table(tmp1,'VariableNames',variables);
data.diurnal.augSept= array2table(tmp2,'VariableNames',variables);
data.diurnal.juneJulyStd= array2table(tmp1std,'VariableNames',variables);
data.diurnal.augSeptStd= array2table(tmp2std,'VariableNames',variables);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   DIURNAL PLOTS OF RECO w/wout inhibition

for ii=1:2
    switch ii
        case 1 % June/July
            plotData=data.diurnal.juneJuly;
            plotDataStd=data.diurnal.juneJulyStd;
            countz= countJuneJulyIso;
            period='JuneJuly';
            baseline=5.1;
        case 2 % Aug/Sept
            plotData=data.diurnal.augSept;
            plotDataStd=data.diurnal.augSeptStd;
            countz= countAugSeptIso;
            period='AugSept';
            baseline=4.2;
    end
    LineW=3;
    
    fig1=figure;
    hold on
    
    dayBreak=find(plotData.SWin>0);
    
    % plot the Reco data
    clear p q
    
    cInd=3; % the SW quadrant presented in Wehr et al.
    ydata=plotData.wehrReco_IFP_SW(1:24);
    ydataStd=plotDataStd.wehrReco_IFP_SW(1:24);
    ydata(1:dayBreak(1))=baseline;
    ydata(end-3:end-1)=baseline;
    ydataStd(1:2)=nanmean(ydataStd);
    ydataStd(end-3:end)=nanmean(ydataStd);
    q(cInd)=shadedErrorBar(1:2:24,ydata(1:2:24),ydataStd(1:2:24)/((countz/12)^0.5),'--');
    set(q(cInd).patch,'facecolor',color.lightGray);
    set(q(cInd).mainLine,'color',color.darkGray);
    set(q(cInd).edge,'color','none');
    
    cInd=4; % all quadrants from Wehr et al.
    ydata=nanmean([plotData.wehrReco_IFP_SW,plotData.wehrReco_IFP_SE,...
        plotData.wehrReco_IFP_NW,plotData.wehrReco_IFP_NE]');
    ydata(1:dayBreak(1))=baseline;
    ydata(end-2:end)=baseline;
    q(cInd)=shadedErrorBar(1:2:24,ydata(1:2:24),ydataStd(1:2:24)/((countz/12)^0.5),'--',1);
    set(q(cInd).patch,'facecolor',color.saddleBrown);
    set(q(cInd).mainLine,'color',color.darkGray);
    set(q(cInd).edge,'color','none');
    
    
    cInd=1;
    q(cInd)=shadedErrorBar(1:24,plotData.RecoREddyProcNT,plotDataStd.RecoREddyProcNT/((60)^0.5),'--',0);
    set(q(cInd).patch,'facecolor',color.darkGreen);
    set(q(cInd).mainLine,'color',color.darkGray);
    set(q(cInd).edge,'color','none');
    
    cInd=2;
    q(cInd)=shadedErrorBar(1:24,plotData.RecoREddyProcNT_wInhib,plotDataStd.RecoREddyProcNT_wInhib/((60)^0.5),'--');
    set(q(cInd).patch,'facecolor',color.darkBlue);
    set(q(cInd).mainLine,'color',color.darkGray);
    set(q(cInd).edge,'color','none');
    
    if ii==1
        leg1=legend([q(1).patch,q(2).patch,q(3).patch,q(4).patch],...
            {'Standard, night-time, no inhibition','Standard, night-time, with inhibition',...
            'Isotopic (SW quadrants)','Isotopic (all quadrants)'});
        set(leg1,'FontSize',18,'box','off','location','NorthWest')
        l1=line([dayBreak(1) dayBreak(1)],[-50 6.8]);
    else
        l1=line([dayBreak(1) dayBreak(1)],[-50 10]);
        
    end
    l2=line([dayBreak(end) dayBreak(end)],[-50 10]);
    l3=line([0 25],[0 0]);
    set(l1,'color','k','LineStyle','--')
    set(l2,'color','k','LineStyle','--')
    set(l3,'color',[0.8 0.8 0.8],'LineStyle','-')
    
    xlim([1 24])
    ylim([0 9])
    
    xlabel('Hour of Day')
    ylabel('R_{eco} (\mumol m^{-2} s^{-1})')
    
    set(gca,'box','off')
    set(gca,'FontSize',22)
    set(gca,'XTick',[6,12,18])
    set(gca,'XTickLabel',{'6:00 am','12:00','6:00 pm'})
    
    if saveFigures==1
        filename=strcat('./figures/',period,'_',cSiteShort,'_RECO');
        saveas(fig1,filename,'png');
    end
    
    
end

%%   DIURNAL PLOTS OF RECO w/wout inhibition

for ii=1:2
    switch ii
        case 1 % June/July
            plotData=data.diurnal.juneJuly;
            plotDataStd=data.diurnal.juneJulyStd;
            countz= countJuneJulyIso;
            period='JuneJuly';
            baseline=5.1;
        case 2 % Aug/Sept
            plotData=data.diurnal.augSept;
            plotDataStd=data.diurnal.augSeptStd;
            countz= countAugSeptIso;
            period='AugSept';
            baseline=4.2;
    end
    LineW=3;
    
    fig1=figure;
    hold on
    
    dayBreak=find(plotData.SWin>0);
    
    % plot the Reco data
    clear p q
    
    cInd=3; % the SW quadrant presented in Wehr et al.
    ydata=plotData.wehrReco_IFP_SW(1:24);
    ydataStd=plotDataStd.wehrReco_IFP_SW(1:24);
    ydata(1:dayBreak(1))=baseline;
    ydata(end-3:end-1)=baseline;
    ydataStd(1:2)=nanmean(ydataStd);
    ydataStd(end-3:end)=nanmean(ydataStd);
    q(cInd)=shadedErrorBar(1:2:24,ydata(1:2:24),ydataStd(1:2:24)/((countz/12)^0.5),'--');
    set(q(cInd).patch,'facecolor',color.lightGray);
    set(q(cInd).mainLine,'color',color.darkGray);
    set(q(cInd).edge,'color','none');
    
    cInd=4; % all quadrants from Wehr et al.
    ydata=nanmean([plotData.wehrReco_IFP_SW,plotData.wehrReco_IFP_SE,...
        plotData.wehrReco_IFP_NW,plotData.wehrReco_IFP_NE]');
    ydata(1:dayBreak(1))=baseline;
    ydata(end-2:end)=baseline;
    q(cInd)=shadedErrorBar(1:2:24,ydata(1:2:24),ydataStd(1:2:24)/((countz/12)^0.5),'--',1);
    set(q(cInd).patch,'facecolor',color.saddleBrown);
    set(q(cInd).mainLine,'color',color.darkGray);
    set(q(cInd).edge,'color','none');
    
    
    cInd=1;
    q(cInd)=shadedErrorBar(1:24,plotData.RecoREddyProcNT,plotDataStd.RecoREddyProcNT/((60)^0.5),'--',0);
    set(q(cInd).patch,'facecolor',color.darkGreen);
    set(q(cInd).mainLine,'color',color.darkGray);
    set(q(cInd).edge,'color','none');
    
    cInd=2;
    q(cInd)=shadedErrorBar(1:24,plotData.RecoREddyProcNT_wInhib,plotDataStd.RecoREddyProcNT_wInhib/((60)^0.5),'--');
    set(q(cInd).patch,'facecolor',color.darkBlue);
    set(q(cInd).mainLine,'color',color.darkGray);
    set(q(cInd).edge,'color','none');
    
    if ii==1
        leg1=legend([q(1).patch,q(2).patch,q(3).patch,q(4).patch],...
            {'Standard, night-time, no inhibition','Standard, night-time, with inhibition',...
            'Isotopic (SW quadrant)','Isotopic (all quadrants)'});
        set(leg1,'FontSize',18,'box','off','location','NorthWest')
        l1=line([dayBreak(1) dayBreak(1)],[-50 6.8]);
    else
        l1=line([dayBreak(1) dayBreak(1)],[-50 10]);
        
    end
    l2=line([dayBreak(end) dayBreak(end)],[-50 10]);
    l3=line([0 25],[0 0]);
    set(l1,'color','k','LineStyle','--')
    set(l2,'color','k','LineStyle','--')
    set(l3,'color',[0.8 0.8 0.8],'LineStyle','-')
    
    xlim([1 24])
    ylim([0 9])
    
    xlabel('Hour of Day')
    ylabel('R_{eco} (\mumol m^{-2} s^{-1})')
    
    set(gca,'box','off')
    set(gca,'FontSize',22)
    set(gca,'XTick',[6,12,18])
    set(gca,'XTickLabel',{'6:00 am','12:00','6:00 pm'})
    
    if saveFigures==1
        filename=strcat('./figures/',period,'_',cSiteShort,'_Reco');
        saveas(fig1,filename,'png');
    end
    
    
end


