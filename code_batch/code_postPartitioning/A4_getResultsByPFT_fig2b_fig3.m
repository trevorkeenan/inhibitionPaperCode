% this script will...

% load the list of FLUXNET sites

% identify the unique PFTs

% for each PFT,
%       get the list of sites
%       load the site data
%       and calculate PFT mean responses and statistics

% T. Keenan, November 2018

close all
clear all

saveFigures=0;
generateSitePlots=0;    % this will generate the seasonal plot for each site
% and save them in PFT specific folders

airTsoilT='airT'; % options, 'airT','soilT' define the temperature used for partitioning

addpath('./functions')

% load the list of Fluxnet sites
sites=readtable('../data_FLUXNET2015release3/siteinfo_fluxnet2015_Aug17original.csv','Delimiter',',');

siteInfo = sites;

% get unique PFTs
uniquePFTs=unique(sites.classid);

% define some colores
color.gray=[0.2 0.2 0.2];
color.darkGray=[0.2 0.2 0.2];
color.lightGray=[0.7 0.7 0.7];
color.barkBrown=[0.396 0.263 0.129]; % 39.6% red, 26.3% green and 12.9% blue
color.saddleBrown=[139 69 19]/255;
color.blue=[0 0 1];
color.darkBlue=[0 0 0.5];
color.deepBlue=[0 0 0.2];
color.green=[0 1 0];
color.darkGreen=[0 0.5 0];
color.deepGreen=[0 0.2 0];

% loop through each PFT and get site
for ii=1:length(uniquePFTs)
    cPFT=uniquePFTs{ii};
    
    % find sites of this PFT
    cPFTsiteNames=sites.mysitename(strcmp(sites.classid,cPFT));
    
    % loop through sites, load data and aggregate
    for jj=1:length(cPFTsiteNames)
        close all
        cSite=cPFTsiteNames{jj};
        
        filename=strcat('./data_meanMonthly/Tier1only_gs95/',cSite,'_',airTsoilT,'.mat');
        
        try % there is no data if the site is Tier 2
            cSiteData=load(filename);
            
            pftData{ii}(jj,:,:)=table2array(cSiteData.T);
        catch
            pftData{ii}(jj,:,:)=table2array(cSiteData.T)*NaN;
            disp(strcat('skipping Tier 2 ', ' ', cSite))
        end
        
        % adjust for southern hemisphere season
        if strcmp(cSite(1:2),'AU') || strcmp(cSite(1:2),'AR') || strcmp(cSite(1:2),'ZA') || strcmp(cSite(1:2),'ZM')
            cSiteData.Tmp1=cSiteData.T(1:6,:);
            cSiteData.Tmp2=cSiteData.T(7:end,:);
            cSiteData.T2= vertcat(cSiteData.Tmp2,cSiteData.Tmp1);
            cSiteData.T2.Month=[1:12]';
            cSiteData.T=cSiteData.T2;
        end
        
        % select only sites with 5+ years
        if mean(cSiteData.T.numYears)>=5
            pftDataSelect{ii}(jj,:,:)=table2array(cSiteData.T);
        else
            pftDataSelect{ii}(jj,:,:)=NaN*table2array(cSiteData.T);
            pftData{ii}(jj,:,:)=table2array(cSiteData.T);
        end
    end
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Calculate and plot the mean inhibition per PFT

% get the mean inhibition per PFT
for ii=1:length(uniquePFTs)
    
    cPFTdata=pftData{ii};
    cPFTdataSelect=pftDataSelect{ii};
    
    % get the perc inhibition data for this PFT
    cPercInhib=squeeze(cPFTdata(:,:,6));
    cPercInhibSelect=squeeze(cPFTdataSelect(:,:,6));
    
    cMeanInh(ii,1)=nanmedian(cPercInhib(:));
    cMeanInh(ii,2)=nanstd(cPercInhib(:));
    cMeanInh(ii,3)=size(cPFTdata,1);
    
    cMeanInhSel(ii,1)=nanmedian(cPercInhibSelect(:));
    cMeanInhSel(ii,2)=nanstd(cPercInhibSelect(:));
    cMeanInhSel(ii,3)=size(cPFTdataSelect,1);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate the seasonal cycle in inhibition

for ii=1:length(uniquePFTs)
    
    cPFTdataSelect=pftDataSelect{ii};
    
    % count all sites that are not just NaNs
    countx=0;
    for jj=1:size(cPFTdataSelect)
        tmp=squeeze(cPFTdataSelect(jj,:,:));
        if nansum(tmp(:)) ~= 0
            countx=countx+1;
        end
    end
    % get the perc inhibition data for this PFT
    cPercInhibSelect=squeeze(cPFTdataSelect(:,:,6));
    
    seasonalInhSelMean(ii,:)=nanmedian(cPercInhibSelect,1);
    seasonalInhSelStd(ii,:)=nanstd(cPercInhibSelect,1);
    seasonalInhSelNsites(ii,:)=countx;
    
    
end
%%
%   Plot the percent inhibition by PFT including only sites with 5+ years

fig1=figure;
hold on
y = cMeanInhSel(:,1);
[~,I] = sort(y,'descend');
errY = sqrt(cMeanInhSel(:,3));          

ySort=y(I);
errSort=errY(I);
pftSort=uniquePFTs(I);
indX=isnan(ySort);
ySort(indX)=[];
errSort(indX)=[];
pftSort(indX)=[];

nSites=seasonalInhSelNsites;
nSitesSort=nSites(I);
nSitesSort(indX)=[];

disp('data from the barchart')
disp(ySort)
disp('+-')
disp(errSort)

% plot
for ii = 1:length(ySort)
    h(ii)=bar(ii,ySort(ii));
    
    set(h(ii),'FaceColor',[0.5 0.5 0.5]);
end
ebar=errorbar(1:length(ySort),ySort,errSort,'k.','LineWidth',2);


fontSizeX=14;
% annotate the statistics
currentPos=get(gca,'Position');
currentPos2=currentPos;
currentPos2(1)=currentPos2(1);  % left/right
currentPos2(2)=currentPos2(2);
%     currentPos2(3)=currentPos2(3)-0.3;
currentPos2(4)=currentPos2(4)-0.675; % up down

annotation(fig1,'textbox',currentPos2,...
    'String','n = ',...
    'FontSize',fontSizeX,'FitBoxToText','on', 'LineStyle','none');

currentPos2(4)=currentPos2(4)-0.05; % up down
currentPos2(1)=currentPos2(1)+0.05; % up down

for ijk = 1:11
    annotation(fig1,'textbox',currentPos2,...
        'String',num2str(nSitesSort(ijk)),...
        'FontSize',fontSizeX,'FitBoxToText','on', 'LineStyle','none');
    currentPos2(1)=currentPos2(1)+0.05;  % left/right
    if ijk ==2
        currentPos2(1)=currentPos2(1)+0.01;  % left/right
    end
    if ijk ==3
        currentPos2(1)=currentPos2(1)+0.01;  % left/right
    end
    if ijk ==4
        currentPos2(1)=currentPos2(1);%+0.015;  % left/right
    end
    if ijk ==5
        currentPos2(1)=currentPos2(1)+0.01;  % left/right
    end
    if ijk ==6
        currentPos2(1)=currentPos2(1)+0.01;  % left/right
    end
    if ijk ==7
        currentPos2(1)=currentPos2(1)+0.01;  % left/right
    end
    if ijk ==8
        currentPos2(1)=currentPos2(1)+0.01;  % left/right
    end
    if ijk ==9
        currentPos2(1)=currentPos2(1)+0.01;  % left/right
    end
    if ijk ==10
        currentPos2(1)=currentPos2(1)+0.005;  % left/right
    end
    
end


ylabel('Inhibition (%)')
set(gca,'box','off','fontsize',18)
xlim([0 length(pftSort)+1])
ylim([-5 30])

set(gca,'XTickLabel',pftSort,'XTick',1:length(pftSort))
xticklabel_rotate([],45,[],'Fontsize',14)

if saveFigures==1
    filename=strcat('./figures/inhibitionByPFTselect/inhibitionByPFT','_',airTsoilT,'_tier1only_gs95');
    saveas(fig1,filename,'png');
    
end


%%
fig1=figure;
hold on

cPFTind=5;
p(3)= shadedErrorBar(1:12,smooth(seasonalInhSelMean(cPFTind,:)),seasonalInhSelStd(cPFTind,:)/sqrt(seasonalInhSelNsites(cPFTind)),'--',0);
set(p(3).patch,'facecolor',color.lightGray);
set(p(3).mainLine,'color',color.darkGray);
set(p(3).edge,'color','none');

cPFTind=3;
p(1)= shadedErrorBar(1:12,(seasonalInhSelMean(cPFTind,:)),seasonalInhSelStd(cPFTind,:)/sqrt(seasonalInhSelNsites(cPFTind)),'--',0);
set(p(1).patch,'facecolor',color.darkGreen);
set(p(1).mainLine,'color',color.deepGreen);
set(p(1).edge,'color','none');

cPFTind=6;
p(2)= shadedErrorBar(1:12,smooth(seasonalInhSelMean(cPFTind,:),3),seasonalInhSelStd(cPFTind,:)/sqrt(seasonalInhSelNsites(cPFTind)),'--',0);
set(p(2).patch,'facecolor',color.darkBlue);
set(p(2).mainLine,'color',color.deepBlue);
set(p(2).edge,'color','none');
set(p(2).patch,'FaceAlpha',0.5);


l1=legend([p(1).patch,p(2).patch,p(3).patch],{'Deciduous Broadleaf Forests','Evergreen Needleleaf Forests','Evergreen Broadleaf Forests'});
set(l1,'box','off','location','NorthEast')
loc1=get(l1,'position'); % [left bottom width height]
% nudge it right a bit
loc2=loc1;
loc2(1)=loc2(1)+0.00005;
set(l1,'position',loc2)


ylabel('Inhibition (%)')
xlabel('Month')
set(gca,'XTick',1:12)
set(gca,'box','off','fontsize',18)

ylim([0 30])
xlim([0.5 12.5])


if saveFigures==1
    filename=strcat('./figures/inhibitionByPFTselect/DBFvsENFvsEBF','_',airTsoilT,'_tier1only_gs95');
    saveas(fig1,filename,'png');
    
end


