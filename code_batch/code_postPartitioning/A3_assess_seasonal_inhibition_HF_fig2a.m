% this script will...

% load the parameters estimated by REddyProc
% compare RrefDay with RrefNight
% and extract monthly values for plotting later

% now using only growing season data for calculating % inhibition

% the results here are not dependent on whether the input is with or
% without inhibition, as it just compares R_night to R_day

% T. Keenan, November 2018

clearvars

saveFigures=1;
plotSiteFigures=1;
noFAPAR=1;      % set to 1 if not using FAPAR

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

% run just on Harvard Forest (US-Ha1)
ii =156 % harvard forest is site # 156, though this depends on # of sites processed
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
    
    
    cSite=sites(ii);
    cSiteShort=cSite{1}(1:6);
    disp(cSiteShort)
    
    filename=strcat(home, cSiteShort,'_devFAPAR');
    hourlyData2=load(filename);
    hourlyData=hourlyData2.T;
    
    fAPAR=hourlyData.fAPAR;
    
    % remove the first year fAPAR
    indX=hourlyData.Year==2000;
    hourlyData.fAPAR(indX)=NaN;
    
end

%%
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
    y=param.(current);
    y(GPP95<1)=NaN;
    param.(current)=y;
    
    y=param.(current);
    y(y==-9999)=NaN;
    param.(current)=y;
end


% look at R_ref to R_night ratio for summer months
%%
Rref=param.R_ref;
RrefNight=param.R_night;
qc=param.qc;

%%
Rday_mean=nanmean(Rref);
Rday_ste=nanstd(Rref)/(sum(~isnan(Rref))^0.5);
Rnight_mean=nanmean(RrefNight);
Rnight_ste=nanstd(RrefNight)/(sum(~isnan(Rref))^0.5);

% get seasonal patterns of Rref and Rref_Night

%% get mean monthly values for each

paramNo20112012 = param;
if strcmp(cSiteShort,'US-Ha1')
    indX=param.Year==1991 | param.Year==1992; % '91, '92 are bad data for HF. No data for 91 and 92 SW_IN is corrupt
    paramNo20112012(indX,:)=[];
    GPP95(indX)=[];
    indX=param.Year==2010 | param.Year==2011; % '10,'11 are odd years.
    paramNo20112012(indX,:)=[];
    GPP95(indX)=[];
end

% clearn the param file by removing large outliers
paramNo20112012.R_night(paramNo20112012.R_night>40)=NaN;
paramNo20112012.R_ref(paramNo20112012.R_ref>40)=NaN;
paramNo20112012.R_night(paramNo20112012.R_night<=0)=NaN;
paramNo20112012.R_ref(paramNo20112012.R_ref<=0)=NaN;

paramNo20112012.percentInhibition=100*(paramNo20112012.R_night-paramNo20112012.R_ref)./paramNo20112012.R_night;

% clear the percent inhibition by removing outliers
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
monthlyParam.Rref=accumarray(paramNo20112012.Month,paramNo20112012.R_ref,[],@nanmedian,NaN);
monthlyParam.RrefStd= accumarray(paramNo20112012.Month,paramNo20112012.R_ref,[],@nanstd,NaN);
monthlyParam.RrefNight= accumarray(paramNo20112012.Month,paramNo20112012.R_night,[],@nanmedian,NaN);
monthlyParam.RrefNightStd= accumarray(paramNo20112012.Month,paramNo20112012.R_night,[],@nanstd,NaN);
monthlyParam.percentInh= accumarray(paramNo20112012.Month,paramNo20112012.percentInhibition,[],@nanmedian,NaN);
monthlyParam.percentInhStd= accumarray(paramNo20112012.Month,paramNo20112012.percentInhibition,[],@nanstd,NaN);

% get the monthly 95th percentile GPP
monthlyParam.GPP95= accumarray(paramNo20112012.Month,GPP95,[],@nanmean,NaN);

monthlyParam.GPP95_growingSeason=monthlyParam.GPP95;
monthlyParam.GPP95_growingSeason(monthlyParam.GPP95_growingSeason<1)=NaN;

if noFAPAR==0
    monthlyParam.fAPARmean= accumarray(hourlyData.Month,hourlyData.fAPAR,[],@nanmedian,NaN);
    monthlyParam.fAPARstd= accumarray(hourlyData.Month,hourlyData.fAPAR,[],@nanstd,NaN);
end

% load the fAPAR data for Harvard Forest
homeFAPARdata='../modis-land-product-subset-master/Rcode/formattedData/US-Ha1_MOD15A2_Fpar_1km.csv';

datainFAPAR = csvread(homeFAPARdata,9,2);

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
areaMedian=nanmedian(datainFAPAR(:,3:end),2);
% remove 2000, 2016 to just use full years
indX=datainFAPAR(:,1)==2000 | datainFAPAR(:,1)==2016;
areaMedian(indX)=[];
datainFAPAR(indX,:)=[];

monthlyParam.fAPARdoy= accumarray(datainFAPAR(:,2),datainFAPAR(:,2),[],@nanmedian,NaN);
monthlyParam.fAPARmean= accumarray(datainFAPAR(:,2),areaMedian,[],@nanmedian,NaN);
monthlyParam.fAPARstd= accumarray(datainFAPAR(:,2),areaMedian,[],@nanstd,NaN);

numYears=length(unique(paramNo20112012.Year));

% define some colors
color.gray=[0.2 0.2 0.2];
color.darkGray=[0.2 0.2 0.2];
color.lightGray=[0.7 0.7 0.7];
color.barkBrown=[0.396 0.263 0.129]; % 39.6% red, 26.3% green and 12.9% blue
color.blue=[0 0 1];
color.green=[0 1 0];


%% Plot the seasonal cycle of Rref Day and Night
%       and the difference between them (Inhibition)
if plotSiteFigures==1
    
    scrsz = get(0,'ScreenSize');
    fig1 =figure('Position',[1 scrsz(4) scrsz(3)/2.5 scrsz(4)/2.5]);
    
    
    % set up the top panel axis
    h1 = axes('Parent',fig1);
    set(gca,'YAxisLocation','left');
    p = get(h1, 'pos');  %This is a 4-element vector [left, bottom, width, height]
    px(1) = p(1)+0.08;
    px(2) = p(2)+(p(4)*(1-0.75));
    px(3) = p(3)*0.8;
    px(4) = p(4)*0.75;
    set(h1, 'pos', px);
    hold on
    set(gca,'XTick',[])
    
    title(strcat(cSiteShort,{': '},'#years = ',num2str(length(unique(param.Year)))))
    
    
    % plot a shaded background mean seasonal cycle
    p=shadedErrorBar(1:12,monthlyParam.RrefNight,monthlyParam.RrefNightStd/(numYears^0.5),'--',0.1);
    set(p.patch,'facecolor',color.gray);
    set(p.mainLine,'color',color.gray);
    set(p.edge,'color','none');
    p=shadedErrorBar(1:12,monthlyParam.Rref,monthlyParam.RrefStd/(numYears^0.5),'--',0.1);
    set(p.patch,'facecolor',color.gray);
    set(p.mainLine,'color',color.gray);
    set(p.edge,'color','none');
    
    % plot the growing season mean seasonal cycle
    % identify times when GPP is not zero
    indX = ~isnan(monthlyParam.GPP95_growingSeason);
    monthz=1:12;
    p1=shadedErrorBar(monthz(indX),monthlyParam.RrefNight(indX),monthlyParam.RrefNightStd(indX)/(numYears^0.5),'b--');
    p2=shadedErrorBar(monthz(indX),monthlyParam.Rref(indX),monthlyParam.RrefStd(indX)/(numYears^0.5),'g--');
    
    set(p1.edge,'color','none');
    set(p2.edge,'color','none');
    
    % add a legend
    maxYtop=max(max(monthlyParam.Rref(3:11),monthlyParam.RrefNight(3:11)))+1;
    if ~isnan(maxYtop)
        l1=legend([p1.patch,p2.patch],{'Night','Day'});
        set(l1,'location','NorthWest','box','off')
        
        ylim([0 maxYtop])
        xlim([0.5 12.5])
    end
    ylabel({'R_{ref} (\mumol m^{-2} s^{-1})'})
    set(gca,'FontSize',18)
    
    % set up the second axis
    h2 = axes('Parent',fig1);
    set(gca,'YAxisLocation','right')
    px2 = get(h1, 'pos');  %This is a 4-element vector [left, bottom, width, height]
    px2(1) = px2(1);
    px2(2) = px2(2)-0.21;
    px2(3) = px2(3);
    px2(4) = 0.2;
    set(h2, 'pos', px2);
    hold on
    set(gca,'XTick',[],'FontSize',18)
    
    % plot fAPAR
    if noFAPAR==0
        pFAPAR=shadedErrorBar(1:12,monthlyParam.fAPARmean,monthlyParam.fAPARstd/(numYears^0.5),'--',0.2);
        set(pFAPAR.patch,'facecolor',color.barkBrown);
        set(pFAPAR.mainLine,'color',color.barkBrown);
        set(pFAPAR.edge,'color','none');
        xlim([0.5 12.5])
        ylim([0 100])
        
        ylabel('fAPAR (%)')
    end
    % set up the top panel fAPAR axis
    h0 = axes('Parent',fig1);
    set(gca,'color','none');
    px = get(h2, 'pos');  %This is a 4-element vector [left, bottom, width, height]
    set(h0, 'pos', px);
    hold on
    set(gca,'XTick',[])
    
    % plot the difference between the two Rref estimates
    p=shadedErrorBar(1:12,monthlyParam.percentInh,monthlyParam.percentInhStd/(numYears^0.5),'--');
    set(p.patch,'facecolor',color.lightGray);
    set(p.mainLine,'color','k');
    set(p.edge,'color','none');
    
    ylabel({'Inhibition (%)'})
    xlabel('Month')
    set(gca,'XTick',1:12)
    
    ylim([0 50])
    xlim([0.5 12.5])
    set(gca,'YTick',[0,15,30,45])
    set(gca,'FontSize',18)
    
    % add a legend
    if ~isnan(maxYtop)
        if noFAPAR==0
            l1=legend([p.patch,pFAPAR.patch],{'Inhibition','fAPAR'});
        else
            l1=legend([p.patch],{'Inhibition'});
        end
        set(l1,'location','NorthWest','box','off','FontSize',10)
    end
    
    
    if saveFigures==1
        filename=strcat('./figures/nightVsDay_All/',cSiteShort,'_RrefAndRnight_monthly_inhibition1','_',airTsoilT);
        saveas(fig1,filename,'png');
    end
    
    
    %% Plot the seasonal cycle of Rref Day and Night 2
    %       and the difference between them (Inhibition)
    
    scrsz = get(0,'ScreenSize');
    fig1 =figure; %('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
    
    
    % set up the top panel axis
    h1 = axes('Parent',fig1);
    set(gca,'YAxisLocation','left');
    p = get(h1, 'pos');  %This is a 4-element vector [left, bottom, width, height]
    px(1) = p(1)+0.08;
    px(2) = p(2)+(p(4)*(1-0.75));
    px(3) = p(3)*0.8;
    px(4) = p(4)*0.75;
    set(h1, 'pos', px);
    hold on
    set(gca,'XTick',[])
    
    % plot a shaded background mean seasonal cycle
    p=shadedErrorBar(1:12,monthlyParam.RrefNight2,monthlyParam.RrefNightStd2/(numYears^0.5),'--',0.1);
    set(p.patch,'facecolor',color.gray);
    set(p.mainLine,'color',color.gray);
    set(p.edge,'color','none');
    p=shadedErrorBar(1:12,monthlyParam.Rref2,monthlyParam.RrefStd2/(numYears^0.5),'--',0.1);
    set(p.patch,'facecolor',color.gray);
    set(p.mainLine,'color',color.gray);
    set(p.edge,'color','none');
    
    % plot the growing season mean seasonal cycle
    indX = ~isnan(monthlyParam.GPP95_growingSeason);
    
    p1=shadedErrorBar(3:11,monthlyParam.RrefNight2(3:11),monthlyParam.RrefNightStd2(3:11)/(numYears^0.5),'b--');
    p2=shadedErrorBar(3:11,monthlyParam.Rref2(3:11),monthlyParam.RrefStd2(3:11)/(numYears^0.5),'g--');
    set(p1.edge,'color','none');
    set(p2.edge,'color','none');
    
    % add a legend
    
    maxYtop=nanmax(nanmax(monthlyParam.Rref2(3:11),monthlyParam.RrefNight2(3:11)))+1;
    if ~isnan(maxYtop)
        l1=legend([p1.patch,p2.patch],{'\itR\rm_{ref}^{ N}','\itR\rm^{ D}_{ref}'});
        set(l1,'location','NorthWest','box','off','FontSize',16)
        
        currentPos=get(l1,'Position');
        currentPos2=currentPos;
        currentPos2(1)=currentPos2(1)-0.02;  % left/right
        set(l1,'Position',currentPos2);
        
        
        % plot vertical lines
        plot([3 3], [0 maxYtop],'--','color',color.lightGray);
        plot([11 11], [0 maxYtop],'--','color',color.lightGray);
        
        ylim([1.8 maxYtop])
    end
    ylabel({'R_{ref} (\mumol m^{-2} s^{-1})'})
    set(gca,'FontSize',18)
    xlim([0.5 12.5])
    
    % add text for growing season and dormant season to the graph
    string=strcat('Growing Season');
    % annotate the statistics
    currentPos=get(gca,'Position');
    currentPos2=currentPos;
    currentPos2(1)=currentPos2(1)+0.14;  % left/right
    currentPos2(2)=currentPos2(2);
    %     currentPos2(3)=currentPos2(3)-0.3;
    currentPos2(4)=currentPos2(4)-0.57; % up down
    
    annotation(fig1,'textbox',...
        currentPos2,...
        'String',string,...
        'FontSize',14,...
        'FitBoxToText','on', 'LineStyle','none');
    
    % add text for growing season and dormant season to the graph
    string=strcat({'Dormant';'Season'});
    % annotate the statistics
    currentPos=get(gca,'Position');
    currentPos2=currentPos;
    currentPos2(1)=currentPos2(1);  % left/right
    currentPos2(2)=currentPos2(2);
    %     currentPos2(3)=currentPos2(3)-0.3;
    currentPos2(4)=currentPos2(4)-0.54; % up down
    
    annotation(fig1,'textbox',...
        currentPos2,...
        'String',string,...
        'FontSize',14,...
        'FitBoxToText','on', 'LineStyle','none');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set up the second axis
    h2 = axes('Parent',fig1);
    set(gca,'YAxisLocation','right')
    px2 = get(h1, 'pos');  %This is a 4-element vector [left, bottom, width, height]
    px2(1) = px2(1);
    px2(2) = px2(2)-0.21;
    px2(3) = px2(3);
    px2(4) = 0.2;
    set(h2, 'pos', px2);
    hold on
    set(gca,'XTick',[],'FontSize',18)
    
    % plot fAPAR
    if noFAPAR==0
        pFAPAR=shadedErrorBar(1:12,monthlyParam.fAPARmean,monthlyParam.fAPARstd/(numYears^0.5),'--',0.2);
        set(pFAPAR.patch,'facecolor',color.barkBrown);
        set(pFAPAR.mainLine,'color',color.barkBrown);
        set(pFAPAR.edge,'color','none');
        xlim([0.5 12.5])
        ylim([0 100])
        
        ylabel('fAPAR (%)')
    end
    
    xdata=12*(monthlyParam.fAPARdoy/361);
    xdata(isnan(xdata))=[];
    ydata=monthlyParam.fAPARmean;
    ydata(isnan(ydata))=[];
    zdata=monthlyParam.fAPARstd;
    zdata(isnan(zdata))=[];
    
    pFAPAR=shadedErrorBar(xdata,ydata,zdata/(15^0.5),'--',0.2);
    set(pFAPAR.patch,'facecolor',color.barkBrown);
    set(pFAPAR.mainLine,'color',color.barkBrown);
    set(pFAPAR.edge,'color','none');
    xlim([0.5 12.5])
    ylim([0 100])
    
    ylabel('fAPAR (%)')
    
    % set up the top panel fAPAR axis
    h0 = axes('Parent',fig1);
    set(gca,'color','none');
    px = get(h2, 'pos');  %This is a 4-element vector [left, bottom, width, height]
    set(h0, 'pos', px);
    hold on
    set(gca,'XTick',[])
    
    % plot the difference between the two Rref estimates
    p=shadedErrorBar(1:12,monthlyParam.percentInh2,monthlyParam.percentInhStd2/(numYears^0.5),'--');
    set(p.patch,'facecolor',color.lightGray);
    set(p.mainLine,'color','k');
    set(p.edge,'color','none');
    
    ylabel({'Inhibition (%)'})
    xlabel('Month')
    set(gca,'XTick',1:12)
    
    ylim([0 50])
    xlim([0.5 12.5])
    set(gca,'YTick',[0,15,30,45])
    set(gca,'FontSize',18)
    
    if ~isnan(maxYtop)
        
        l1=legend([p.patch,pFAPAR.patch],{'Inhibition','fAPAR'});
        set(l1,'location','NorthWest','box','off','FontSize',14)
        loc1=get(l1,'position'); % [left bottom width height]
        
        % nudge it up a bit
        loc2=loc1;
        loc2(2)=loc2(2)+0.025;
        set(l1,'position',loc2)
        
    end
    set(fig1,'PaperPositionMode', 'auto');
    
    if saveFigures==1
        filename=strcat('./figures/inhibitionHF/',cSiteShort,'_RrefAndRnight_monthly_inhibition3','_',airTsoilT);
        saveas(fig1,filename,'png');
    end
    
end

