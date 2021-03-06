% this script will...

% load the bias data generated by calculateAnnualTotals_DT and _NT

% and plot a bar chart of the biases

% T. Keenan, November 2018

saveFigures=1;
close all

biasDT=load('./data_inter/bias_DT.mat','bias');
biasDT=biasDT.bias;
biasNT=load('./data_inter/bias_NT.mat','bias');
biasNT=biasNT.bias;


addpath('./functions')

%%
%   Plot the percent bias by flux
fig1=figure;
hold on

% set the night time data
barData=100*[biasNT.GPPgs,biasNT.Reco,biasNT.RecoGs,biasNT.RecoGsDt];
errBarData=100*[biasNT.GPPgsdt,biasNT.Reco2_std,biasNT.Reco2gs_std,biasNT.Reco2gsdt_std]/sqrt(biasNT.count);
% now add the daytime method biases
barData2=100*abs([biasDT.GPP,biasDT.Reco,biasDT.Reco2gs_mean,biasDT.Reco2gsnt_mean]);
errBarData2=100*[biasDT.GPP2_std,biasDT.Reco2_std,biasDT.Reco2gs_std,biasDT.Reco2gsnt_std]/sqrt(biasDT.count);

numbars=1+length(barData)+length(barData2); % adding one for an empty bar

% plot
for ii = 1:length(barData)
    h(ii)=bar(ii,barData(ii));
    
    set(h(ii),'FaceColor',[0.3 0.3 0.3]);
end
ebar=errorbar(1:length(barData),barData,errBarData,'k.','LineWidth',2);


% plot
for ii = 1+(length(barData)+1):numbars
    h(ii)=bar(ii,barData2(ii-5));
    
    set(h(ii),'FaceColor',[0.8 0.8 0.8]);
end
ebar=errorbar((1+length(barData)+1):numbars,barData2,errBarData2,'k.','LineWidth',2);

xlabelz={'\it{F}\rm_p (Annual)','\it{F}\rm_r (Annual)','\it{F}\rm_r (GS)','\it{F}\rm_r (GS, DT)',...
            ' ',...
            '\it{F}\rm_p (Annual)','\it{F}\rm_r (Annual)','\it{F}\rm_r (GS)','\it{F}\rm_r (GS, NT)'};

% plot vertical lines
plot([5 5], [0 20],'--','color',[0.3 0.3 0.3]);

% add text for the Day-time and Night-time methods
string=strcat('Night-time method');
currentPos=get(gca,'Position');
currentPos2=currentPos;
currentPos2(1)=currentPos2(1);  % left/right
currentPos2(2)=currentPos2(2);    
%     currentPos2(3)=currentPos2(3)-0.3;
currentPos2(4)=currentPos2(4); % up down

annotation(fig1,'textbox',...
currentPos2,...
'String',string,...
'FontSize',14,...
'FitBoxToText','on', 'LineStyle','none');
    
% add text for the Day-time and Night-time methods
string=strcat('Day-time method');
currentPos=get(gca,'Position');
currentPos2=currentPos;
currentPos2(1)=currentPos2(1)+0.33;  % left/right
currentPos2(2)=currentPos2(2);    
%     currentPos2(3)=currentPos2(3)-0.3;
currentPos2(4)=currentPos2(4); % up down

annotation(fig1,'textbox',...
currentPos2,...
'String',string,...
'FontSize',14,...
'FitBoxToText','on', 'LineStyle','none');


ylabel('| Bias | (%)')

set(gca,'box','off','fontsize',18)
xlim([0 numbars+1])
set(gca,'XTickLabel',xlabelz,'XTick',1:numbars)
xticklabel_rotate([],45,[],'Fontsize',14)

if saveFigures==1
    filename=strcat('./figures/bias_barchart');
    saveas(fig1,filename,'png');
    
end

%%
%   Plot the percent bias by flux
fig1=figure;
hold on
% Create second Y axes on the right.
a2 = gca(fig1);
set(a2,'YAxisLocation', 'Right')
set(a2,'YTick',[0,5,10,15,20],...
    'YTickLabel', ...
    {'0','-5','-10','-15','-20'})
% Hide second plot.
set(a2, 'color', 'none')
set(a2, 'XTick', 1:9)
set(a2, 'YLim', [0 20])
plot(1:9,0,'.')

% Create second Y axes on the right.
a1 = axes('YAxisLocation', 'Left');
% Hide second plot.
set(a1, 'color', 'none')
set(a1, 'XTick', 1:9)

hold all

% set the night time data
barData=100*[biasNT.GPPgs,biasNT.Reco,biasNT.RecoGs,biasNT.RecoGsDt];
errBarData=100*[biasNT.GPPgsdt,biasNT.Reco2_std,biasNT.Reco2gs_std,biasNT.Reco2gsdt_std]/sqrt(biasNT.count);
disp('Night time Bias')
disp('GPPgs, Reco, Reco gs, Reco gs dt')
disp(barData)
disp('+/-')
disp(errBarData)

% now add the daytime method biases
barData2=100*abs([biasDT.GPP,biasDT.Reco,biasDT.Reco2gs_mean,biasDT.Reco2gsnt_mean]);
errBarData2=100*[biasDT.GPP2_std,biasDT.Reco2_std,biasDT.Reco2gs_std,biasDT.Reco2gsnt_std]/sqrt(biasDT.count);
disp('Day time Bias')
disp('GPP, Reco, Reco gs, Reco gs nt')
disp(barData2)
disp('+/-')
disp(errBarData2)

numbars=1+length(barData)+length(barData2); % adding one for an empty bar

% plot
for ii = 1:length(barData)
    h(ii)=bar(ii,barData(ii));
    
    set(h(ii),'FaceColor',[0.3 0.3 0.3]);
end
ebar=errorbar(1:length(barData),barData,errBarData,'k.','LineWidth',2);


% plot the day-time columns
for ii = 1+(length(barData)+1):numbars
    h(ii)=bar(ii,barData2(ii-5));
 
    set(h(ii),'FaceColor',[0.8 0.8 0.8]);
end
ebar=errorbar((1+length(barData)+1):numbars,barData2,errBarData2,'k.','LineWidth',2);

xlabelz={'\it{F}\rm_p (Annual)','\it{F}\rm_r (Annual)','\it{F}\rm_r (GS)','\it{F}\rm_r (GS, DT)',...
            ' ',...
            '\it{F}\rm_p (Annual)','\it{F}\rm_r (Annual)','\it{F}\rm_r (GS)','\it{F}\rm_r (GS, NT)',' '};

% plot vertical lines
plot([5 5], [0 20],'--','color',[0.3 0.3 0.3]);

% add text for the Day-time and Night-time methods
string=strcat('Night-time method');
currentPos=get(gca,'Position');
currentPos2=currentPos;
currentPos2(1)=currentPos2(1);  % left/right
currentPos2(2)=currentPos2(2);    
%     currentPos2(3)=currentPos2(3)-0.3;
currentPos2(4)=currentPos2(4); % up down

annotation(fig1,'textbox',...
currentPos2,...
'String',string,...
'FontSize',14,...
'FitBoxToText','on', 'LineStyle','none');
    
% add text for the Day-time and Night-time methods
string=strcat('Day-time method');
currentPos=get(gca,'Position');
currentPos2=currentPos;
currentPos2(1)=currentPos2(1)+0.4;  % left/right
currentPos2(2)=currentPos2(2);    
%     currentPos2(3)=currentPos2(3)-0.3;
currentPos2(4)=currentPos2(4); % up down

annotation(fig1,'textbox',...
currentPos2,...
'String',string,...
'FontSize',14,...
'FitBoxToText','on', 'LineStyle','none');


ylabel('Bias (%)')

set(a2,'fontsize',18)
set(a1,'fontsize',18)
 %xlim([0 numbars+1])
set(a1,'XTickLabel',xlabelz,'XTick',1:numbars)
set(a2,'XTickLabel',xlabelz,'XTick',1:numbars)
% %xticklabel_rotate([],45,[],'Fontsize',14)
  a1.XTickLabelRotation = 45;
  a2.XTickLabelRotation = 45;

if saveFigures==1
    filename=strcat('./figures/bias_barchart_notAbsolute');
    saveas(fig1,filename,'png');
    
end



