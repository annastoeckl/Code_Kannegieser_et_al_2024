%this plots all the flightORI_data from flight_orientation.m from individuals in the same
%condition (or at least conditions that CAN be pooled)


clear all
close all

%load data: rotated data
[filenames, path]=uigetfile('flightORI*.mat','multiselect','on');
%if only one selected turn filenames into cell
if iscell(filenames)==0
    temp=filenames;clear filenames;
    filenames=cell(1,1);filenames{1}=temp;
end

Fs = 200;            % Sampling frequency
rotation_snips_all=[];
heading_snips_all=[];
animalID_all=[];
snips_groups=[];
heading_snips_start=cell(length(filenames),1);
heading_snips_middle=cell(length(filenames),1);
heading_snips_end=cell(length(filenames),1);
rotSpeed_touchNot=nan(length(filenames),2);
animalID=nan(length(filenames),1);
patternID=nan(length(filenames),1);

%define bins for circular histograms
binInt=44.9999999;
circBins=(0:binInt:360)-binInt/2;%circular histogram bins - centered on 0
binCenters=round(circBins(1:end-1)+binInt/2);

for i=1:length(filenames)
    cd(path)
    load(filenames{i})
    disp(filenames{i})

    %extract general parameters
    indTName=strfind(filenames{i},'_T');
    animal=str2num(filenames{i}(indTName+2:find(filenames{i}(indTName+2:end)=='_',1,'first')+indTName));
    animalID(i)=animal;
    
    if size(unique(patternR(:,1)),1)>4 || sum(~isnan(patternR(:)))==0 %this if the case for the cross pattern or no pattern
patternID(i)=-1;
    else
        patternID(i)=+1;
    end

    if isempty(heading_snips)==0
    %this compares rotation speed while touching the flower and not touching
    rotSpeed_touchNot(i,:)=[nanmedian(abs(rotSpeed_smooth(isnan(allprobR(2:end,1))))) nanmedian(abs(rotSpeed_smooth(isnan(allprobR(2:end,1))==0)))];
    
    turnCumSum_all(i)=turnCumSum(1)/turnCumSum(2);
    
    firstHeading(i)=heading_snips(1,2);%this is the very first time the animal approaches the flower
    
    heading_snips_all=[heading_snips_all;heading_snips];
    animalID_all=[animalID_all;repmat(animalID(i),size(heading_snips,1),1)];

    %collect headings for each animal to display them as an average
    %distribution of the population (averaging the distribution for each
    %animal). This cannot read nan values.
    if numel(heading_snips(isnan(heading_snips(:,2))==0,2))>1
    heading_snips_start{i}=heading_snips(isnan(heading_snips(:,2))==0,2);
    else heading_snips_start{i}=[];
    end
        
    if numel(heading_snips(isnan(heading_snips(:,4))==0,4))>2
    heading_snips_middle{i}=heading_snips(isnan(heading_snips(:,4))==0,4);
    else heading_snips_middle{i}=[];
    end
    
    if numel(heading_snips(isnan(heading_snips(:,6))==0,6))>2
    heading_snips_end{i}=heading_snips(isnan(heading_snips(:,6))==0,6);
    else heading_snips_end{i}=[];
    end
    
    %remove extreme values and nan
    rotation_snips(rotation_snips>1000)=nan;
    rotation_snips(sum(isnan(rotation_snips),2)>0,:)=[];
    
    rotation_snips_all=[rotation_snips_all;rotation_snips];
    snips_groups=[snips_groups;i*ones(size(heading_snips,1),1)];
    end
end

%% PLOT
spread=0.2;
favCol=[0 0.4 0.8];

f1=figure;
indCols=[1 2];
for u=1:length(indCols)
    v=violin(rotSpeed_touchNot(:,indCols(u)),'x',u,'facecolor',favCol,'medc','k');
    n=numel(rotSpeed_touchNot(:,indCols(u)));
    plot(spread*rand(n,1)+u*ones(n,1)-0.5*spread,rotSpeed_touchNot(:,indCols(u)),'.','MarkerSize',12,'color','k');hold on
end
b=boxplot(rotSpeed_touchNot(:,indCols),'color','w','PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol','','labels',{'no contact','contact'},'LabelOrientation','horizontal');
ylabel('avg agular velocity (°/s)');

% 
% f2=figure('position',[300   500    1200    400]);
% subplot(1,3,1)
% polarhistogram(deg2rad(heading_snips_all(:,2)),deg2rad([0:10:360]),'facecolor',favCol)
% set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
% title('start')
% 
% subplot(1,3,2)
% polarhistogram(deg2rad(heading_snips_all(:,4)),deg2rad([0:10:360]),'facecolor',favCol)
% set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
% title('middle')
% 
% subplot(1,3,3)
% polarhistogram(deg2rad(heading_snips_all(:,6)),deg2rad([0:10:360]),'facecolor',favCol)
% set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
% title('end')
% 
%% heading of the very first approach of each animal
veryFirst=nan(length(heading_snips_start),1);
for k=1:length(heading_snips_start)
    try
    veryFirst(k)=heading_snips_start{k}(1);
    catch
      veryFirst(k)=nan;  
    end
end

f0 = figure;
obj1=CircHist(veryFirst,circBins, 'areAxialData', true, 'threshStatsPlot',0.01);
obj1.polarAxs.ThetaAxis.MinorTickValues = []; % remove dotted tick-lines

f0.Position([3,4]) = [500,400]; % Figure dimensions


f0=figure('Position',[500  500  850  400]);

subplot(1,3,1);
for n=1:size(heading_snips_all,1)
    polarplot([0; deg2rad(heading_snips_all(n,2))],[0; heading_snips_all(n,1)],'k');hold on;
end
set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
title('path vectors')

subplot(1,3,2);
for n=1:size(heading_snips_all,1)
    polarplot([0; deg2rad(heading_snips_all(n,4))],[0; heading_snips_all(n,3)],'k');hold on;
end
set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
title('path vectors')

subplot(1,3,3);
for n=1:size(heading_snips_all,1)
    polarplot([0; deg2rad(heading_snips_all(n,6))],[0; heading_snips_all(n,5)],'k');hold on;
end
set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
title('path vectors')

%could do the same for the second set of info in the data

f0.Position([3,4]) = [1200,600]; % Figure dimensions
print(f0,'apprMiddleEndHeading.eps','-dpdf','-r300','-painters')

    %sort data for statistical analysis
    if sum(patternID>0)==0
        error('sectors not defined yet')
      
        sectorRatio=.5;
    else
        sector_scores1=heading_snips_all(:,[2,4,6])>150 | heading_snips_all(:,[2,4,6])<-150;
        sector_scores2=heading_snips_all(:,[2,4,6])>0 & heading_snips_all(:,[2,4,6])<30;
        sector_scores3=heading_snips_all(:,[2,4,6])<0 & heading_snips_all(:,[2,4,6])>-30;
        sector_scores=(sector_scores1+sector_scores2+sector_scores3);
        sectorRatio=.33;
    end

    counts=nan(size(animalID,1),4);
    for n=1:length(animalID)
    counts(n,:)=[nansum(sector_scores(animalID_all==animalID(n),1)) ...
                nansum(sector_scores(animalID_all==animalID(n),2)) ...
                nansum(sector_scores(animalID_all==animalID(n),3)) nansum(animalID_all(~isnan(heading_snips_all(:,2)))==animalID(n))];
    end

    tb21 = table(counts(:,1),counts(:,4)-counts(:,1),counts(:,4),repmat(sectorRatio,size(counts,1),1),animalID,'VariableNames',{'patternCount','bgCount','allCounts','threshold','animalID'});
    tb22 = table(counts(:,2),counts(:,4)-counts(:,2),counts(:,4),repmat(sectorRatio,size(counts,1),1),animalID,'VariableNames',{'patternCount','bgCount','allCounts','threshold','animalID'});    
    tb23 = table(counts(:,3),counts(:,4)-counts(:,3),counts(:,4),repmat(sectorRatio,size(counts,1),1),animalID,'VariableNames',{'patternCount','bgCount','allCounts','threshold','animalID'});    
    
    currPath=pwd;inds=strfind(currPath,'\');
    writetable(tb21,[currPath(inds(end)+1:end),'_apprHeadingCounts.csv'])
    writetable(tb22,[currPath(inds(end)+1:end),'_midlHeadingCounts.csv'])
    writetable(tb23,[currPath(inds(end)+1:end),'_deprHeadingCounts.csv'])
    

%% 
f4=figure;
indCols=[1 3 5];
for u=1:length(indCols)
    v=violin(heading_snips_all(:,indCols(u)),'x',u,'facecolor',favCol,'medc','k');
    n=numel(heading_snips_all(:,indCols(u)));
    plot(spread*rand(n,1)+u*ones(n,1)-0.5*spread,heading_snips_all(:,indCols(u)),'.','MarkerSize',12,'color','k');hold on
    
end
b=boxplot(heading_snips_all(:,indCols),'color','w','PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol','','labels',{'start','middle','end'},'LabelOrientation','horizontal');
ylabel('mean vector length');

f5=figure;
%this calculates heading between 0 and 90 -- so that distribution is NOT
%centered on 0
% headingNew=abs(heading_snips_all);
% for t=1:size(headingNew,2);headingNew(headingNew>90)=180-headingNew(headingNew>90);end
%this calculates heading between 0 and +-90 -- so that distribution IS
%centered on 0
headingNew=heading_snips_all;
for t=1:size(headingNew,2);headingNew(headingNew>90)=180-headingNew(headingNew>90);
    headingNew(headingNew<-90)=180+headingNew(headingNew<-90);end

indCols=[2 4 6];
for u=1:length(indCols)
    v=violin(headingNew(:,indCols(u)),'x',u,'facecolor',favCol,'medc','k');
    n=numel(headingNew(:,indCols(u)));
    plot(spread*rand(n,1)+u*ones(n,1)-0.5*spread,headingNew(:,indCols(u)),'.','MarkerSize',12,'color','k');hold on
end
b=boxplot(headingNew(:,indCols),'color','w','PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol','','labels',{'start','middle','end'},'LabelOrientation','horizontal');
ylabel('mean heading (°)');
ylim([-90 90])

f6=figure;
indCols=[1 2 3];
for u=1:length(indCols)
    v=violin(rotation_snips_all(:,indCols(u)),'x',u,'facecolor',favCol,'medc','k');
    n=numel(rotation_snips_all(:,indCols(u)));
    plot(spread*rand(n,1)+u*ones(n,1)-0.5*spread,rotation_snips_all(:,indCols(u)),'.','MarkerSize',12,'color','k');hold on
end
b=boxplot(rotation_snips_all(:,indCols),'color','w','PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol','','labels',{'start','middle','end'},'LabelOrientation','horizontal');
ylabel('avg rotational velocity (°/s)');

f7=figure;hold on;
v=violin(turnCumSum_all','facecolor',favCol,'medc','k');hold on;
n=numel(turnCumSum_all);
plot(spread*rand(n,1)+1-0.5*spread,turnCumSum_all,'.','MarkerSize',12,'color','k');hold on
% b=boxplot(turnCumSum_all,'color','w','PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol','','LabelOrientation','horizontal');
ylabel('turning index (left>1, right<1)');
set(gca,'yscale','log')
ylim([0.9*nanmin(turnCumSum_all) 1.1*nanmax(turnCumSum_all)])

tbTurnInd=table(animalID,turnCumSum_all','VariableNames',{'animalID','turningIndex'})
writetable(tbTurnInd,[currPath(inds(end)+1:end),'_turnIndex.csv'])

%% plot simple histogram, rather than average histogram
%heading distribution for touchdown positions and touchdown headings
allsnips=[];
for u=1:length(heading_snips_start)
    allsnips=[allsnips;heading_snips_start{u}];
end

f5 = figure;
subplot(1,2,1)
heading_snips_touchd_vec=[];
h=polarhistogram(deg2rad(allsnips),deg2rad(circBins));
landing_pos=h.Values;

set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','clockwise');
title('approach headings')
 
%this is to make a dot distribution

% radhist=cdot_plot(heading_snips_touchd_vec+180,binInt);
%  cdot_label(0:30:330, radhist);

 

% f8=figure;hold on;
% %headin from -90 to 90
% headingNew=firstHeading;
% for t=1:size(headingNew,2);headingNew(headingNew>90)=180-headingNew(headingNew>90);
%     headingNew(headingNew<-90)=180+headingNew(headingNew<-90);end
%
% v=violin(headingNew','facecolor',favCol,'medc','k');hold on;
% n=numel(headingNew);
% plot(spread*randn(n,1)+1,headingNew,'.','MarkerSize',12,'color','k');hold on
% b=boxplot(headingNew,'color','w','PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol','','LabelOrientation','horizontal');
% ylabel('first touch heading (°)');
% ylim([-90 90])