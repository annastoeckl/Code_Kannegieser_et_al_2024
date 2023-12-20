function plot_heat_map(plot_type,pattern_type)
% analyses the distribution of proboscis contacts on the flower, relative  to the pattern. 
% It plots a heatmap of proboscis contacts and generates histograms of the 
% contact distributions along the cardinal axes of the pattern, 
% and in the outer and inner thirds of the pattern. 
% It also calculates contact scores, which count and compare the number of proboscis contacts 
% along either cardinal axis of the flower and pattern, 
% as well as within the inner and outer thirds of the pattern. 
% For line patterns, it compares these scores to a hypothetical pattern 
% rotated 90 degrees from the originial, for cross patterns 45 degrees rotated.

%please select the plot type
% 'norm' = normalised (to max.) individual data is averaged
% 'sum' = all data is summed

%can plot different patterns:
%"line" patterns    which have either 4 points (single "line") or 12 points
%                   cross) and are labelled consecutively clockwise
%"broken" line      which consist of a single line, and a pattern2R, which
%                   is a second line indicating the gap
%circle             also indicated by 4 points, but these are not measured
%                   clockwise, but top, bottom, left, right (potentially
%                   the last two swapped) -- NEED to indicate for plot
%                   type! otherwise plots patterns as "line"

clear all, close all

if exist('plot_type')==0 || isempty(plot_type)
    plot_type='norm';
end

if exist('pattern_type')==0
    pattern_type='line';
end

%load data
% [filenames, path]=uigetfile('all*Data*.mat','multiselect','on');
[filenames, path]=uigetfile('*.mat','multiselect','on');
%if only one selected turn filenames into cell
if iscell(filenames)==0
    temp=filenames;clear filenames;
    filenames=cell(1,1);filenames{1}=temp;
end

scores=nan(length(filenames),4);%initialise scores
allMatrix=nan(200,200,length(filenames));
allMatrixNorm=nan(200,200,length(filenames));
patterns=cell(length(filenames),1);

for i=1:length(filenames)
    cd(path)
    load(filenames{i})
    %     disp(filenames{i})
    %load general parameters
    indTName=strfind(filenames{i},'_T');
    animal(i)=str2num(filenames{i}(indTName+2:find(filenames{i}(indTName+2:end)=='_',1,'first')+indTName));

    startInd=strfind(filenames{i},'_M');
    endInd=strfind(filenames{i}(startInd+2:end),'_');
    patterns{i}=filenames{i}(startInd+2:startInd+endInd(1));%this used to be str2num and double, but doesn't work for all patterns

    if sum(isnan(allprobR(:,1))==0)>50
        %delete nan entries from proboscis coordinates
        allprobR(isnan(allprobR(:,1)),:)=[];
        %convert all proboscis coordinates to the same framework
        %subtract midpoint
        temp=diff(outlineR);temp(2,:)=[];
        midpoint=[max(abs(temp(:,1)))/2+min(outlineR(:,1))  max(abs(temp(:,2)))/2+min(outlineR(:,2))];

        probNorm=allprobR-repmat(midpoint,size(allprobR,1),1);
        outlineNorm=outlineR-repmat(midpoint,size(outlineR,1),1);

        patternNorm=patternR-repmat(midpoint,size(patternR,1),1);

        %pattern2R is if we have 2 patterns (i.e. for a broken line: line
        %and "empty" part
        if ~exist('pattern2R'); pattern2R=[nan nan];end
        pattern2Norm=pattern2R-repmat(midpoint,size(pattern2R,1),1);

        %extra check to make sure pattern is also on 0 midpoint (which in
        %some cases it is not)
        temp=diff(patternNorm);temp(2,:)=[];midpointPattern=[max(abs(temp(:,1)))/2+min(patternNorm(:,1))  max(abs(temp(:,2)))/2+min(patternNorm(:,2))];
        if sum(midpoint)==0 && sum(midpointPattern)~=0
        patternNorm=patternNorm-repmat(midpointPattern,size(patternNorm,1),1);
        pattern2Norm=pattern2R-repmat(midpointPattern,size(pattern2R,1),1);
        end

        %load into matrix
        %add 100 to each dimension so that the origin is at (100 / 100) and we don't have negative values.
        %This looks good from eyeballing but has to change when the size of flower or field of view changes
        tempMatrix=nan(200,200,size(allprobR,1));
        for j=1:size(allprobR,1)
            if round(probNorm(j,2))+100<200 && round(probNorm(j,2))+100>0 && round(probNorm(j,1))+100 < 200 && round(probNorm(j,1))+100 >0 %make this for safety so points outside flower do not crash programme
                tempMatrix(round(probNorm(j,2))+100,round(probNorm(j,1))+100,j)=1;
            else
                tempMatrix(1,1,j)=nan;
            end
        end
        %sum over all proboscis touches. nansum turns those into 0. this is
        %good for Gaussian kernel, which otherwise doesnt work
        allMatrix(:,:,i)=imgaussfilt(nansum(tempMatrix,3),1);
        %normalise the matrix to 99% maximum, to not get too caught up with
        %outliers
        temp=allMatrix(:,:,i);scaling=quantile(temp(:),0.98);
        if scaling==0; scaling=1;end
        allMatrixNorm(:,:,i)=allMatrix(:,:,i)/scaling;
        %count entriess for each animal along cardinal directions and along
        %diagonals. rotate the image, to make sure that the same number of
        %pixels is sampled
        tempData=nansum(tempMatrix,3);
        tempData(tempData==0)=nan;

        %calculate the contacts in 1mm around the cardinal and diagonal axes.
        %for that, calculate how many pixels are 1mm
        %the outlineR, which is generally parallel to the coordinate system to start with, might become turned... so need to
        %take the euclidian distance between points rather than the actual distance
        temp1=sqrt((outlineNorm(1,1)-outlineNorm(2,1))^2+(outlineNorm(1,2)-outlineNorm(2,2))^2);
        temp2=sqrt((outlineNorm(3,1)-outlineNorm(4,1))^2+(outlineNorm(3,2)-outlineNorm(4,2))^2);
        lengthCircle=nanmean([temp1,temp2]);

        pixPerMM=lengthCircle/38;%38 mm flowers
        minInd=round(100-2*pixPerMM);maxInd=round(100+2*pixPerMM);

        %cardinal scores
        scores(i,1)=nansum(nansum(tempData([1:minInd maxInd:200],minInd:maxInd)));%X
        scores(i,2)=nansum(nansum(tempData(minInd:maxInd,[1:minInd maxInd:200])));%Y
        %rotate matrix by 45 degrees
        tempDataR=imrotate(tempData,45,'crop');%appends additional rows and cols
        scores(i,3)=nansum(nansum(tempDataR([1:minInd maxInd:200],minInd:maxInd)));%diag1
        scores(i,4)=nansum(nansum(tempDataR(minInd:maxInd,[1:minInd maxInd:200])));%diag2

        %pattern contact scores (not around axes but actually on pattern)
        %calculate contacts in pattern polygon
        patternIndex=checkCircle(patternNorm);
        [outlineCircle(:,1),outlineCircle(:,2)]=generateCircle(outlineNorm);
        polyShapePattern=polyshape(patternIndex(:,1),patternIndex(:,2));
        scoresA(i,1)=nansum(inpolygon(probNorm(:,1),probNorm(:,2),patternIndex(:,1),patternIndex(:,2)));
        scoresA(i,2)=nansum(inpolygon(probNorm(:,1),probNorm(:,2),outlineCircle(:,1),outlineCircle(:,2)))-scoresA(i,1);%all prob contacts - in pattern contacts
        scoresA(i,3)=area(polyshape(patternIndex(:,1),patternIndex(:,2)));%area pattern
        scoresA(i,4)=area(polyshape(outlineCircle(:,1),outlineCircle(:,2)))-scoresA(i,3);%area background

    end
end

%check if either pattern is a circle -- traced like outline, so 2nd entry
%has both x and y large step, unlike the line or cross stimuli
%if it is circle, generate circular pattern
patternFinal=checkCircle(patternNorm+100);%get patterns into right dimensions with +100
pattern2Final=checkCircle(pattern2Norm+100);


%calculate the circle of the flower
[xunitC,yunitC]=generateCircle(outlineNorm+100)


%% plot heatmap
f1=figure;hold on;
cmap=colormap('inferno');
if strcmp(plot_type,'norm')
    sumMatrixData=nanmean(allMatrixNorm,3); %this is averaging over normed entries
    imagesc(sumMatrixData,[0 1]);
elseif strcmp(plot_type,'sum')
    sumMatrixData=nansum(allMatrix,3); %this is summing or averaging over all
    imagesc(sumMatrixData,[0 quantile(sumMatrixData(:),.98)]);
end

plot(patternFinal(:,1),patternFinal(:,2),'w')
plot(pattern2Final(:,1),pattern2Final(:,2),'w')
plot(xunitC,yunitC,'w')
xlim([0 200]);ylim([0 200])
ylabel('pixels');xlabel('pixels');
axis square
colorbar

print(f1,'heatmap.eps','-dpdf','-r300','-painters','-bestfit')

%% plot scores
% if ~strcmp(pattern_type,'circle')
    %for analysis of CROSS patterns, both arms should be the same, therefore
    %add together

    if (size(unique(patternFinal(:,1)),1)>4 && size(unique(patternFinal(:,1)),1)<20) || size(unique(patternFinal(:,1)),1)==1 %this if the case for the cross pattern
        finalScores=[scores(:,1)+scores(:,2) scores(:,3)+scores(:,4)];
        scoreLabel={'cardinal','diagonal'};

    %for analysis of the CIRCLE pattern, use the contacts in circle vs
    %background, scale by the area of the circle foreground and background
    elseif size(unique(patternFinal(:,1)),1)>20  %this if the case for the circle pattern
        finalScores=[scoresA(:,1)./scoresA(:,3) scoresA(:,2)./scoresA(:,4)];
        scoreLabel={'foregr','backgr'};

    else

        %for analysis of single LINE patterns, there is only a vertical arm, and it
        %should be compared to the horizontal probes, not the diagonal
        finalScores=[scores(:,1) scores(:,2)];
        scoreLabel={'vertical','horizontal'};

    end

    % plot the summed contacts for both flower axes

    %sum all contacts per animal along the X and Y-axis
    X_integrated = permute(nanmean(allMatrixNorm,2),[1 3 2]);
    Y_integrated = permute(nanmean(allMatrixNorm,1),[2 3 1]);

    %integrate only outer 3rd of flower, and inner third
    Y_integrated_out = permute(nanmean(allMatrixNorm([1:66,134:200],:,:),1),[2 3 1]);
    Y_integrated_in = permute(nanmean(allMatrixNorm([67:133],:,:),1),[2 3 1]);
    X_integrated_out = permute(nanmean(allMatrixNorm(:,[1:66,134:200],:),2),[1 3 2]);
    X_integrated_in = permute(nanmean(allMatrixNorm(:,[67:133],:),2),[1 3 2]);


    %histgram of contacts summed along the Y-axis
    f2=figure('Position',[500   200   1120   840]);

    subplot(2,2,1)
    plot(1:size(allMatrixNorm,2),Y_integrated);
    ylabel('avg. contacts Y-all/ pixel');

    subplot(2,2,2)
    meanData=nanmean(Y_integrated,2);errorData=nanstd(Y_integrated,1,2)/sqrt(size(allMatrixNorm,3));
    shadedErrorBar_anna(1:size(allMatrixNorm,2),meanData./nanmax(meanData),errorData./nanmax(meanData))
    ylabel('avg. contacts Y-all/ pixel');

    subplot(2,2,3)
    meanData=nanmean(Y_integrated_out,2);errorData=nanstd(Y_integrated_out,1,2)/sqrt(size(allMatrixNorm,3));
    shadedErrorBar_anna(1:size(allMatrixNorm,2),meanData./nanmax(meanData),errorData./nanmax(meanData))
    ylabel('avg. contacts Y-out/ pixel');

    subplot(2,2,4)
    meanData=nanmean(Y_integrated_in,2);errorData=nanstd(Y_integrated_in,1,2)/sqrt(size(allMatrixNorm,3));
    shadedErrorBar_anna(1:size(allMatrixNorm,2),meanData./nanmax(meanData),errorData./nanmax(meanData))
    ylabel('avg. contacts Y-in/ pixel');

    print(f2,'contacts_X.eps','-dpdf','-r300','-painters','-bestfit')

    % plot the pattern contacts as line histogram collapsed on y axis
    %(integrated along X)
    f3=figure('Position',[500   200   1120   840]);

    subplot(2,2,1)
    plot(1:size(allMatrixNorm,2),X_integrated);
    ylabel('avg. contacts X-all/ pixel');

    subplot(2,2,2)
    meanData=nanmean(X_integrated,2);errorData=nanstd(X_integrated,1,2)/sqrt(size(allMatrixNorm,3));
    shadedErrorBar_anna(1:size(allMatrixNorm,2),meanData./nanmax(meanData),errorData./nanmax(meanData))
    ylabel('avg. contacts X-all/ pixel');

    subplot(2,2,3)
    meanData=nanmean(X_integrated_out,2);errorData=nanstd(X_integrated_out,1,2)/sqrt(size(allMatrixNorm,3));
    shadedErrorBar_anna(1:size(allMatrixNorm,2),meanData./nanmax(meanData),errorData./nanmax(meanData))
    ylabel('avg. contacts X-out/ pixel');

    subplot(2,2,4)
    meanData=nanmean(X_integrated_in,2);errorData=nanstd(X_integrated_in,1,2)/sqrt(size(allMatrixNorm,3));
    shadedErrorBar_anna(1:size(allMatrixNorm,2),meanData./nanmax(meanData),errorData./nanmax(meanData))
    ylabel('avg. contacts X-in/ pixel');

    print(f3,'contacts_Y.eps','-dpdf','-r300','-painters','-bestfit')

    %plot only the x-contacts all summed for paper
    f22=figure('position',[300 200 500 600]);
    meanData=nanmean(Y_integrated,2);errorData=nanstd(Y_integrated,1,2)/sqrt(size(allMatrixNorm,3));
    shadedErrorBar_anna([1:size(allMatrixNorm,2)]/pixPerMM-0.5*size(allMatrixNorm,2)/pixPerMM,meanData./nanmax(meanData),errorData./nanmax(meanData))
    ylabel('avg. contacts Y-all/ pixel');xlabel('position (mm)');
    xlim([-20 20]);ylim([0 1.1]);

    %get the half width
    indHalf1=nan(size(Y_integrated,2),1);indHalf2=nan(size(Y_integrated,2),1);
    for n=1:size(Y_integrated,2)
        temp=Y_integrated(:,n)/nanmax(Y_integrated(:,n));
        if ~nansum(~isnan(Y_integrated(:,n)))==0
            indHalf1(n)=find(temp>0.5,1,'first');
            indHalf2(n)=find(temp(find(temp==1,1,'first'):end)<0.5,1,'first')+find(temp==1,1,'first');
        else
            indHalf1(n)=nan;indHalf2(n)=nan;
        end
    end
    halfWidth=(indHalf2-indHalf1)/pixPerMM;

    %get the half width of the inner contacts
    indHalf1=nan(size(Y_integrated_in,2),1);indHalf2=nan(size(Y_integrated_in,2),1);
    for n=1:size(Y_integrated_in,2)
        temp=Y_integrated_in(:,n)/nanmax(Y_integrated_in(:,n));
        if ~nansum(~isnan(Y_integrated_in(:,n)))==0
            indHalf1(n)=find(temp>0.5,1,'first');
            indHalf2(n)=find(temp(find(temp==1,1,'first'):end)<0.5,1,'first')+find(temp==1,1,'first');
        else
            indHalf1(n)=nan;indHalf2(n)=nan;
        end
    end
    halfWidth_in=(indHalf2-indHalf1)/pixPerMM;

    print(f22,'HW_contacts_Xall.eps','-dpdf','-r300','-painters','-bestfit')

    %% plot pattern contact scores

    f4=figure('Position',[500 200 800 400]);
    subplot(1,2,1)
    if size(finalScores,1)>1
        boxplot(finalScores);hold on;
        plot(1+0.01*randn(1,size(finalScores,1)),finalScores(:,1),'k.','Markersize',10);
        plot(2+0.01*randn(1,size(finalScores,1)),finalScores(:,2),'k.','Markersize',10);
        set(gca,'xticklabel',scoreLabel)
    else
        plot([1 2],finalScores,'ro');xlim([0 3]);
        set(gca,'xtick',[1 2],'xticklabel',scoreLabel)
    end
    [p,h,stats]=signrank(finalScores(:,1),finalScores(:,2));
    title(sprintf('n=%2.0f, p=%1.3f, signrank=%1.3u',sum(isnan(finalScores(:,1))==0),p,stats.signedrank));
    ylabel('summed contacts');
    %     ylim([0 450])

    subplot(1,2,2)
    if size(finalScores,1)>1
        boxplot(finalScores);hold on;
        for u=1:size(finalScores,1)
            plot([1+0.01*randn(1,1),2+0.01*randn(1,1)],[finalScores(u,1),finalScores(u,2)],'k.-','Markersize',10);
        end
        set(gca,'xticklabel',scoreLabel)
    else
        plot([1 2],finalScores,'ro');xlim([0 3]);
        set(gca,'xtick',[1 2],'xticklabel',scoreLabel)
    end
    ylabel('summed contacts');
    %     ylim([0 450])

    print(f4,'contactScores_abs.eps','-dpdf','-r300','-painters','-bestfit')

    %% plot pattern contact scores relative to "background"

    finalScores(finalScores==0)=1;
    relScores=finalScores./repmat(finalScores(:,2),1,size(finalScores,2));
    ylimRoundScaleMax=round(max(relScores(:,1)));
    ylimRoundScaleMin=ylimRoundScaleMax*10;

    f5=figure('Position',[500 200 400 600]);
    subplot(1,2,1)
    if size(relScores,1)>1
        boxplot(relScores);hold on;
        for u=1:size(relScores,1)
            plot([1+0.01*randn(1,1),2+0.01*randn(1,1)],[relScores(u,1),relScores(u,2)],'k.-','Markersize',10);
        end
        set(gca,'xticklabel',scoreLabel)
    else
        plot([1 2],relScores,'ro');xlim([0 3]);hold on;
        set(gca,'xtick',[1 2],'xticklabel',scoreLabel)
    end
    ylabel('rel. contact time');
    plot([0 3],[1 1],'-k');

    [p,h,stats]=signrank(relScores(:,1),1);
    title(sprintf('n=%2.0f, p=%1.3f, signrank=%1.3u',sum(isnan(relScores(:,1))==0),p,stats.signedrank));
    ylabel('rel. contact time');
    ylim([0 ceil(max(relScores(:,1))/ylimRoundScaleMax)*ylimRoundScaleMax])

    subplot(1,2,2)
    if size(relScores,1)>1
        boxplot(relScores);hold on;
        for u=1:size(relScores,1)
            plot([1+0.01*randn(1,1),2+0.01*randn(1,1)],[relScores(u,1),relScores(u,2)],'k.-','Markersize',10);
        end
        set(gca,'xticklabel',scoreLabel)
    else
        plot([1 2],relScores,'ro');xlim([0 3]);
        set(gca,'xtick',[1 2],'xticklabel',scoreLabel)
    end
    ylabel('rel. contact time');
    plot([0 3],[1 1],'-k');
    set(gca,'yscale','log')

    [p,h,stats]=signrank(relScores(:,1),1);
    title(sprintf('n=%2.0f, p=%1.3f, signrank=%1.3u',sum(isnan(relScores(:,1))==0),p,stats.signedrank));
    ylabel('rel. contact time');
    set(gca,'yscale','log')
    ylim([floor(min([1; relScores(:,1)])*ylimRoundScaleMin)/ylimRoundScaleMin ceil(max(relScores(:,1))/ylimRoundScaleMax)*ylimRoundScaleMax])
    %     ylim([0 450])

    print(f5,'contactScores_rel.eps','-dpdf','-r300','-painters','-bestfit')

    %% perform a X square test for every individual to see if they significantly touched along the axis
    % for i=1:size(finalScores,1)
    %     [h,chiXp(i), chi2stat,df] = prop_test([finalScores(i,1) finalScores(i,1)] , [sum(finalScores(i,:)) sum(finalScores(i,:))], false);
    % end

    %% data summary
    % patternName=filenames{i}(strfind(filenames{i},'M'):find(filenames{i}(strfind(filenames{i},'M'):end)=='_',1,'first')+strfind(filenames{i},'M')-2);
    % save(['heatmap_data_',patternName,'.mat'],'animal','finalScores','patterns','chiXp');

    %% save data in excel files
    currPath=pwd;inds=strfind(currPath,'\');

    %all contact scores along entire pattern and 90 degree shifted axis
    all_contacts_tbl=table(animal',finalScores(:,1),finalScores(:,2),'VariableNames',{'animal','vertContacts','horzContacts'})
    writetable(all_contacts_tbl,[currPath(inds(end)+1:end),'_allContacts.csv'])
    % all_contacts_data=all_contacts_tbl{:,:};save([currPath(inds(end)+1:end),'_allContacts.mat'],'all_contacts_data');

    % contacts along the Z axis summed over collapsed Y axis of the flower in 2 outer and middle third
    %sum all contacts per animal along the X and Y-axis - not normalised
    %matrix
    X_integrated_sum = permute(nansum(allMatrix,2),[1 3 2]);
    Y_integrated_sum = permute(nansum(allMatrix,1),[2 3 1]);

    thirdLength=round(ceil(lengthCircle)/3);indThird=[1:thirdLength:3*thirdLength+1];%calculate thirds
    indThird=indThird+(200-indThird(end))/2;%center third indices
    for l=2:length(indThird)
        X_third_scores(:,l-1)=round(nansum(Y_integrated_sum(indThird(l-1):indThird(l),:),1));%x axis scores come from integrating Y
        Y_third_scores(:,l-1)=round(nansum(X_integrated_sum(indThird(l-1):indThird(l),:),1));%and vice versa
    end


    %contacts in the centre of the pattern and on the line along the
    %X-axis, integrated over the middle third of the flower
    Y_integrated_sum = permute(nansum(allMatrixNorm([67:133],:,:),1),[2 3 1]);
    ind_mid_contr=[min(patternNorm(:,1)), 0, max(patternNorm(:,1))]+100;
    if sum(~isnan(ind_mid_contr))==3
        X_contrast_scores_midthird=round(Y_integrated_sum(round(ind_mid_contr),:));
        %save data in table
        X_contacts_contrast_tbl=table(animal',X_contrast_scores_midthird(1,:)',X_contrast_scores_midthird(2,:)',X_contrast_scores_midthird(3,:)','VariableNames',{'animal','l_line','mid','r_line'})
        writetable(X_contacts_contrast_tbl,[currPath(inds(end)+1:end),'_XContacts_contrast.csv'])
        % X_contacts_contrast_data=X_contacts_contrast_tbl{:,:};save([currPath(inds(end)+1:end),'_XContacts_contrast.mat'],'X_contacts_contrast_data');
    end
    %half width of contacts along the pattern

    HW_table=table(animal',halfWidth,halfWidth_in,'VariableNames',{'animal','halfWidth','halfWidth_inhalfWidth_in'});
    writetable(HW_table,[currPath(inds(end)+1:end),'_HWcontacts.csv'])

end


function patternFinal=checkCircle(patternNorm)
%check if pattern is a circle -- which has 1 entry, where the distance in x
%and y is large, while for oriented lines and crosses, always one of the
%dimensions is close to 0 - if circle generate circle pattern from entries
if sum(~isnan(patternNorm))>=2
    diffP(:,1)=diff(patternNorm(:,1));diffP(:,2)=diff(patternNorm(:,2));
    if size(patternNorm,1)==4 && nansum(mink(abs(diffP(:)),3))>6
        [xunitC,yunitC]=generateCircle(patternNorm)
        patternFinal=[xunitC' yunitC'];
    else
        patternFinal=[patternNorm;patternNorm(1,:)];%close the pattern
    end
else
    patternFinal=patternNorm;
end
end

function [xunitC,yunitC]=generateCircle(outlineData)
temp=diff(outlineData);temp(2,:)=[];
midpoint=[max(abs(temp(:,1)))/2+min(outlineData(:,1))  max(abs(temp(:,2)))/2+min(outlineData(:,2))];
[D, alpha]=euclid_dist(outlineData);
lengthCircle=nanmean(D([1,3],:));
th = 0:pi/50:2*pi;
xunitC = lengthCircle/2 * cos(th) + midpoint(1);
yunitC = lengthCircle/2 * sin(th) + midpoint(2);
end

% %check if pattern is a circle
% if strcmp(pattern_type,'circle')
%     temp=diff(patternFinal);temp(2,:)=[];
%     midpoint=[max(abs(temp(:,1)))/2+min(patternFinal(:,1))  max(abs(temp(:,2)))/2+min(patternFinal(:,2))];
%     [D, alpha]=euclid_dist(patternFinal);
%     lengthCircle=nanmean(D([1,3],:));
%     th = 0:pi/50:2*pi;
%     xunitC = lengthCircle/2 * cos(th) + midpoint(1);
%     yunitC = lengthCircle/2 * sin(th) + midpoint(2);
%     patternFinal=[xunitC' yunitC'];
% else
%     patternFinal=[patternFinal;patternFinal(1,:)];
%     pattern2Final=[pattern2Final;pattern2Final(1,:)];
% end