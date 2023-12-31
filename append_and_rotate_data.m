% This script loads loads the tracked data from the .mat files generated 
% with load_excel_file_proboscis_tracks, adds together data from multiple files 
% (if multiple trials were recorded for the same animal and condition), 
% and rotates the data to align the pattern's main axis with the cardinal coordinates. 
% In addition to identifying trials from multiple files, it also identifies individual 
% approaches of the moths to the flower within one datafile (videofile), as flight bouts.
% it saves the output as allRotData_filename.mat, which provides the input
% for all subsequent analysis

clear all
close all

%load data
%select all data from the same animal same pattern (=multiple flights)
filenames=uigetfile('data*.mat','multiselect','on');

%if only one selected turn filenames into cell
if iscell(filenames)==0
    temp=filenames;clear filenames;
    filenames=cell(1,1);filenames{1}=temp;
end

allCounts=[];
chanceL=zeros(1,length(filenames));
dataR=[];

allprobR=[];%proboscis coordinates
allThoraxR=[];%thorax coordinates
allHeadR=[];%head coordinates
allDataR=[];%all coordinates
allR=[];%rotation coordinates for all datafiles

timeBouts=[]; %timespoints (beginning, end) of individual flight bouts

figure;hold on;
temp=colormap('jet');
cols=temp(1:ceil(length(temp)/length(filenames)):length(temp),:);

for j=1:length(filenames)
    load(filenames{j})
    
    %some older datafiles do not have 'dataMatrix', which combines all
    %labelled indices. in this case, re-create it from individual vectors
    temp=diff(outline);temp(2,:)=[];
    midpoint=[max(abs(temp(:,1)))/2+min(outline(:,1))  max(abs(temp(:,2)))/2+min(outline(:,2))];
    
    if exist('dataMatrix')==0
        labels=cell(1,1);labels{1}='Proboscis';
        dataMatrix=[prob];
        prob(prob<0)=nan;%remove negative entries
        if exist('head')==1 && exist('thorax')==1
            if sum(~isnan(head))>0 && sum(~isnan(thorax))>0 %addition for Alexas March data. often head was not digitized, only thorax
                dataMatrix=[dataMatrix thorax head];
                labels{2}='Thorax';labels{3}='Head';
                head(head<0)=nan;thorax(thorax<0)=nan;%remove negative entries
            else
                head=[nan nan;nan nan];
                thorax=[nan nan;nan nan];
            end
        else
            head=[nan nan;nan nan];
            thorax=[nan nan;nan nan];
        end
    end
    dataMatrix(dataMatrix<0)=nan;%remove negative entries
    
    
    %% rotate the pattern so that its arms align the the x-y axes
    %check if this is a circle pattern - which shouldn't be rotated
    %(because the circle is indicated by 4 points so rotation makes no
    %sense)
    indStart=strfind(filenames{j},'_M');
    indCircle=strfind(filenames{j}(indStart+1:find(filenames{j}(indStart+1:end)=='_',1,'first')+indStart+1),'C');
    
    %remove double entries from pattern pointer
    % [uniquePoints, ind1,ind2]=unique(round(cage/5)*5,'rows','stable');
    % cage=cage(ind1,:);%delete the overlapping entries
    
    if sum(isnan(pattern(:)))~=0 || ~isempty(indCircle)
        %this is the case if there was no pattern, so no rotation needed
        chanceL=nan;
        allCounts=nan;
        
        %move proboscis data and flower indicators to the same position
        %actually just put to origin! %because consecutive flights might
        %not have the flower turned the same way
        %midpoint of the flower needed, it is already calculated above
        
        %subtract from flower and proboscis data
        probR=prob-midpoint;
        patternR=pattern-midpoint;%need to subtract midpoint in case it is a circle
        if exist('pattern2'); pattern2R=pattern2-midpoint; else pattern2R=[nan nan];end
        outlineR=outline-midpoint;
        thoraxR=thorax-midpoint;
        headR=head-midpoint;
        dataR=dataMatrix-repmat(midpoint,1,size(dataMatrix,2)/2);
        
        
        %save appended proboscis and  head and thorax coordinates
        %!!!for now these are all separate but can eventually be used from
        %allDataR
        
        allThoraxR=[allThoraxR;thoraxR];
        allHeadR=[allHeadR;headR];
        allDataR=[allDataR;dataR];
        allR=nan;
        
        % remove proboscis entries outside of the flower face (which are generated by DLC)
        %allProbR
        temp=diff(outlineR);temp(2,:)=[];
        lengthCircle=nanmean(nanmax(abs(temp)),2);
        
        th = 0:pi/50:2*pi;
        xunitC = lengthCircle/2 * cos(th);
        yunitC = lengthCircle/2 * sin(th);
        %     plot(xunitC,yunitC,'-');
        
        in = inpolygon(probR(:,1),probR(:,2),xunitC,yunitC);
        probR(in==0,:)=repmat([nan nan],sum(in==0),1);
        
        
        allprobR=[allprobR;probR];
        
        %plot control image
        %calculate the circle of the flower
        %         temp=diff(outlineR);temp(2,:)=[];
        %         midpoint=[max(abs(temp(:,1)))/2+min(outlineR(:,1))  max(abs(temp(:,2)))/2+min(outlineR(:,2))];
        %         lengthCircle=nanmean(nanmax(abs(temp)),2);
        %
        %         th = 0:pi/50:2*pi;
        %         xunitC = lengthCircle/2 * cos(th) + midpoint(1);
        %         yunitC = lengthCircle/2 * sin(th) + midpoint(2);
        %
        %         plot(probR(:,1),probR(:,2),'*','color',cols(j,:));
        %         plot(xunitC,yunitC,'-','color',cols(j,:));
    else
        
        %% pattern rotation
        % use each consecutive point on the pattern to calculate rotation
        %double the first entry in the end, to use the entire pattern for
        %this, and not miss the last side
        anglePattern=[pattern;pattern(1,:)];
        
        %use only the long sides for rotation
        temp=sum(abs(diff(anglePattern,1)),2);
        
        for i=1:size(anglePattern,1)-1
            alpha(i)=atand((anglePattern(i,2)-anglePattern(i+1,2))/(anglePattern(i,1)-anglePattern(i+1,1)));
        end
        
        %this is to remove the shorter sides of the pattern from
        %calculation. but cannot just remove short sides, or will exclude
        %smaller patterns (like 2x2mm)
        alpha(temp<0.5*max(temp))=[];
        
        %put all angles in a 90 degree framework (mod 90 doesnt work, because it
        %turns small negative angles into large pos ones
        %         alpha(1)=[];%take out the first value, can be compromised
        alpha(abs(alpha)>45)=mod(alpha(abs(alpha)>45),90);
        alpha(abs(alpha)>45)=alpha(abs(alpha)>45)-90;
        
        %turn the "faster" way
        if nanmean(alpha)>45
            alpha=alpha-90;
        end
        
        %rotation matrix
        R=[cosd(-nanmean(alpha)) -sind(-nanmean(alpha)); sind(-nanmean(alpha)) cosd(-nanmean(alpha))];
        
        %rotate all coordinates
        patternR=[R*pattern']';
        if exist('pattern2'); pattern2R=[R*pattern2']'; else pattern2R=[nan nan];end
        probR=[R*prob']';
        outlineR=[R*outline']';
        thoraxR=[R*thorax']';
        headR=[R*head']';
        for u=1:size(dataMatrix,2)/2
            dataR(:,u*2-1:u*2)=[R*dataMatrix(:,u*2-1:u*2)']';
        end
        
        %check if pattern is turned so that longer side is in line with Y
        %axis, if not, rotate another 90 degrees to do so
        check=max(abs(diff(patternR,1)));
        if check(1)>check(2)
            alpha=90;
            R90=[cosd(-nanmean(alpha)) -sind(-nanmean(alpha)); sind(-nanmean(alpha)) cosd(-nanmean(alpha))];
            patternR=[R90*patternR']';
            pattern2R=[R90*pattern2R']';
            probR=[R90*probR']';
            outlineR=[R90*outlineR']';
            thoraxR=[R90*thoraxR']';
            headR=[R90*headR']';
            for u=1:size(dataMatrix,2)/2
                dataR(:,u*2-1:u*2)=[R90*dataR(:,u*2-1:u*2)']';
            end
        end
        
        % in some videos, the outline of the flower was not oriented the
        % same way as the pattern. define the length of the circle with euclidian distance between outline points and
        % not with absolute distance along axes. also, be aware those
        % points are not rotated along same axes as pattern necessarily.
        
        
        %move data indicators to the same position
        %actually just put to origin! %because consecutive flights might
        %not have the flower turned the same way
        %calculate the midpoint of the flower
        temp=diff(outlineR);temp(2,:)=[];
        midpoint=[max(abs(temp(:,1)))/2+min(outlineR(:,1))  max(abs(temp(:,2)))/2+min(outlineR(:,2))];
        
        %subtract from flower and proboscis data
        probR=probR-midpoint;
        patternR=patternR-midpoint;
        pattern2R=pattern2R-midpoint;
        outlineR=outlineR-midpoint;
        thoraxR=thoraxR-midpoint;
        headR=headR-midpoint;
        dataR=dataR-repmat(midpoint,1,size(dataMatrix,2)/2);
        
        %save appended proboscis
        %!!!for now these are all separate but can eventually be used from
        %allDataR
        
        allThoraxR=[allThoraxR;thoraxR];
        allHeadR=[allHeadR;headR];
        allR=[allR;[R(1,:) size(probR,1)]];%save rotation coordinate and index where rotation happened
        allDataR=[allDataR;dataR];
        
        
        
        %% remove proboscis entries outside of the flower face (which are generated by DLC)
        %allProbR - since all data points are rotated according to the
        %pattern, the outlineR, which is generally parallel to the
        %coordinate system to start with, might become turned... so need to
        %take the euclidian distance between points rather than the actual
        %distance
        temp1=sqrt((outlineR(1,1)-outlineR(2,1))^2+(outlineR(1,2)-outlineR(2,2))^2);
        temp2=sqrt((outlineR(3,1)-outlineR(4,1))^2+(outlineR(3,2)-outlineR(4,2))^2);
        lengthCircle=nanmean([temp1,temp2]);
        
        th = 0:pi/50:2*pi;
        xunitC = lengthCircle/2 * cos(th);
        yunitC = lengthCircle/2 * sin(th);
        %     plot(xunitC,yunitC,'-');
        
        in = inpolygon(probR(:,1),probR(:,2),xunitC,yunitC);
        probR(in==0,:)=repmat([nan nan],sum(in==0),1);
        
        
        allprobR=[allprobR;probR];
    end
    
    %% find individual flight bouts (when pause was used to record several approaches in one file)
    %use the head and thorax data if available
    if size(thoraxR,1)>2 %if there is data in the thorax file
        temp=diff(isnan(thoraxR(:,1)));%find blocks of data and nan to find out where the flight bouts are.
        startInds=find(temp==-1);
        endInds=find(temp==1);
        if numel(endInds)<numel(startInds) %if the end is missing because the file ends with data
            endInds=[endInds;size(thoraxR,1)];
        elseif numel(endInds)>numel(startInds) %if the beginning is missing because file starts with data
            startInds=[1;startInds];
        end
        
        %check how close the start and end inds are together. if less than
        %10 frames, make 1 bout out of it.
        if numel(startInds)>1
        indDiff=startInds(2:end)-endInds(1:end-1);
        startInds([false; indDiff<10])=[];
        endInds([indDiff<10;false])=[];
        end
        
        %check also if coordinates are jumping within these bouts. This could also be an indicator
        %for cut together videos.
        for o=1:numel(startInds)
            %check for jumps in both x and y coordinates
            temp=abs(diff(thoraxR(startInds(o):endInds(o),1)))>50+abs(diff(thoraxR(startInds(o):endInds(o),2)))>50;
            if sum(temp(:))>0
                disp('there are data sequences clipped together!')
            end
        end
    else
        temp=diff(isnan(probR(:,1)));%if only proboscis data exists take the entire block of data (hard to differentiate within)
        startInds=find(temp==-1,1,'first');
        endInds=find(temp==1,1,'last');
        if numel(endInds)<numel(startInds) %if the end is missing because the file ends with data
            endInds=size(probR,1);
        elseif numel(endInds)>numel(startInds) %if the beginning is missing because file starts with data
            startInds=1;
        elseif isempty(startInds) && isempty(endInds)
            startInds=1;endInds=size(probR,1);
        end
    end
    startingTime=size(allprobR,1)-size(probR,1);%get the starting time frame when not in first data iteration
    timeBouts=[timeBouts;startInds+startingTime endInds+startingTime];
    
    
    %% remove "stationary" points, where the tracking algorithm seems to fill
    % in points when the animal leaves the frame
    
    %use all Data to find these, and the head and thorax points
    statDataInds=diff(round(allDataR),1)==0;
    if sum(isnan(head(:)))<numel(head) && sum(isnan(thorax(:)))<numel(thorax)
        %go through data and count all stretches of no difference and remember
        %indices
        count=0;bodyInd=3;%index for which body part to use for this
        indEnd=[];indStart=[];statLength=[];
        for t=2:size(statDataInds,1)
            if statDataInds(t,bodyInd)==1 && statDataInds(t,bodyInd)-statDataInds(t-1,bodyInd)==0
                if t==2
                    count=count+1;
                    indStart(count)=1;
                    statLength(count)=0;
                end
                statLength(count)=statLength(count)+1;
            elseif statDataInds(t,bodyInd)-statDataInds(t-1,bodyInd)==-1
                %then stationary strip ends
                count=count+1;
                indEnd(count)=t;
            elseif statDataInds(t,bodyInd)-statDataInds(t-1,bodyInd)==1
                %then stationary strip starts
                count=count+1;
                indStart(count)=t;
                statLength(count)=0;
                
            end
        end
        
        if length(indEnd)<length(indStart)
            indEnd(end+1)=length(statDataInds);
        end
        
        %delete all stationary entries shorter than 20 frames
        indStart(statLength<20)=[];indEnd(statLength<20)=[];statLength(statLength<20)=[];
        %make the stationary parts nan
        if isempty(indStart)~=0
            for t=1:length(indStart)
                allDataR(indStart(t):indEnd(t),:)=nan;
                allprobR(indStart(t):indEnd(t),:)=nan;
                allHeadR(indStart(t):indEnd(t),:)=nan;
                allThoraxR(indStart(t):indEnd(t),:)=nan;
            end
        end
        
    end
    
    %% plot control image
    
    subplot(1,3,1);hold on;
    plot(probR(:,1),probR(:,2),'*','color',cols(j,:));
    plot([patternR(:,1);patternR(1,1)],[patternR(:,2);patternR(1,2)],'-','color',cols(j,:));
    plot([pattern2R(:,1);pattern2R(1,1)],[pattern2R(:,2);pattern2R(1,2)],'-','color',cols(j,:));
    %calculate the circle of the flower
    temp=diff(outlineR);temp(2,:)=[];
    midpoint=[max(abs(temp(:,1)))/2+min(outlineR(:,1))  max(abs(temp(:,2)))/2+min(outlineR(:,2))];
    %calculate circle length as absolute difference between outline points. Because they might not be rotated according to
    %cardinal axes
    [D, alpha]=euclid_dist(outlineR);
    lengthCircle=nanmean(D([1,3],:));
    
    th = 0:pi/50:2*pi;
    xunitC = lengthCircle/2 * cos(th) + midpoint(1);
    yunitC = lengthCircle/2 * sin(th) + midpoint(2);
    plot(xunitC,yunitC,'-','color',cols(j,:));
    xlim([-100 100]);ylim([-100 100]);
    axis square
    title('proboscis')
    
    subplot(1,3,2);hold on;
    plot(thoraxR(:,1),thoraxR(:,2),'o','color',cols(j,:));hold on;
    plot([patternR(:,1);patternR(1,1)],[patternR(:,2);patternR(1,2)],'-','color',cols(j,:));
    plot([pattern2R(:,1);pattern2R(1,1)],[pattern2R(:,2);pattern2R(1,2)],'-','color',cols(j,:));
    plot(xunitC,yunitC,'-','color',cols(j,:));
    xlim([-300 300]);ylim([-300 300]);
    axis square
    title('thorax')
    
    subplot(1,3,3);hold on;
    plot(headR(:,1),headR(:,2),'*','color',cols(j,:));hold on;
    plot([patternR(:,1);patternR(1,1)],[patternR(:,2);patternR(1,2)],'-','color',cols(j,:));
    plot([pattern2R(:,1);pattern2R(1,1)],[pattern2R(:,2);pattern2R(1,2)],'-','color',cols(j,:));
    plot(xunitC,yunitC,'-','color',cols(j,:));
    xlim([-300 300]);ylim([-300 300]);
    axis square
    title('head')
    
    clear dataR probR headR thoraxR
end

%display number of contacts
disp(['number of contacts: ',num2str(sum(isnan(allprobR(:,1))==0))]);
%display number of flight bouts and length
disp(['number of flight bouts: ',num2str(size(timeBouts,1))]);
flight_bout_length=(timeBouts(:,2)-timeBouts(:,1))/200;%in s, provided Fs = 200
t=table([1:size(timeBouts,1)]',flight_bout_length,'variableNames',{'bout_nr','bouth_length_s'})


indStart=strfind(filenames{j},'_T');
save(['allRotData_',filenames{1}(indStart+1:end)],'allprobR','allThoraxR','allHeadR','patternR','pattern2R','outlineR','allR','allDataR','labels','timeBouts')