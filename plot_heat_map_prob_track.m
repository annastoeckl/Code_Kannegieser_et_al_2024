function plot_heat_map_prob_track(plot_type,pattern_type)

% analyses the distribution of proboscis contacts on the flower, relative  to the pattern. 
% It plots a heatmap of proboscis contacts and the paths as line plots. 
% It also calculates contact scores, which count and compare the number of proboscis contacts 
% along either cardinal axis of the flower and pattern, and within the pattern and on the 
% background, and includes scales for the area of pattern and background. 
% From the proboscis tracks, the duration, speed and tortuosity of probing are extracted and plotted, 
% as well as the mean direction and vector length of the probing tracks, relative to the orientation of the pattern.
% Moreover, the body orientation (defined by the head-thorax axis) and position of the proboscis, 
% head and thorax relative to the long axis of a line pattern are extracted and plotted.

%please select the plot type
% 'norm' = normalised (to max.) individual data is averaged
% 'sum' = all data is summed

%can plot different patterns:
%"line" patterns    which have either 4 points (single "line") or 12 points
%                   cross) and are labelled consecutively clockwise
%"broken" line      which consist of a single line, and a pattern2R, which
%                   is a second line indicating the gap
%circle             also indicated by 4 pointheats, but these are not measured
%                   clockwise, but top, bottom, left, right (potentially
%                   the last two swapped) -- NEED to indicate for plot
%                   type! otherwise plots patterns as "line"
%use interpolated data for contact scores!!!
close all

if exist('plot_type')==0 || isempty(plot_type)
    plot_type='norm';
end

if exist('pattern_type')==0
    pattern_type='line';
end

%load data
[filenames, path]=uigetfile('*ata*.mat','multiselect','on');%load all "data" formats (normal and allRotData)
%if only one selected turn filenames into cell
if iscell(filenames)==0
    temp=filenames;clear filenames;
    filenames=cell(1,1);filenames{1}=temp;
end

%sample rate
Fs=200;%hz - should be true for all videos
rad_patternCount=2;%the radius of the corridor around the cardinal axes of the pattern for checking pattern contacts
gaussFiltKernel=1;
scores=nan(length(filenames),4);%initialise scores
allMatrix=nan(200,200,length(filenames));
allMatrixNorm=nan(200,200,length(filenames));
patterns=cell(length(filenames),1);
colourInd=cell(length(filenames),1);
%make colour code for each bee
col=hsv;
ind_col=interp1(col,size(col,1)/length(filenames):size(col,1)/length(filenames):size(col,1),[],'extrap');
ind_col(ind_col<0)=0;
%spread of boxplot original datapoints and colour of violins
spread=0.2;
favCol=[0 0.4 0.8];

%summary data matrices
path_dur=[];%approach, depr, individualNr - walking time in ms
path_tort=[];%approach, depr, individualNr - tortuosity
path_length=[];%approach, depr, individualNr - path length in px
path_r=[];%approach, individualNr - path mean vec length
path_dir=[];%approach, individualNr - path mean dir
path_speed=[];%average speed of proboscis movement
allAMatrix=[];%all entries for each flight bout of contacts on flower appr
allScoresA=[];
allScoresLine=[];
allXYScoresA=[];
animalID_all=[];%collects all animal IDs for each label
animalID=ones(length(filenames),1);
label=[]; %collects all filenames for each data entry
colourID_all=[];
allBodyAngles=[];
allxBodyPos=[];
allxProbHist=[];
probCenterDistAngl=cell(length(filenames),1);%collect distance and angle from centre over time for all animals
PCA=nan(length(filenames),8);

%initialise video


f1=figure('Position',[500  500  850  400]);hold on;

for i=1:length(filenames)
    cd(path)
    load(filenames{i})
    disp(['processing ',filenames{i}])
    %load general parameters
    indTName=strfind(filenames{i},'_T');
    animalID(i)=str2num(filenames{i}(indTName+2:find(filenames{i}(indTName+2:end)=='_',1,'first')+indTName));
    
    startInd=strfind(filenames{i},'_M');
    endInd=strfind(filenames{i}(startInd+2:end),'_');
    patterns{i}=filenames{i}(startInd+2:startInd+endInd(1));%this used to be str2num and double, but doesn't work for all patterns
    colourInd{i}=str2num(patterns{i});
    if isempty(colourInd{i});colourInd{i}=99;end %if there is no number for pattern width, give own category of pattern
    
    if sum(isnan(allprobR(:,1))==0)>50
        
        %find start and end points of each proboscis contact within each
        %flight bout (timeBouts)
        %touchBouts when proboscis had contact
        touchBouts=nan(size(timeBouts,1),2);
        for j=1:size(timeBouts,1)
            %if no feedBouts, try to generate them here
            if sum(isnan(allprobR(timeBouts(j,1):timeBouts(j,2),1)))>=2
                touchBouts(j,1)=find(isnan(allprobR(timeBouts(j,1):timeBouts(j,2),1)),1,'first')+timeBouts(j,1);
                touchBouts(j,2)=find(isnan(allprobR(timeBouts(j,1):timeBouts(j,2),1)),1,'last')+timeBouts(j,1);
            else %then the animal was tracked entirely, without indicating where it fed
                touchBouts(j,:)=[nan nan];
            end
        end
        
        %bee nomenclature has timeBout when bee lands, and feedBout when
        %bee feeds nectar. hawkmoth experiments don't have nextar, so
        %timeBouts are each approach flight, and touchBouts are when
        %proboscis has contact.
        %transfer into nomensclature by using timeBout for all contact
        %time, and feedBout for the whole time the proboscis is on the
        %flower (so feeding right at the end of the timebout)
        oldTimeBouts=timeBouts;%preserve original
        timeBouts=touchBouts;
        feedBouts=[timeBouts(:,2) timeBouts(:,2)+1];
        
        %convert all proboscis coordinates to the same framework
        %subtract midpoint
        temp=diff(outlineR);temp(2,:)=[];
        midpoint=[max(abs(temp(:,1)))/2+min(outlineR(:,1))  max(abs(temp(:,2)))/2+min(outlineR(:,2))];
        %         probNorm=allprobR-repmat(midpoint,size(allprobR,1),1);
        
        patternNorm=patternR-repmat(midpoint,size(patternR,1),1);
        pattern2Norm=pattern2R-repmat(midpoint,size(pattern2R,1),1);
        outlineNorm=outlineR-repmat(midpoint,size(outlineR,1),1);
        
        %calculate the dimensions of the flower pattern and the flower area
        %calculate how many pixels are 1mm
        %the outlineR, which is generally parallel to the coordinate system to start with, might become turned... so need to
        %take the euclidian distance between points rather than the actual distance
        [D, alpha]=euclid_dist(outlineR);
        lengthCircleFlower=nanmean(D([1,3],:));
        pixPerMM=lengthCircleFlower/38;%38 mm flowers
        
        %to remove proboscis entries outside of the flower face (which are generated by DLC)
        %take the euclidian distance between points rather than the actual
        %distance -- removed for now, as bbees dont do that much...
        temp1=sqrt((outlineNorm(1,1)-outlineNorm(2,1))^2+(outlineNorm(1,2)-outlineNorm(2,2))^2);
        temp2=sqrt((outlineNorm(3,1)-outlineNorm(4,1))^2+(outlineNorm(3,2)-outlineNorm(4,2))^2);
        lengthCircle=nanmean([temp1,temp2]);
        
        th = 0:pi/50:2*pi;
        xunitC = lengthCircle/2 * cos(th);
        yunitC = lengthCircle/2 * sin(th);
        
        %extract approaches and departures
        approaches=cell(size(timeBouts,1),1);
        plotTrackData=cell(size(timeBouts,1),1);
        tempAMatrix=nan(200,200,size(timeBouts,1));
        scoresA=nan(size(timeBouts,1),5);%scores cardinal diagonal position components
        scoresInside=nan(size(timeBouts,1),2);%inside line
        scoresOutside=nan(size(timeBouts,1),2);%outside line
        appr_length=nan(1,size(timeBouts,1));
        appr_tort=nan(1,size(timeBouts,1));
        appr_speed=nan(1,size(timeBouts,1));
        rmean=nan(2,size(timeBouts,1));
        dirmean=nan(2,size(timeBouts,1));
        xBodyPos_all=nan(size(timeBouts,1),6);%proboscis, head
        xProbHistInd=-35:1:35;
        xProbPosHist=nan(size(timeBouts,1),length(xProbHistInd));
        bodyAngle_all=nan(size(timeBouts,1),4);%approach angle and lector length r, average probing angle and r
        
        
        for u=1:size(timeBouts,1)
            if ~isnan(feedBouts(u,1))
                %this is a hack from the bumblebee script.
                
                %The probing was split into approach and departure (separated by feeding).
                %we use only "approach here for all, so set feeding to end of proboscis contacts in bout
                if feedBouts(u,1)<size(allprobR,1)
                    approaches{u}=allprobR(timeBouts(u,1):feedBouts(u,1),:)-midpoint;
                else
                    approaches{u}=[nan nan];
                end
                %remove points outside of circle
                in = inpolygon(approaches{u}(:,1),approaches{u}(:,2),xunitC,yunitC);
                if sum(in)~=0
                    approaches{u}(in==0,:)=repmat([nan nan],sum(in==0),1);%set all proboscis data outside of circle to 0
                    %calculate tortuosity (NOT normalised path step length)
                    tempAppr=approaches{u};tempAppr(isnan(tempAppr(:,1)),:)=[];
                    appr_length(u)=nansum(euclid_dist(tempAppr))/pixPerMM;
                    appr_tort(u)=nansum(euclid_dist(tempAppr))/euclid_dist(tempAppr(1,:),tempAppr(end,:));
                    appr_speed(u)=nanmedian(euclid_dist(tempAppr))/pixPerMM*Fs;
                    
                    
                    %calculate from interpolated tracks
                    %interpolate every single snippet, otherwise
                    %interpolation creates fake tracks when data jumps
                    if round(size(approaches{u},1)/5)>=2
                        %create vector index of continuous segments of proboscis probing
                        startInd=find(diff(isnan(approaches{u}(:,1)))==-1);
                        if isempty(startInd);startInd=1;end
                        endInd=find(diff(isnan(approaches{u}(:,1)))==1);if isempty(endInd);endInd=startInd;end
                        if endInd(1)<startInd(1);startInd=[1;startInd];end %if starts right away
                        if length(endInd)>length(startInd);endInd(end)=[];end %if end too much
                        if length(startInd)>length(endInd);endInd=[endInd;length(approaches{u})];end %if end is missing
                        lengthSnips=endInd-startInd;
                        XYint=[];
                        
                        %snips have to be longer than 20
                        for l=1:length(lengthSnips)
                            if lengthSnips(l)>20
                                
                                try
                                    %put a tiny amount of noise on data pionts to
                                    %make them distinct
                                    temp=approaches{u}(startInd(l):endInd(l),:)+0.001*randn(size(approaches{u}(startInd(l):endInd(l),:)));
                                    intBin=3;
                                    XYint_snip=interparc(round(size(temp,1)/intBin),temp(~isnan(temp(:,1)),1),temp(~isnan(temp(:,1)),2),'pchip');
                                    %spline generates fake data, pchip does not
                                catch
                                    XYint_snip=[nan nan];
                                end
                                XYint=[XYint;XYint_snip;nan nan];%add a nan between all snippets to not create fake track
                            end
                        end
                        
                    end
                    
                    %calculate angle of tracks
                    intAngles=5;
                    if sum(~isnan(XYint))>2 %if the interpolation worked
                        intBin=3;
                        angleProbData=XYint(1:intAngles:end,:);%take only every 5th entry for general directions
                    else
                        intBin=3;
                        angleProbData=approaches{u}(1:intAngles*intBin:end,:);%take only every 5th entry for general directions
                    end
                    
                    [dist, alpha]=euclid_dist(angleProbData);%this uses atan2d method
                    %                     alpha=90-alpha;%to make the Y-axis the 0 axis
                    
                    if sum(~isnan(alpha))>3
                        %                         figure;subplot(1,2,1)
                        %                         plot(angleProbData(:,1),angleProbData(:,2));axis equal
                        %
                        %                         subplot(1,2,2)
                        %                         polarhistogram(deg2rad(alpha),12)
                        %                         set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
                        %
                        %calculate the mean vector length from the average
                        %angles
                        %do this for -90 to +90 as well, because 0 and 180 are the
                        %same, and if animal goes back and forwards, if
                        %goes "no vector" on average, even though has very
                        %strong pattern aligned vector -- double angles and
                        %half after
                        tempAlpha=alpha(~isnan(alpha));tempDist=dist(~isnan(alpha));
                        rmean(1,u)=circ_r(deg2rad(tempAlpha),tempDist/nanmax(tempDist));
                        dirmean(1,u)=rad2deg(circ_mean(deg2rad(tempAlpha),tempDist/nanmax(tempDist)));
                        %calculate the mean r and direction of a bimodal
                        %distribution, by doubling the angle data, and then
                        %halfing the dirmean afterwards again
                        rmean(2,u)=circ_r(deg2rad(tempAlpha*2),tempDist/nanmax(tempDist));
                        dirmean(2,u)=rad2deg(circ_mean(deg2rad(tempAlpha*2),tempDist/nanmax(tempDist)))/2;
                    else
                        rmean(:,u)=nan;
                        dirmean(:,u)=nan;
                    end
                    
                    %% extract body angle and proboscis contacts relative to pattern during probing
                    %calculate angle between head and thorax. do not need
                    %to smooth tracks, because the noise here doesn't
                    %affect overall angle
                    [~,angleBody]=euclid_dist(allHeadR(timeBouts(u,1):feedBouts(u,1),:),allThoraxR(timeBouts(u,1):feedBouts(u,1),:));
                    
                    %extract xposition of proboscis.
                    xProbPos=allprobR(timeBouts(u,1):feedBouts(u,1),1)-midpoint(:,1);
                    xHeadPos=allHeadR(timeBouts(u,1):feedBouts(u,1),1)-midpoint(:,1);
                    xThoraxPos=allThoraxR(timeBouts(u,1):feedBouts(u,1),1)-midpoint(:,1);
                    %turn into left and right relative to flower pattern
                    %(based on body orientation)
                    xProbPos(abs(angleBody)>90)=-xProbPos(abs(angleBody)>90);%turning points for angles are up and down, there is goes +-0 and +-180
                    xHeadPos(abs(angleBody)>90)=-xHeadPos(abs(angleBody)>90);%turning points for angles are up and down, there is goes +-0 and +-180
                    xThoraxPos(abs(angleBody)>90)=-xThoraxPos(abs(angleBody)>90);%turning points for angles are up and down, there is goes +-0 and +-180
                    indProbFirst=find(~isnan(xProbPos(:,1)),1,'first');%find first prob contact
                    indProbLast=indProbFirst+0.25*Fs;%find first prob contact
                    if length(xProbPos)<indProbLast; indProbLast=length(xProbPos);end
                    indAngleFirst=indProbFirst-0.25*Fs;%find first prob contact
                    if indAngleFirst<1; indAngleFirst=1;end
                    %turn body angle into left / right of pattern coodinates
                    angleBody(angleBody>0 & angleBody<=90)=-angleBody(angleBody>0 & angleBody<=90);
                    angleBody(angleBody<0 & angleBody>=-90)=-angleBody(angleBody<0 & angleBody>=-90);
                    angleBody(angleBody>90)=-(angleBody(angleBody>90)-180);
                    angleBody(angleBody<-90)=-(angleBody(angleBody<-90)+180);
                    
                    %save data
                    bodyAngle_all(u,:)=rad2deg([circ_mean(deg2rad(angleBody(indAngleFirst:indProbFirst))) circ_r(deg2rad(angleBody(indAngleFirst:indProbFirst)))...
                        circ_mean(deg2rad(angleBody(~isnan(angleBody)))) circ_r(deg2rad(angleBody(~isnan(angleBody))))]);
                    xBodyPos_all(u,:)=[nanmean(xProbPos(indProbFirst:indProbLast)) nanmean(xProbPos) ...
                        nanmean(xHeadPos(indProbFirst:indProbLast)) nanmean(xHeadPos) nanmean(xThoraxPos(indProbFirst:indProbLast)) nanmean(xThoraxPos)]/pixPerMM;
                    xProbPosHist(u,:)=histc(xProbPos/pixPerMM,xProbHistInd);%calculated in mm
                    
                    %% make the matrix of all contacts in approach
                    tempMatrix=nan(200,200,size(approaches{u},1));
                    for j=1:size(approaches{u},1)
                        if round(approaches{u}(j,2))+100<200 && round(approaches{u}(j,1))+100 < 200 %make this for safety so points outside flower do not crash programme
                            tempMatrix(round(approaches{u}(j,2))+100,round(approaches{u}(j,1))+100,j)=1;
                        else
                            tempMatrix(1,1,j)=nan;
                        end
                    end
                    %                 tempAMatrix(:,:,u)=nansum(tempMatrix,3);%collect data for each flight bout
                    tempAMatrix(:,:,u)=imgaussfilt(nansum(tempMatrix,3),gaussFiltKernel);
                    %calculate the contacts in 1mm around the cardinal and diagonal axes.
                    %for that, calculate how many pixels are 1mm
                    tempData=nansum(tempMatrix,3);tempData(tempData==0)=nan;
                    minInd=round(100-2*pixPerMM);maxInd=round(100+2*pixPerMM);
                   
                    
                    %% calculate contacts in and out of pattern
                    
                    %cardinal scores - around the line (+-2mm)
                    %in 0 and 90 orientation
                    scoresA(u,1)=nansum(nansum(tempData([1:minInd maxInd:200],minInd:maxInd)));%X
                    scoresA(u,2)=nansum(nansum(tempData(minInd:maxInd,[1:minInd maxInd:200])));%Y
                    %in 45 deg orientation (only for cross patterns)
                    tempDataR=imrotate(tempData,45,'crop');%appends additional rows and cols
                    scoresA(u,3)=nansum(nansum(tempDataR([1:minInd maxInd:200],minInd:maxInd)));%diag1
                    scoresA(u,4)=nansum(nansum(tempDataR(minInd:maxInd,[1:minInd maxInd:200])));%diag2
                    
                    %use interpolated data for contact scores, to get an
                    %idea where the proboscis was, but not how long
                    if sum(~isnan(XYint))>2 %if the interpolation worked
                        probPlacementData=XYint;%take only every 5th entry for general directions
                    else
                        probPlacementData=approaches{u}(1:intBin:end,:);%take only every 5th entry for general directions
                    end
                    plotTrackData{u}=probPlacementData;
                    % !!ATTENTION!! Here I use all prob placement data, not
                    % interpolated
                    probPlacementData=approaches{u}; %if I want how long, use all proboscis positions
                    
                    %calculate contacts at pattern, in central part and at
                    %rim edges of pattern of pattern -- with buffer space (rad_patternCount)
                    %make a polygon out of the line pattern (pattern width + buffer space)
                    bufferSpace=rad_patternCount*pixPerMM;
                    lengthPattern=lengthCircleFlower/2;%make the length a bit longer because sometimes pattern that was traced does not extend to edge of flower
                    wholeLine=polyshape([[-bufferSpace; -bufferSpace; bufferSpace;  bufferSpace], [-lengthPattern; lengthPattern; lengthPattern; -lengthPattern;]]);
                    %calculate contacts in polygon
                    scoresA(u,1)=nansum(inpolygon(probPlacementData(:,1),probPlacementData(:,2),wholeLine.Vertices(:,1),wholeLine.Vertices(:,2)));
                    wholeLine90 = rotate(wholeLine,90);
                    scoresA(u,2)=nansum(inpolygon(probPlacementData(:,1),probPlacementData(:,2),wholeLine90.Vertices(:,1),wholeLine90.Vertices(:,2)));
                    scoresA(u,5)=nansum(~isnan(probPlacementData(:,1)));%all prob contacts
                    %this line goes from rim (4mm) to other side rim - so
                    %whole line - 4mm at each ends
                    %                     yPoints=lengthCircleFlower/2-bufferSpace;
                    yPoints=lengthCircleFlower/4;%this is splitting inside and outside half-half
                    insideLine=polyshape([[-bufferSpace; -bufferSpace; bufferSpace;  bufferSpace], [-yPoints; yPoints; yPoints; -yPoints;]]);
                    %calculate contacts in polygon
                    scoresInside(u,1)=nansum(inpolygon(probPlacementData(:,1),probPlacementData(:,2),insideLine.Vertices(:,1),insideLine.Vertices(:,2)));
                    insideLine90 = rotate(insideLine,90);
                    scoresInside(u,2)=nansum(inpolygon(probPlacementData(:,1),probPlacementData(:,2),insideLine90.Vertices(:,1),insideLine90.Vertices(:,2)));
                    %calculate contacts in edges of pattern, made from
                    %subtraction of shapes
                    outsideLine=subtract(wholeLine,insideLine);
                    scoresOutside(u,1)=nansum(inpolygon(probPlacementData(:,1),probPlacementData(:,2),outsideLine.Vertices(:,1),outsideLine.Vertices(:,2)));
                    outsideLine90 = rotate(outsideLine,90);
                    scoresOutside(u,2)=nansum(inpolygon(probPlacementData(:,1),probPlacementData(:,2),outsideLine90.Vertices(:,1),outsideLine90.Vertices(:,2)));
                    
                else
                    approaches{u}=[nan nan];
                    appr_length(u)=nan;
                    appr_length(u)=nan;
                    appr_speed(u)=nan;
                    rmean(:,u)=nan;
                    dirmean(:,u)=nan;
                end
                
                %deleted calculating departures, because they don't exist
                %in this framework (only in bees)
                %(then departures = nan)
                
            else
                approaches{u}=[nan nan];
                appr_tort(u)=nan;
                appr_length(u)=nan;
                appr_speed(u)=nan;
            end
        end
        %calculate duration of approch and departure
        apprDur=(feedBouts(:,1)-timeBouts(:,1))/Fs;
        %         if length(apprDur)~=length(appr_tort)
        %             disp(i)
        %         end
        
        %calculate distance of proboscis and angle from center over time
        [distCentre, ~]=euclid_dist(allprobR,[0,0]);
        [~, anglCentre]=euclid_dist(allThoraxR,[0,0]);%use thorax for angle to centre
        probCenterDistAngl{i}(:,1)=distCentre./pixPerMM;
        probCenterDistAngl{i}(:,2)=anglCentre;
        
        %plot individual walking paths
        subplot(1,2,1);hold on;%approaches
        for u=1:size(timeBouts,1)
            if ~isempty(plotTrackData{u})
                %remove outside circle data
                in = inpolygon(plotTrackData{u}(:,1),plotTrackData{u}(:,2),xunitC,yunitC);
                plotTrackData{u}(in==0,:)=repmat([nan nan],sum(in==0),1);
                if ~isempty(plotTrackData{u})
                    %             plot(approaches{u}(:,1),approaches{u}(:,2),'-','color',ind_col(i,:))
                    plot(plotTrackData{u}(:,1), plotTrackData{u}(:,2),'-','color',ind_col(i,:))%interpolated tracks
                end
            end
        end
        
        
        %collect the summary data
        path_dur=[path_dur;apprDur i*ones(size(apprDur,1),1)];
        path_length=[path_length;appr_length' i*ones(size(appr_length',1),1)];
        path_tort=[path_tort;appr_tort' i*ones(size(appr_tort',1),1)];
        path_speed=[path_speed;appr_speed' i*ones(size(appr_speed',1),1)];
        path_r=[path_r;rmean' i*ones(size(rmean',1),1)];
        path_dir=[path_dir;dirmean' i*ones(size(dirmean',1),1)];
        allAMatrix=cat(3,allAMatrix,tempAMatrix);
        allScoresA=[allScoresA;scoresA i*ones(size(scoresA,1),1)];
        allScoresLine=[allScoresLine;scoresInside scoresOutside i*ones(size(scoresInside,1),1)];
        label=[label;repmat(filenames(i),size(apprDur,1),1)];
        animalID_all=[animalID_all;repmat(animalID(i),size(apprDur,1),1)];
        colourID_all=[colourID_all;repmat(colourInd{i},size(apprDur,1),1)];
        allBodyAngles=[allBodyAngles;bodyAngle_all];
        allxBodyPos=[allxBodyPos;	xBodyPos_all];
        allxProbHist=[allxProbHist;xProbPosHist];
        
        %collect data for PCA (average data per animal)
        PCA(i,1)=nansum(scoresA(:,1))/nansum(scoresA(:,5));%coverage (total number of contacts on pattern vs all contacts)
        PCA(i,2)=nansum(scoresInside(:,1))/nansum(scoresOutside(:,1));%uniformity (number of contacts on inside of pattern vs contacts outer rim )
        
        %sum all contacts per animal along Y-axis after normalisation
        Y_integrated = permute(nanmean(tempAMatrix,1),[2 3 1]);
        temp_Y=Y_integrated./repmat(nanmax(Y_integrated,[],1),size(Y_integrated,1),1);
        all_temp_Y=nanmean(temp_Y,2);all_temp_Y=all_temp_Y/nanmax(all_temp_Y);
        %find half width points
        index1 = find(all_temp_Y >= 0.5, 1, 'first');index2 = find(all_temp_Y >= 0.5, 1, 'last');
        PCA(i,3)=(index2-index1)/pixPerMM;%concentration of contacts along pattern axis (half width of contact distribution)
        PCA(i,4)=nanmean(nanmean(xProbPosHist,2));%lefty or righty (median of the proboscis histogram which is prob. position rel. to pattern viewed from body axis)
        PCA(i,5)=nanmean(apprDur);%duration of bout contact
        PCA(i,6)=nanmean(appr_length);%length of each bout
        PCA(i,7)=rad2deg(circ_mean(deg2rad(dirmean(2,~isnan(dirmean(2,:)))')));%mean direction of probing
        PCA(i,8)=nanmean(path_r(:,2));%mean vector length of probing
    end
end

%% individual tracks
%plot the circles over the tracks
% plot heatmap all contacts
patternFinal0=checkCircle(patternNorm);%get patterns into right dimensions with +100
pattern2Final0=checkCircle(pattern2Norm);

subplot(1,2,1);hold on;
title('approaches')
plot(xunitC,yunitC,'-k');
plot(patternFinal0(:,1),patternFinal0(:,2),'k')
plot(pattern2Final0(:,1),pattern2Final0(:,2),'k')
xlim([-100 100]);ylim([-100 100])
axis square

%% plot the summary stats for duration and tortuosity
f2=figure('Position',[500  500  850  400]);
labeltext={'all contacts'};

subplot(1,3,1)
indCols=[1];
for u=1:length(indCols)
    v=violin(path_dur(:,indCols(u)),'x',u,'facecolor',favCol,'medc','k');
    n=numel(path_dur(:,indCols(u)));
    plot(spread*rand(n,1)+u*ones(n,1)-0.5*spread,path_dur(:,indCols(u)),'.','MarkerSize',12,'color','k');hold on
end
b=boxplot(path_dur(:,indCols),'color','w','PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol','','labels',labeltext,'LabelOrientation','horizontal');
ylabel('bout contact time (s)');
ylim([0 5])

subplot(1,3,2)
indCols=[1];
for u=1:length(indCols)
    v=violin(path_length(:,indCols(u)),'x',u,'facecolor',favCol,'medc','k');
    n=numel(path_length(:,indCols(u)));
    plot(spread*rand(n,1)+u*ones(n,1)-0.5*spread,path_length(:,indCols(u)),'.','MarkerSize',12,'color','k');hold on
end
b=boxplot(path_length(:,indCols),'color','w','PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol','','labels',labeltext,'LabelOrientation','horizontal');
ylabel('path length (mm)');
ylim([0 700])

% subplot(1,3,3)
% indCols=[1];
% for u=1:length(indCols)
%     v=violin(path_tort(:,indCols(u)),'x',u,'facecolor',favCol,'medc','k');
%     n=numel(path_tort(:,indCols(u)));
%     plot(spread*rand(n,1)+u*ones(n,1)-0.5*spread,path_tort(:,indCols(u)),'.','MarkerSize',12,'color','k');hold on
% end
% b=boxplot(path_tort(:,indCols),'color','w','PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol','','labels',labeltext,'LabelOrientation','horizontal');
% set(gca,'yscale','log')
% ylabel('tortuosity');
% axy=gca;
% ylim([1 100])


%% plot heatmaps and analysis

% plot heatmap all contacts
patternFinal=checkCircle(patternNorm+100);%get patterns into right dimensions with +100
pattern2Final=checkCircle(pattern2Norm+100);
cmap=colormap('inferno');

f3=figure('Position',[500  500  850  400]);hold on;
title('approaches')
sumMatrixData=nansum(allAMatrix,3); %this is summing or averaging over all
%     imagesc(sumMatrixData,[0 quantile(sumMatrixData(:),.99)]);
imagesc(sumMatrixData);
plot(patternFinal(:,1),patternFinal(:,2),'w')
plot(pattern2Final(:,1),pattern2Final(:,2),'w')
plot(xunitC+100,yunitC+100,'w')
xlim([0 200]);ylim([0 200])
ylabel('pixels');xlabel('pixels');
axis square
colorbar


%% plot scores
%for analysis of CROSS patterns, both arms should be the same, therefore
%add together
if size(unique(patternR(:,1)),1)>4 || sum(~isnan(patternR(:)))==0 %this if the case for the cross pattern or no pattern
    finalScoresA=[allScoresA(:,1)+allScoresA(:,2) allScoresA(:,3)+allScoresA(:,4)];
    scoreLabel={'cardinal','diagonal'};
else
    %for analysis of single LINE patterns, there is only a vertical arm, and it
    %should be compared to the horizontal probes, not the diagonal
    finalScoresA=[allScoresA(:,1) allScoresA(:,2)];
    scoreLabel={'vertical','horizontal'};
end

%exclude all scores that do not match condition (i.e. duration on the
%flower
% finalScoresD(path_dur(:,2)<1,:)=nan;

f4=figure('Position',[500 200 800 500]);
subplot(1,3,1)
% analyseScore=finalScoresA;%this is all data
analyseScore=grpstats(allScoresA(:,[1 2]),allScoresA(:,5),'sum'); %this is pooled data

if size(analyseScore,1)>1
    boxplot(analyseScore);hold on;
    for u=1:size(analyseScore,1)
        plot([1+0.01*randn(1,1),2+0.01*randn(1,1)],[analyseScore(u,1),analyseScore(u,2)],'k.-','Markersize',10);
    end
    set(gca,'xticklabel',scoreLabel)
else
    plot([1 2],finalScoresA,'ro');xlim([0 3]);
    set(gca,'xtick',[1 2],'xticklabel',scoreLabel)
end
[p,h,stats]=signrank(analyseScore(:,1),analyseScore(:,2));
title(sprintf('whole line n=%2.0f, p=%1.3f, rank=%1.3u',sum(isnan(analyseScore(:,1))==0),p,stats.signedrank));
ylabel('summed contacts');

subplot(1,3,2)
% analyseScore=allScoresLine(:,[1 2]);%this is all data
analyseScore=grpstats(allScoresLine(:,[1 2]),allScoresLine(:,5),'sum') %this is pooled data

if size(analyseScore,1)>1
    boxplot(analyseScore(:,[1 2]));hold on;
    for u=1:size(analyseScore,1)
        plot([1+0.01*randn(1,1),2+0.01*randn(1,1)],[analyseScore(u,1),analyseScore(u,2)],'k.-','Markersize',10);
    end
    set(gca,'xticklabel',scoreLabel)
else
    plot([1 2],analyseScore(:,[1 2]),'ro');xlim([0 3]);
    set(gca,'xtick',[1 2],'xticklabel',scoreLabel)
end
[p,h,stats]=signrank(analyseScore(:,1),analyseScore(:,2));
title(sprintf('inside line n=%2.0f, p=%1.3f, rank=%1.3u',sum(isnan(analyseScore(:,1))==0),p,stats.signedrank));
ylabel('summed contacts');

analyseScoreIn=analyseScore;

subplot(1,3,3)
% analyseScore=allScoresLine(:,[3 4]);%this is all data
analyseScore=grpstats(allScoresLine(:,[3 4]),allScoresLine(:,5),'sum') %this is pooled data

if size(allScoresLine,1)>1
    boxplot(analyseScore(:,[1 2]));hold on;
    for u=1:size(analyseScore,1)
        plot([1+0.01*randn(1,1),2+0.01*randn(1,1)],[analyseScore(u,1),analyseScore(u,2)],'k.-','Markersize',10);
    end
    set(gca,'xticklabel',scoreLabel)
else
    plot([1 2],analyseScore(:,[1 2]),'ro');xlim([0 3]);
    set(gca,'xtick',[1 2],'xticklabel',scoreLabel)
end
[p,h,stats]=signrank(analyseScore(:,1),analyseScore(:,2));
title(sprintf('outside line n=%2.0f, p=%1.3f, rank=%1.3u',sum(isnan(analyseScore(:,1))==0),p,stats.signedrank));
ylabel('summed contacts');

currPath=pwd;inds=strfind(currPath,'\');
line_scores_tbl=table(analyseScoreIn(:,1),analyseScoreIn(:,2),analyseScore(:,1),analyseScore(:,2),unique(allScoresLine(:,5)),'VariableNames',{'Y_in','X_in','Y_out','X_out','animal'})
writetable(line_scores_tbl,[currPath(inds(end)+1:end),'_line_scores.csv'])

%% plot direction analysis of probing tracks
f5=figure('Position',[500  500  850  400]);
labeltext={'all contacts'};

subplot(1,3,1)
indCols=[1]; %circ averages for single mode 1, axial distribution 2
for u=1:length(indCols)
    v=violin(path_r(:,indCols),'x',u,'facecolor',favCol,'medc','k');
    n=numel(path_r(:,indCols));
    plot(spread*rand(n,1)+u*ones(n,1)-0.5*spread,path_r(:,indCols),'.','MarkerSize',12,'color','k');hold on
end
b=boxplot(path_r(:,indCols),'color','w','PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol','','labels',labeltext,'LabelOrientation','horizontal');
ylabel('vector strength');
% ylim([0 1000])

subplot(1,3,2)
binInt=12;%degrees, circular histogram bin intervals - for Alex
% binInt=10;%degrees, circular histogram bin intervals - for Alex
circBins=(0:binInt:360)-binInt/2;%circular histogram bins - centered on 0
binCenters=round(circBins(1:end-1)+binInt/2);

polarhistogram(deg2rad(path_dir(:,indCols)),deg2rad(circBins))
set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
title('path directions')

subplot(1,3,3);
for n=1:size(path_dir,1)
    polarplot([0; deg2rad(path_dir(n,indCols))],[0; path_r(n,indCols)],'k');hold on;
end
set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
title('path vectors')

%save the orientation of tracks data
%calculate the number of entries in sectors for statistical analysis
    if size(unique(patternR(:,1)),1)>4 || sum(~isnan(patternR(:)))==0 %this if the case for the cross pattern or no pattern
        error('sectors not defined yet')
      
        sectorRatio=.5;
    else
        sector_scores1=path_dir(:,1)>150 | path_dir(:,1)<-150;
        sector_scores2=path_dir(:,1)>0 & path_dir(:,1)<30;
        sector_scores3=path_dir(:,1)<0 & path_dir(:,1)>-30;
        sector_scores=(sector_scores1+sector_scores2+sector_scores3);
        sectorRatio=.33;
    end

    % for statistical analysis, consider only proboscis track vectors with
    % longer than 0.2 vector strength
    for n=1:length(animalID)
    counts(n,:)=[nansum(sector_scores(animalID_all==animalID(n) & path_r(:,1)>0.2,1)) ...
                nansum(animalID_all(~isnan(path_dir(:,1)))==animalID(n))];
    end

    tb21 = table(counts(:,1),counts(:,2)-counts(:,1),counts(:,2),repmat(sectorRatio,size(counts,1),1),animalID,'VariableNames',{'patternCount','bgCount','allCounts','threshold','animalID'});
    
    currPath=pwd;inds=strfind(currPath,'\');
    writetable(tb21,[currPath(inds(end)+1:end),'_probTrackCounts.csv'])

    %save path direction and straightness
    tb212 = table(path_dir(:,1),path_r(:,1),animalID_all,pattern,'VariableNames',{'vectorDir','vectorR','animalID'});
    currPath=pwd;inds=strfind(currPath,'\');
    writetable(tb212,[currPath(inds(end)+1:end),'_trackVectors.csv'])

%% plot body orientation and probing position relative to stripe
f6=figure('Position',[200  200  1000  600]);
labeltext={'proboscis','head', 'thorax'};

subplot(2,3,1)
indCols=[1,3,5];%these are the first 250ms position data
for u=1:length(indCols)
    v=violin(allxBodyPos(:,indCols(u)),'x',u,'facecolor',favCol,'medc','k');
    n=numel(allxBodyPos(:,indCols(u)));
    plot(spread*rand(n,1)+u*ones(n,1)-0.5*spread,allxBodyPos(:,indCols(u)),'.','MarkerSize',12,'color','k');hold on
end
b=boxplot(allxBodyPos(:,indCols),'color','w','PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol','','labels',labeltext,'LabelOrientation','horizontal');
ylabel('x-position (mm)');
title('first 250ms')

subplot(2,3,2)
indCols=[2,4,6];%these are the average position data
for u=1:length(indCols)
    v=violin(allxBodyPos(:,indCols(u)),'x',u,'facecolor',favCol,'medc','k');
    n=numel(allxBodyPos(:,indCols(u)));
    plot(spread*rand(n,1)+u*ones(n,1)-0.5*spread,allxBodyPos(:,indCols(u)),'.','MarkerSize',12,'color','k');hold on
end
b=boxplot(allxBodyPos(:,indCols),'color','w','PlotStyle','compact','whisker',1.5,'jitter',0.1,'Symbol','','labels',labeltext,'LabelOrientation','horizontal');
ylabel('x-position (mm)');
title('bout averages')
% ylim([0 1000])

subplot(2,3,3)
relProbPosHist=allxProbHist./repmat(nanmax(nanmean(allxProbHist,1)),size(allxProbHist,1),1);
relProbPosHist=relProbPosHist./repmat(nanmax(nanmean(relProbPosHist,1)),size(relProbPosHist,1),1);
plotProbHist=relProbPosHist;%plotProbHist=allxProbHist;
s1=shadedErrorBar_anna(xProbHistInd,nanmean(plotProbHist,1),nanstd(plotProbHist,1,1)/sqrt(size(plotProbHist,1)-1))
hold on;plot([0 0],[0 nanmax(s1.patch.Vertices(:))],'k--')
% ylim([0 nanmax(s1.patch.Vertices(:,2))])
ylim([0 1.1])
xlim([-20 20])
xlabel('prob position rel. pattern (mm)')
ylabel('number of occurances')



subplot(2,3,4);
for n=1:size(path_dir,1)
    polarplot([0; deg2rad(allBodyAngles(n,1))],[0; allBodyAngles(n,2)],'k');hold on;
end
set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
title('first 250ms body vectors')

subplot(2,3,5);
for n=1:size(path_dir,1)
    polarplot([0; deg2rad(allBodyAngles(n,3))],[0; allBodyAngles(n,4)],'k');hold on;
end
set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
title('average body vectors')

%save median body position relative to pattern data
%this is for all trials, positions only during probing in trial
body_pos_midline_tbl=table(allxBodyPos(:,2),allxBodyPos(:,4),allxBodyPos(:,6),'VariableNames',{'med_prob_pos_mm','med_head_pos_mm','med_thorax_pos_mm'})
writetable(body_pos_midline_tbl,'allRotData_med_body_pos_pattern.csv')

%generate value for each animal, weighted by length of trial
allAnimalBodyPos=nan(length(animalID),3);

%hack different patterns into diffreent numbers to separate animal IDs
%where same animal indicator, but different set of experiments
% tempPatterns(find(contains(patterns,temp{y})))=scales(y);

%use colourID, which is taken from the pattern variable, and a different
%number for all patterns

tempAnimalID=animalID.*cell2mat(colourInd);

for u=1:length(tempAnimalID)
    inds = find(animalID_all.*colourID_all==tempAnimalID(u));

% allAnimalBodyPos(u,:)=nansum(allxBodyPos(inds,[2 4 6]).*path_dur(inds,1))/nansum(path_dur(inds,1)); 
allAnimalBodyPos(u,:)=nanmedian(allxBodyPos(inds,[1 3 5])); %not weighted

end
body_pos_midline_anim_tbl=table(allAnimalBodyPos(:,1),allAnimalBodyPos(:,2),allAnimalBodyPos(:,3),'VariableNames',{'med_prob_pos_mm','med_head_pos_mm','med_thorax_pos_mm'})
writetable(body_pos_midline_anim_tbl,'allAnimals_med_body_pos_pattern.csv')

%% plot distance and angle to centre over time
% nrAnimalss=length(probCenterDistAngl);
% 
% figure('position',[50 50 1200 800]);hold on;
% for i=1:nrAnimalss
%     if ~isempty(probCenterDistAngl{i})
%         ax(1)=subplot(2,1,1);hold on;
%         plot(1/Fs:1/Fs:length(probCenterDistAngl{i}(:,1))/Fs,probCenterDistAngl{i}(:,1),'-','color',ind_col(i,:));xlabel('time (s)');ylabel('prob. distance to centre (mm)')
%         ax(2)=subplot(2,1,2);hold on;
%         plot(1/Fs:1/Fs:length(probCenterDistAngl{i}(:,2))/Fs,probCenterDistAngl{i}(:,2),'.-.','color',ind_col(i,:));xlabel('time (s)');ylabel('thorax angle to centre (deg)')
%     end
%     linkaxes(ax,'x')
%     legend(num2str(animalID(~cellfun(@isempty,probCenterDistAngl))))
% end
% 
% figure('position',[50 50 1200 nrAnimalss*300]);hold on;
% for i=1:nrAnimalss
%     if ~isempty(probCenterDistAngl{i})
%         ax(1)=subplot(nrAnimalss,2,(i*2)-1);hold on;
%         plot(1/Fs:1/Fs:length(probCenterDistAngl{i}(:,1))/Fs,probCenterDistAngl{i}(:,1),'-','color',ind_col(i,:));set(gca,'XTick',[]);
%         ax(2)=subplot(nrAnimalss,2,(i*2));hold on;
%         plot(1/Fs:1/Fs:length(probCenterDistAngl{i}(:,2))/Fs,probCenterDistAngl{i}(:,2),'.-.','color',ind_col(i,:));set(gca,'XTick',[]);
%         linkaxes(ax,'x')
%         legend(num2str(animalID(i)))
%     end
% end
% xticks(ax(1),'auto');xticks(ax(2),'auto');
% set(get(ax(1),'Xlabel'),'String','time (s)');set(get(ax(2),'Xlabel'),'String','time (s)');
% set(get(ax(1),'Ylabel'),'String','prob. distance to centre (mm)');
% set(get(ax(2),'Ylabel'),'String','prob. distance to centre (mm)');

%% data summary

% save('contact_data.mat','path_dur','path_length','path_tort','path_r','path_dir','finalScoresA','animalID_all','label')

% patternName=filenames{i}(strfind(filenames{i},'M'):find(filenames{i}(strfind(filenames{i},'M'):end)=='_',1,'first')+strfind(filenames{i},'M')-2);
% save(['heatmap_data_',patternName,'.mat'],'animal','finalScores','patterns','chiXp');

print(f1,'probTracks.jpg','-djpeg','-r300')
print(f2,'pathData.jpg','-djpeg','-r300')
print(f3,'Heatmap.jpg','-djpeg','-r300')
print(f4,'patternScores.jpg','-djpeg','-r300')
print(f5,'orientationTracks.jpg','-djpeg','-r300')
print(f6,'probPos.jpg','-djpeg','-r300')

print(f1,'probTracks.eps','-dpdf','-r300','-painters','-bestfit')
print(f2,'Approach_Depart_pathData.eps','-dpdf','-r300','-painters','-bestfit')
print(f3,'walkingHeatmap.eps','-dpdf','-r300','-painters','-bestfit')
print(f4,'patternScores.eps','-dpdf','-r300','-painters','-bestfit')
print(f5,'orientationTracks.eps','-dpdf','-r300','-painters','-bestfit')
print(f6,'probPos.eps','-dpdf','-r300','-painters','-bestfit')

%% save walking paths contact counts along patterns for stat. analysis
currPath=pwd;inds=strfind(currPath,'\');

%save data for PCA to compare different contact pattern types
% PCAtable=table(PCA(:,1),PCA(:,2),PCA(:,3),PCA(:,4),PCA(:,5),PCA(:,6),PCA(:,7),PCA(:,8),'VariableNames',{'coverage','uniformity','concentration','left_right','duration','length','direction','vector_length'})
% writetable(PCAtable,[currPath(inds(end)+1:end),'_PCA.csv'])
% save([currPath(inds(end)+1:end),'_PCA.mat'],PCA);

path_distCentr

%save contact length and contact time per animal
dur_length_tbl=table(animalID,grpstats(path_dur(:,1),path_dur(:,2),'mean'),grpstats(path_length(:,1),path_length(:,2),'mean'),grpstats(path_tort(:,1),path_tort(:,2),'mean'),grpstats(path_speed(:,1),path_speed(:,2),'mean'),'VariableNames',{'animal','duration_s','length_mm','tortuosity','speed_mms'})
writetable(dur_length_tbl,[currPath(inds(end)+1:end),'_dur_length_tort_speed.csv'])
% dur_length_data=dur_length_tbl{:,:};save([currPath(inds(end)+1:end),'_dur_length.mat'],'dur_length_data');



end

function patternFinal=checkCircle(patternNorm)
%check if pattern is a circle -- which has 1 entry, where the distance in x
%and y is large, while for oriented lines and crosses, always one of the
%dimensions is close to 0
if sum(~isnan(patternNorm))>=2
    diffP(:,1)=diff(patternNorm(:,1));diffP(:,2)=diff(patternNorm(:,2));
    if size(patternNorm,1)==4 && nansum(mink(abs(diffP(:)),3))>6
        temp=diff(patternNorm);temp(2,:)=[];
        midpoint=[max(abs(temp(:,1)))/2+min(patternNorm(:,1))  max(abs(temp(:,2)))/2+min(patternNorm(:,2))];
        [D, alpha]=euclid_dist(patternNorm);
        lengthCircle=nanmean(D([1,3],:));
        th = 0:pi/50:2*pi;
        xunitC = lengthCircle/2 * cos(th) + midpoint(1);
        yunitC = lengthCircle/2 * sin(th) + midpoint(2);
        patternFinal=[xunitC' yunitC'];
    else
        patternFinal=[patternNorm;patternNorm(1,:)];%close the pattern
    end
else
    patternFinal=patternNorm;
end
end
%%
