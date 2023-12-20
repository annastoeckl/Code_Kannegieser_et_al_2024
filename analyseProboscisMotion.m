%analyse the movement of the proboscis and the body in relation to each other, 
% and to the pattern. It generates Fast Fourier Transforms of body and proboscis movements, 
% relative to the pattern axes, and relative to the animals' body axis.

clear all
close all

%some analysis was performed for patterns were animal rotation was resicted
%to aligned with the pattern (0) or perpendicular to it (90). 
%select rotation group (0,90,nan is no restriction)
patternRot=nan;

%load data
[filenames, path]=uigetfile('allRot*.mat','multiselect','on');
%if only one selected turn filenames into cell
if iscell(filenames)==0
    temp=filenames;clear filenames;
    filenames=cell(1,1);filenames{1}=temp;
end

%check if filenames are on list of good animals (those that touched the
%pattern)
if exist('goodAnimals.txt')
    Tbl=readtable('goodAnimals.txt');
    if ~isempty(Tbl)
        goodLabels=Tbl{:,1};

        for i=1:length(filenames)
            for j=1:length(goodLabels)
                temp(j)=~isempty(strfind(filenames{i},goodLabels{j}));
            end
            indgood(i)=nansum(temp);
        end
        filenames(indgood==0)=[];
    end
end

Fs = 200;            % Sampling frequency
binEdges=[0:0.25:20];
margin=90;%margin around 0 and 90 degrees 
%this is the minimum length of proboscis touch bout required for FFT analysis
boutLength=10; 
fmin=0.5;fmax=75;%min and max frequency for fourier analysis -
fmax=100;
fNew=fmin:0.2:fmax;
fSumMin=4;fSumMax=15;%min and max frequency for sum analysis
minLengthFFT=1*Fs;%minimum length of data we need for FFT to give sensible results

%initialise parameters
turnVelocity=nan(length(filenames),3);%median turning velocity and lower and upper quartile
% prefTurnDir=nan(length(filenames),4);%preferred sign of turning, and signtest against distribution with 0 median (p zval sign)
medM_probNorm0=nan(length(filenames),2);
medM_probNorm90=nan(length(filenames),2);
histM_probNorm0=nan(length(filenames),length(binEdges),2);
histM_probNorm90=nan(length(filenames),length(binEdges),2);
specM_probNorm=nan(length(filenames),length(fNew),2);
specM_prob=nan(length(filenames),length(fNew),2);
specM_head=nan(length(filenames),length(fNew),2);
specM_thorax=nan(length(filenames),length(fNew),2);
specM_probNorm0=nan(length(filenames),length(fNew),2);
specM_probNorm90=nan(length(filenames),length(fNew),2);
specM_probDist=nan(length(filenames),length(fNew),3);
specM_probDist0=nan(length(filenames),length(fNew),3);
specM_probDist90=nan(length(filenames),length(fNew),3);
allIsolProbPos=nan(120,120,length(filenames));
allOrigProbPos=nan(120,120,length(filenames));
allIsolProbPos0=nan(120,120,length(filenames));
allIsolProbPos90=nan(120,120,length(filenames));

specM_prob0=nan(length(filenames),length(fNew),2);
specM_prob90=nan(length(filenames),length(fNew),2);
specM_head0=nan(length(filenames),length(fNew),2);
specM_head90=nan(length(filenames),length(fNew),2);
specM_thorax0=nan(length(filenames),length(fNew),2);
specM_thorax90=nan(length(filenames),length(fNew),2);

absMov_prob0=cell(length(filenames),1);
absMov_head0=cell(length(filenames),1);
absMov_thorax0=cell(length(filenames),1);
absMov_prob90=cell(length(filenames),1);
absMov_head90=cell(length(filenames),1);
absMov_thorax90=cell(length(filenames),1);

specM_probAbs=nan(length(filenames),length(fNew),2);
specM_thoraxAbs=nan(length(filenames),length(fNew),2);

meanAngles=nan(length(filenames),1);
meanBodyDist=nan(length(filenames),3);
meanAngles0=nan(length(filenames),1);
meanAngles90=nan(length(filenames),1);
meanBodyDist0=nan(length(filenames),3);
meanBodyDist90=nan(length(filenames),3);

xpos_bout_all=[];

for i=1:length(filenames)
    cd(path)
    load(filenames{i})
    disp(filenames{i})
    %extract file name items
    indT=strfind(filenames{i},'_T');
    %animal name
    animal(i)=str2num(filenames{i}(indT+2:find(filenames{i}(indT+2:end)=='_',1,'first')+indT));
    %pattern type
    startInd=strfind(filenames{i},'_M');
    if ~isempty(strfind(filenames{i},'MC'))%for circle patterns
        startInd=startInd+1;
    end
    endInd=strfind(filenames{i}(startInd+2:end),'_');
    if ~isempty(filenames{i}(startInd+2:startInd+endInd(1))) && ~isempty(str2num(filenames{i}(startInd+2:startInd+endInd(1))))
        patterns(i)=str2num(filenames{i}(startInd+2:startInd+endInd(1)));
    else
        patterns(i)=nan;
    end
    %find the date as unique identifier
    startInd2=strfind(filenames{i},'_20');
    %check that the  entry before next _ is exactly 8 digits long
    if length(startInd2)>=1
        endChar=filenames{i}(startInd2+9);
        startInd2_final=startInd2(endChar=='_');
    end

    if ~isempty(startInd2) && ~isempty(str2num(filenames{i}(startInd2+1:startInd2+8)))
        date(i)=str2num(filenames{i}(startInd2+1:startInd2+8));
    else
        date(i)=nan;
    end

    %since this script is about checking how the movements of the body and
    %proboscis influence pattern touching, movement of the body outside of
    %proboscis movement is not relevant. delete this
    allHeadR(isnan(allprobR))=nan;
    allThoraxR(isnan(allprobR))=nan;
    % calculate pixels per mm
    %the outlineR, which is generally parallel to the coordinate system to start with, might become turned... so need to
    %take the euclidian distance between points rather than the actual distance
    temp1=sqrt((outlineR(1,1)-outlineR(2,1))^2+(outlineR(1,2)-outlineR(2,2))^2);
    temp2=sqrt((outlineR(3,1)-outlineR(4,1))^2+(outlineR(3,2)-outlineR(4,2))^2);
    pixPerMM(i)=nanmean([temp1,temp2])/38;%38 mm flowers

    %use only files of minimum length
    if size(allprobR,1)>500 && sum(isnan(allprobR(:,1))==0)>minLengthFFT

        %% subtract body movement from proboscis movement
        if exist('allHeadR') && exist('allThoraxR')
            %first turn the frame so that the head is always up

            %                         figure;hold on; %control figure to check that all were turned
            %             to have head up
            alphas=nan(size(allHeadR,1),1);
            flipInd=ones(size(allHeadR,1),1);%indicates if animal's position is flipped (head down)
            tempProbR=nan(size(allHeadR));
            tempHeadR=nan(size(allHeadR));
            tempThoraxR=nan(size(allHeadR));

            for j=1:size(allHeadR,1)

                if isnan(allprobR(j,1))==0
                    %turn all head-body axes parallel to y axis
                    %extract angle of body axis (head-thorax) in overall
                    %coordinate system (of camera)
                    alpha=atand((allHeadR(j,2)-allThoraxR(j,2))/(allHeadR(j,1)-allThoraxR(j,1)));
                    %NOTE: that means movement in the direction of the
                    %pattern is seen in the second coordinate,
                    %and movement perpendicular to it in the first
                    %coordinate - axis

                    %turn all head-body axes parallel to y axis
                    alphas(j)=-(-90-alpha);

                    %this is to turn the animals parallel to pattern,
                    %because with the tangens, the angles that come out are
                    %90/-90 at the pattern axis and 0 perpendicular
                    alpha=-(-90-alpha);

                    R=-[cosd(-alpha) -sind(-alpha); sind(-alpha) cosd(-alpha)];
                    %                     R=1;

                    tempProbR(j,:)=[R*allprobR(j,:)']';
                    tempHeadR(j,:)=[R*allHeadR(j,:)']';
                    tempThoraxR(j,:)=[R*allThoraxR(j,:)']';

                    %                     make sure head is always up
                    if tempHeadR(j,2) < tempThoraxR(j,2)
                        R=[cosd(-180) -sind(-180); sind(-180) cosd(-180)];
                        tempProbR(j,:)=[R*tempProbR(j,:)']';
                        tempHeadR(j,:)=[R*tempHeadR(j,:)']';
                        tempThoraxR(j,:)=[R*tempThoraxR(j,:)']';
                        flipInd(j)=-1;%register the turn in the rotation value
                    end
                    %                                         plot([allprobR(j,1);allHeadR(j,1);allThoraxR(j,1)],[allprobR(j,2);allHeadR(j,2);allThoraxR(j,2)],'b-o');
                    %                                         hold on; plot(allHeadR(j,1),allHeadR(j,2),'co');plot(allprobR(j,1),allprobR(j,2),'go');

                end
            end

            %then subtract head position
            tempHeadR(isnan(tempProbR(:,1)),:)=nan;
            probNormR=tempProbR-tempHeadR;
            %when the whole animal is rotated, there can be jumps in the
            %traces (from +x or y to - or vice versa), these should be eliminated
            %as they don't constitute real movement
            probNormR(abs(diff(probNormR(:,2)))>50,:)=nan;
            probNormR(abs(diff(probNormR(:,1)))>50,:)=nan;

            %recreate isolated proboscis tracks for individual time bouts
            %and calculate contact distribution over x axis
            %this is to see if individual animals have individual
            %preferences 
            xpos_bout=nan;
            %             boutInds=-60:0.5:60;
            %             xpos_bout=nan(size(timeBouts,1),size(boutInds,2));
            %             for u=1:size(timeBouts,1)
            %             probNormR_bout=probNormR(timeBouts(u,1):timeBouts(u,2),:);
            %             xpos_bout(u,:)=histc(probNormR_bout(:,1),boutInds);
            %             end

            %sort isolated proboscis data into matrix to plot proboscis
            %positions across hawkmoths - matrix goes from -50 to +50
            %*convert to 0 to 100, and from 0 to 100 in y axis. head is at
            %0/0
            tempMatrix=nan(120,120,size(probNormR,1));
            for j=1:size(probNormR,1)
                if round(probNormR(j,2))<120 && round(probNormR(j,1))+60 < 120 && round(probNormR(j,2))>0 && round(probNormR(j,1))+60>0%make this for safety so points outside flower do not crash programme
                    tempMatrix(round(probNormR(j,2)),round(probNormR(j,1))+60,j)=1;
                else
                    tempMatrix(1,1,j)=nan;
                end
            end
            allIsolProbPos(:,:,i)=nansum(tempMatrix,3);

            %sort original proboscis data (relative to pattern coordinate system) into matrix to plot proboscis
            %positions across hawkmoths - matrix goes from -50 to +50
            %*convert to 0 to 100, and from 0 to 100 in y axis. head is at
            %0/0
            tempMatrix2=nan(120,120,size(allprobR,1));
            for j=1:size(allprobR,1)
                if round(allprobR(j,2))<120 && round(allprobR(j,1))+60 < 120 && round(allprobR(j,2))>0 && round(allprobR(j,1))+60>0%make this for safety so points outside flower do not crash programme
                    tempMatrix(round(allprobR(j,2)),round(allprobR(j,1))+60,j)=1;
                else
                    tempMatrix2(1,1,j)=nan;
                end
            end
            allOrigProbPos(:,:,i)=nansum(tempMatrix2,3);
            %check that turning worked
            %                 figure;subplot(1,2,1);hold on;
            %                 plot(allHeadR(:,1)-allHeadR(:,1),allHeadR(:,2)-allHeadR(:,2),'ro');
            %                 plot(allThoraxR(:,1)-allHeadR(:,1),allThoraxR(:,2)-allHeadR(:,2),'bo');
            %                 subplot(1,2,2);hold on;
            %                 plot(tempHeadR(:,1)-tempHeadR(:,1),tempHeadR(:,2)-tempHeadR(:,2),'ro');
            %                 plot(tempThoraxR(:,1)-tempHeadR(:,1),tempThoraxR(:,2)-tempHeadR(:,2),'bo');

            %calculate body angle and prob-head distance
            %extract angle of body axis (head-thorax) in overall
            %coordinate system (of camera)
            alpha1=atand((allHeadR(:,2)-allThoraxR(:,2))./(allHeadR(:,1)-allThoraxR(:,1)));

            %exact relative angle between head-thorax axis and
            %proboscis tip
            alpha2=atand((allprobR(:,2)-allHeadR(:,2))./(allprobR(:,1)-allHeadR(:,1)));
            probAngle=alpha2-alpha1;
            %recitfy the 180 flips
            probAngle(probAngle>100)=probAngle(probAngle>100)-180;
            probAngle(probAngle<-100)=probAngle(probAngle<-100)+180;

            %extract absolute distance between proboscis and head,
            %and head and thorax
            %added proboscis and thorax, to see if thorax movement jitter
            %from flight is also visible here
            bodyPosDist=[euclid_dist(allprobR,allHeadR),euclid_dist(allHeadR,allThoraxR),euclid_dist(allprobR,allThoraxR)];
        end

        %% save proboscis angle and body position distance (average turning speed (and direction))

        meanAngles(i)=nanstd(probAngle);
        meanBodyDist(i,:)=nanmean(bodyPosDist,1);

        %unwrap heading angles to get smooth curve
        unwrap_alphas=alphas;
        for m=2:length(unwrap_alphas)
            if (unwrap_alphas(m)-unwrap_alphas(m-1))>150
                unwrap_alphas(m:end)=unwrap_alphas(m:end)-180;
            end
            if (unwrap_alphas(m)-unwrap_alphas(m-1))<-150
                unwrap_alphas(m:end)=unwrap_alphas(m:end)+180;
            end
        end


        %% analyse movement

       
        probNormR0=probNormR;probNormR0(alphas>=margin/2 & alphas<=180-margin/2,:)=nan;
        probNormR90=probNormR;probNormR90(alphas>=90+margin/2 | alphas<=90-margin/2,:)=nan;

        %entire proboscis movement in 0 and 90 orientation
        probR0=allprobR;probR0(alphas>=margin/2 & alphas<=180-margin/2,:)=nan;
        probR90=allprobR;probR90(alphas>=90+margin/2 | alphas<=90-margin/2,:)=nan;
        %head movement in 0 and 90 orientation
        headR0=allHeadR;headR0(alphas>=margin/2 & alphas<=180-margin/2,:)=nan;
        headR90=allHeadR;headR90(alphas>=90+margin/2 | alphas<=90-margin/2,:)=nan;

        %thorax movement in 0 and 90 orientation
        thoraxR0=allThoraxR;thoraxR0(alphas>=margin/2 & alphas<=180-margin/2,:)=nan;
        thoraxR90=allThoraxR;thoraxR90(alphas>=90+margin/2 | alphas<=90-margin/2,:)=nan;

        %calculate proboscis angles and distance in 0 and 90 situations
        meanAngles0(i)=nanstd(probAngle(alphas>=margin/2 & alphas<=180-margin/2));
        meanAngles90(i)=nanstd(probAngle(alphas>=90+margin/2 | alphas<=90-margin/2));
        meanBodyDist0(i,:)=nanmean(bodyPosDist((alphas>=margin/2 & alphas<=180-margin/2),:),1);
        meanBodyDist90(i,:)=nanmean(bodyPosDist((alphas>=90+margin/2 | alphas<=90-margin/2),:),1);

        %another way to get proboscis movement - relative to head!!
        specM_probDist(i,:,:) = fftMovement(bodyPosDist,fNew,Fs,boutLength,minLengthFFT);
        specM_probDist0(i,:,:) = fftMovement(bodyPosDist((alphas>=margin/2 & alphas<=180-margin/2),:),fNew,Fs,boutLength,minLengthFFT);
        specM_probDist90(i,:,:) = fftMovement(bodyPosDist((alphas>=90+margin/2 | alphas<=90-margin/2),:),fNew,Fs,boutLength,minLengthFFT);

        %calculate isolated proboscis position for aligned (0) and not
        %aligned (90)
        allIsolProbPos0(:,:,i)=nansum(tempMatrix(:,:,(alphas>=margin/2 & alphas<=180-margin/2)),3);
        allIsolProbPos90(:,:,i)=nansum(tempMatrix(:,:,(alphas>=90+margin/2 | alphas<=90-margin/2)),3);


        %% fourier analysis of movement
        % (in both cardinal directions

        specM_probNorm(i,:,:) = fftMovement(probNormR,fNew,Fs,boutLength,minLengthFFT);
        specM_prob(i,:,:) = fftMovement(allprobR,fNew,Fs,boutLength,minLengthFFT);
        specM_head(i,:,:) = fftMovement(allHeadR,fNew,Fs,boutLength,minLengthFFT);
        specM_thorax(i,:,:) = fftMovement(allThoraxR,fNew,Fs,boutLength,minLengthFFT);

        %make sure that both proboscis data have entries, and are the same
        %length, otherwise keep nan inputs
        if sum(~isnan(probNormR0(:)))~=0 || sum(~isnan(probNormR90(:)))~=0
            if sum(sum(~isnan(probNormR0)))>10
                specM_probNorm0(i,:,:) = fftMovement(probNormR0,fNew,Fs,boutLength,minLengthFFT);
                %test: take the proboscis movement (NOT isolated) in the
                %two movement phases
                specM_prob0(i,:,:) = fftMovement(probR0,fNew,Fs,boutLength,minLengthFFT);
                specM_head0(i,:,:) = fftMovement(headR0,fNew,Fs,boutLength,minLengthFFT);
                specM_thorax0(i,:,:) = fftMovement(thoraxR0,fNew,Fs,boutLength,minLengthFFT);
                %calculate absolute movement
                % [diff(probNormR0(:,1)) diff(probNormR0(:,2))]
                %reconcile absolute movement spectra with X and Y
                %spectra!!!!!!!!!!!
                absMov_prob0{i}=euclid_dist(probR0);
                absMov_head0{i}=euclid_dist(headR0);
                absMov_thorax0{i}=euclid_dist(thoraxR0);
                temp = fftMovement([euclid_dist(probR0) euclid_dist(probR0)],fNew,Fs,boutLength,minLengthFFT);
                specM_probAbs(i,:,1)=temp(:,1);
                temp = fftMovement([absMov_thorax0{i} absMov_thorax0{i}],fNew,Fs,boutLength,minLengthFFT);
                specM_thoraxAbs(i,:,1)=temp(:,1);
            end

            if sum(sum(~isnan(probNormR90)))>10
                specM_probNorm90(i,:,:) = fftMovement(probNormR90,fNew,Fs,boutLength,minLengthFFT);
                %test: take the proboscis movement (NOT isolated) in the
                %two movement phases
                specM_prob90(i,:,:) = fftMovement(probR90,fNew,Fs,boutLength,minLengthFFT);
                specM_head90(i,:,:) = fftMovement(headR90,fNew,Fs,boutLength,minLengthFFT);
                specM_thorax90(i,:,:) = fftMovement(thoraxR90,fNew,Fs,boutLength,minLengthFFT);
                %calculate absolute movement
                absMov_prob90{i}=euclid_dist(probR90);
                absMov_head90{i}=euclid_dist(headR90);
                absMov_thorax90{i}=euclid_dist(thoraxR90);
                temp = fftMovement([absMov_prob90{i} absMov_prob90{i}],fNew,Fs,boutLength,minLengthFFT);
                specM_probAbs(i,:,2)=temp(:,1);
                temp = fftMovement([absMov_thorax90{i} absMov_thorax90{i}],fNew,Fs,boutLength,minLengthFFT);
                specM_thoraxAbs(i,:,2)=temp(:,1);
            end
        end


        %generate 10 dummy data entries with scrambled proboscis for each
        %file
        generateRand=0;
        meanProbPos=1;
        if meanProbPos==1
            numReps=1;
        else 
            numReps=10;
        end

        if generateRand==1
            for r=1:numReps

                %select random isolated proboscis - use only non-nan entries
                ProbNormRnoNan=probNormR(~isnan(allprobR(:,1)),:);
                randomNonNanIndices=1:size(ProbNormRnoNan,1);randomNonNanIndices=randperm(length(randomNonNanIndices));
                tempProb=allprobR;
                tempProb(~isnan(allprobR(:,1)),:)=ProbNormRnoNan(randomNonNanIndices,:);

                if meanProbPos==1
                %use mean proboscis position (x/y) for all frames
                tempProb=repmat([nanmean(probNormR,1)],size(allprobR,1),1);
                end

                %turn back to original position
                for t=1:size(tempProb,1)
                    Rback=-[cosd(alphas(t)) -sind(alphas(t)); sind(alphas(t)) cosd(alphas(t))];
                    tempProbR(t,:)=[Rback*tempProb(t,:)']';%now I did NOT flip heads turned down extra!!
                end
                %add head indices
                tempProbR2=tempProbR+flipInd.*allHeadR;
                %seems like need last flip around the origin Yaxis for all data
                allprobR=-tempProbR2;

                if meanProbPos==1
                save(['randProbMean_',filenames{i}(indT:end-4),'_',num2str(r),'.mat'],'allDataR','allHeadR','allprobR','allR','allThoraxR','labels','outlineR','patternR');
                else
                save(['randProb_',filenames{i}(indT:end-4),'_',num2str(r),'.mat'],'allDataR','allHeadR','allprobR','allR','allThoraxR','labels','outlineR','patternR');
                end
            end
        end
    end

    %% store general parameters

        xpos_bout_all=[xpos_bout_all;xpos_bout i*ones(size(xpos_bout,1),1)];

end




%% save data
MpixPerMM=nanmean(pixPerMM);
currPath=pwd;inds=strfind(currPath,'\');



%% Plot movement spectra

% original body data
f1=figure('Position',[500 500 1200 400]);
subplot(1,3,1)
shadedErrorBar_anna(fNew+10^-10,nanmean(specM_prob(:,:,1),1)/MpixPerMM,nanstd(specM_prob(:,:,1),1,1)/(sqrt(size(specM_prob,1))-1)*2/MpixPerMM,'b',0.5,true);hold on;
shadedErrorBar_anna(fNew+10^-10,nanmean(specM_prob(:,:,2),1)/MpixPerMM,nanstd(specM_prob(:,:,2),1,1)/(sqrt(size(specM_prob,1))-1)*2/MpixPerMM,'c',0.5,true);
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([fmin fmax]);
ylim([0.01 50])
ylabel('amplitude (mm)')
xlabel('frequency (Hz)')
title('proboscis')
myLegend({'parall. pattern', 'perp. pattern'},{'c-','b-'})

subplot(1,3,2)
shadedErrorBar_anna(fNew+10^-10,nanmean(specM_head(:,:,1),1)/MpixPerMM,nanstd(specM_head(:,:,1),1,1)/(sqrt(size(specM_head,1))-1)*2/MpixPerMM,'b',0.5,true);hold on;
shadedErrorBar_anna(fNew+10^-10,nanmean(specM_head(:,:,2),1)/MpixPerMM,nanstd(specM_head(:,:,2),1,1)/(sqrt(size(specM_head,1))-1)*2/MpixPerMM,'c',0.5,true);
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([fmin fmax]);
ylim([0.01 50])
ylabel('amplitude (mm)')
xlabel('frequency (Hz)')
title('head')
myLegend({'parall. pattern', 'perp. pattern'},{'c-','b-'})

subplot(1,3,3)
shadedErrorBar_anna(fNew+10^-10,nanmean(specM_thorax(:,:,1),1)/MpixPerMM,nanstd(specM_thorax(:,:,1),1,1)/(sqrt(size(specM_thorax,1))-1)*2/MpixPerMM,'b',0.5,true);hold on;
shadedErrorBar_anna(fNew+10^-10,nanmean(specM_thorax(:,:,2),1)/MpixPerMM,nanstd(specM_thorax(:,:,2),1,1)/(sqrt(size(specM_thorax,1))-1)*2/MpixPerMM,'c',0.5,true);
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([fmin fmax]);
ylim([0.01 50])
ylabel('amplitude (mm)')
xlabel('frequency (Hz)')
title('thorax')
myLegend({'parall. pattern', 'perp. pattern'},{'c-','b-'})

print -f1 -dpdf -r300 -painters -bestfit absoluteMov.eps

%plot individual frequency spectra
% figure;plot(fNew+10^-10,specM_thorax(1,:,2)/MpixPerMM);
% set(gca,'xscale','log');set(gca,'yscale','log')
% xlim([50 100]);ylim([0.001 1])

%plot movement ratios
f2=figure('Position',[500 500 1200 400]);
subplot(1,3,1);hold on;
% shadedErrorBar_anna(fNew+10^-10,nanmean(specM_probDist(:,:,2),1),nanstd(specM_probDist(:,:,2),1,1)/(sqrt(size(specM_probDist,2))-1)*2,'b',0.5,true);
shadedErrorBar_anna(fNew+10^-10,nanmean(specM_thorax(:,:,1),1)/MpixPerMM,nanstd(specM_thorax(:,:,1)/MpixPerMM,1,1)/(sqrt(size(specM_probDist,1))-1)*2,'y',0.5,true);
shadedErrorBar_anna(fNew+10^-10,nanmean(specM_thorax(:,:,2),1)/MpixPerMM,nanstd(specM_thorax(:,:,2),1,1)/MpixPerMM/(sqrt(size(specM_thorax,1))-1)*2,'b',0.5,true);
shadedErrorBar_anna(fNew+10^-10,nanmean(specM_prob(:,:,1),1)/MpixPerMM,nanstd(specM_prob(:,:,1)/MpixPerMM,1,1)/(sqrt(size(specM_prob,1))-1)*2,'m',0.5,true);
shadedErrorBar_anna(fNew+10^-10,nanmean(specM_prob(:,:,2),1)/MpixPerMM,nanstd(specM_prob(:,:,2)/MpixPerMM,1,1)/(sqrt(size(specM_prob,1))-1)*2,'r',0.5,true);
% shadedErrorBar_anna(fNew+10^-10,nanmean(specM_probDist(:,:,3),1),nanstd(specM_probDist(:,:,3),1,1)/(sqrt(size(specM_probDist,3))-1)*2,'m',0.5,true);
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([fmin fmax]);
ylim([0.005 10])
ylabel('amplitude (mm)')
xlabel('frequency (Hz)')
title('proboscis movement')
% myLegend({'prob - head', 'head - thorax', 'prob - thorax'},{'c-','b-','m-'})
myLegend({'thoraxX', 'thoraxy', 'probX','probY'},{'y-','b-','m-','r-'})

% ratio between prob and head and head and thorax
subplot(1,3,2)
% norm_data1=specM_probDist(:,:,1)./repmat(nanmax(specM_probDist(:,:,1),[],2),1,size(specM_probDist(:,:,1),2));
% norm_data2=specM_probDist(:,:,2)./repmat(nanmax(specM_probDist(:,:,2),[],2),1,size(specM_probDist(:,:,2),2));
norm_dataA1=specM_prob(:,:,1);
norm_dataA2=specM_thorax(:,:,1);%I
norm_dataB1=specM_prob(:,:,2);
norm_dataB2=specM_thorax(:,:,2);%I
% found this in 20230630, but not sure it is original
% norm_data1=specM_probDist(:,:,2)./repmat(nanmax(specM_probDist(:,:,2),[],2),1,size(specM_probDist(:,:,2),2));
% norm_data2=specM_thorax(:,:,1)./repmat(nanmax(specM_thorax(:,:,1),[],2),1,size(specM_thorax(:,:,1),2));%I
ratioMeanA=nanmean(norm_dataA1./norm_dataA2,1);ratiostdA=nanstd(norm_dataA1./norm_dataA2,1,1);
shadedErrorBar_anna(fNew+10^-10,ratioMeanA,ratiostdA/(sqrt(size(specM_prob,1))-1),'m',0.5,true);hold on
ratioMeanB=nanmean(norm_dataB1./norm_dataB2,1);ratiostdB=nanstd(norm_dataB1./norm_dataB2,1,1);
shadedErrorBar_anna(fNew+10^-10,ratioMeanB,ratiostdB/(sqrt(size(specM_prob,1))-1),'r',0.5,true);hold on

plot([1 50],[1 1],'k')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([fmin fmax]);
ylim([.1 10])
ylabel('ratio mov. (prob/thorax)')
xlabel('frequency (Hz)')
title('thorax')
myLegend({'X', 'Y'},{'m-','r-'})

subplot(1,3,3)
% norm_data1=specM_probDist(:,:,1)./repmat(nanmax(specM_probDist(:,:,1),[],2),1,size(specM_probDist(:,:,1),2));
% norm_data2=specM_probDist(:,:,2)./repmat(nanmax(specM_probDist(:,:,2),[],2),1,size(specM_probDist(:,:,2),2));
norm_dataA1=specM_prob(:,:,1);
norm_dataA2=specM_head(:,:,1);%I
norm_dataB1=specM_prob(:,:,2);
norm_dataB2=specM_head(:,:,2);%I
% found this in 20230630, but not sure it is original
% norm_data1=specM_probDist(:,:,2)./repmat(nanmax(specM_probDist(:,:,2),[],2),1,size(specM_probDist(:,:,2),2));
% norm_data2=specM_thorax(:,:,1)./repmat(nanmax(specM_thorax(:,:,1),[],2),1,size(specM_thorax(:,:,1),2));%I
ratioMeanA=nanmean(norm_dataA1./norm_dataA2,1);ratiostdA=nanstd(norm_dataA1./norm_dataA2,1,1);
shadedErrorBar_anna(fNew+10^-10,ratioMeanA,ratiostdA/(sqrt(size(specM_prob,1))-1),'m',0.5,true);hold on
ratioMeanB=nanmean(norm_dataB1./norm_dataB2,1);ratiostdB=nanstd(norm_dataB1./norm_dataB2,1,1);
shadedErrorBar_anna(fNew+10^-10,ratioMeanB,ratiostdB/(sqrt(size(specM_prob,1))-1),'r',0.5,true);hold on

plot([1 50],[1 1],'k')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([fmin fmax]);
ylim([.1 10])
ylabel('ratio mov. (prob/head)')
xlabel('frequency (Hz)')
title('head')
myLegend({'X', 'Y'},{'m-','r-'})

print -f2 -dpdf -r300 -painters -bestfit ratioMovThxProb.eps


%
%% isolated proboscis

%distance between proboscis and head - should show movement of proboscis
%only (second entry is distance between head and thorax, should only
%marginally change)
%!!! Careful! The probDist is already calclated in mm, the isol. proboscis
%and other measures are still in px units

f3=figure('Position',[500 500 1200 400]);
subplot(1,3,1);hold on;
% shadedErrorBar_anna(fNew+10^-10,nanmean(specM_probDist(:,:,2),1),nanstd(specM_probDist(:,:,2),1,1)/(sqrt(size(specM_probDist,2))-1)*2,'b',0.5,true);
shadedErrorBar_anna(fNew+10^-10,nanmean(specM_probDist(:,:,1),1)/MpixPerMM,nanstd(specM_probDist(:,:,1)/MpixPerMM,1,1)/(sqrt(size(specM_probDist,1))-1)*2,'y',0.5,true);
shadedErrorBar_anna(fNew+10^-10,nanmean(specM_probNorm(:,:,1),1)/MpixPerMM,nanstd(specM_probNorm(:,:,1)/MpixPerMM,1,1)/(sqrt(size(specM_probNorm,1))-1)*2,'m',0.5,true);
shadedErrorBar_anna(fNew+10^-10,nanmean(specM_probNorm(:,:,2),1)/MpixPerMM,nanstd(specM_probNorm(:,:,2)/MpixPerMM,1,1)/(sqrt(size(specM_probNorm,1))-1)*2,'r',0.5,true);
% shadedErrorBar_anna(fNew+10^-10,nanmean(specM_probDist(:,:,3),1),nanstd(specM_probDist(:,:,3),1,1)/(sqrt(size(specM_probDist,3))-1)*2,'m',0.5,true);
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([fmin fmax]);
ylim([0.01 10])
ylabel('amplitude (mm)')
xlabel('frequency (Hz)')
title('proboscis movement')
% myLegend({'prob - head', 'head - thorax', 'prob - thorax'},{'c-','b-','m-'})
myLegend({'prob - head', 'probIsolx','probIsoly'},{'y-','m-','r-'})

% ratio between prob X and Y
subplot(1,3,2)
norm_dataA1=specM_prob(:,:,2);
norm_dataA2=specM_prob(:,:,1);%I

% found this in 20230630, but not sure it is original
% norm_data1=specM_probDist(:,:,2)./repmat(nanmax(specM_probDist(:,:,2),[],2),1,size(specM_probDist(:,:,2),2));
% norm_data2=specM_thorax(:,:,1)./repmat(nanmax(specM_thorax(:,:,1),[],2),1,size(specM_thorax(:,:,1),2));%I
ratioMeanA=nanmean(norm_dataA1./norm_dataA2,1);ratiostdA=nanstd(norm_dataA1./norm_dataA2,1,1);
shadedErrorBar_anna(fNew+10^-10,ratioMeanA,ratiostdA/(sqrt(size(specM_prob,1))-1),'k',0.5,true);hold on

plot([1 50],[1 1],'k')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([fmin fmax]);
ylim([.1 10])
ylabel('ratio mov. (prob Y/X)')
xlabel('frequency (Hz)')
title('proboscis')

%total movement extracted from integral of FFT
subplot(1,3,3)
hold on

binWidth=nanmean(diff(fNew));
% binRange=find(fNew>=0.1 & fNew<=3);
binRange=find(fNew>=1 & fNew<=50);
xMov=nansum(specM_probNorm(:,binRange,1),2)/MpixPerMM*binWidth;
yMov=nansum(specM_probNorm(:,binRange,2),2)/MpixPerMM*binWidth;

%reorganise table so that data from same animals is in same row of two
%columns file
uniqueAnimal=unique(animal.*date);
n=length(uniqueAnimal);
dataXvsY_sameAnimal=nan(n,2);%first column is 0, second is 90
%make versatile by allowing different inputs
data1=xMov;
data2=yMov;

dataXvsY_sameAnimal=[data1 data2];
dataXvsY_sameAnimal(((dataXvsY_sameAnimal(:,1)==0 ).* (dataXvsY_sameAnimal(:,2)==0))==1,:)=nan;

%find indices without entries
if sum(~isnan(dataXvsY_sameAnimal(:)))>0
    if size(dataXvsY_sameAnimal,1)>1
        boxplot(dataXvsY_sameAnimal,'labels',{'perp','parallel'},'colors','kk')
        for u=1:size(dataXvsY_sameAnimal,1)
            plot([1+0.01*randn(1,1),2+0.01*randn(1,1)],dataXvsY_sameAnimal(u,:),'k.-','Markersize',10);
        end
    else
        plot([1 2],dataXvsY_sameAnimal,'ro');xlim([0 3]);
        set(gca,'xtick',[1 2],'xticklabel',{'X','Y'})
    end
    xlabel('orientation')
    ylabel('movement amplitude (mm)')
    xlim([0.5 2.5])
    set(gca,'yscale','log')
    ylim([floor(round(min(dataXvsY_sameAnimal(:))*10)/10) ceil(max(dataXvsY_sameAnimal(:))/10^floor(log10(max(dataXvsY_sameAnimal(:)))))*10^floor(log10(max(dataXvsY_sameAnimal(:))))])
end

%paired statistics
pairedData=dataXvsY_sameAnimal;
pairedData(isnan(pairedData(:,1)),:)=[];
pairedData(isnan(pairedData(:,2)),:)=[];

% paired or non-paired stats
if ~isempty(pairedData(:,1)) && ~isempty(pairedData(:,2))
    if ~isempty(pairedData)
        [p1,h,stats1] = signrank(pairedData(:,1), pairedData(:,2));
        title(['X,ppaired = ',num2str(p1),' signedrank = ',num2str(stats1.signedrank)])

    else[p2,h,stats1] = ranksum(pairedData(:,1), pairedData(:,2));
        if exist('p2')
            title(['X,p-nopair = ',num2str(p2)])
        end
    end

    %print adata and animal ID
    minMaxThresh=[fSumMin;fSumMax;nan(length(uniqueAnimal)-2,1)];
    tbl01 = table(uniqueAnimal',dataXvsY_sameAnimal(:,1),dataXvsY_sameAnimal(:,2),minMaxThresh,'VariableNames',{'animalID','relMovX','relMovY','minmaxThresh'});
    currPath=pwd;inds=strfind(currPath,'\');
    writetable(tbl01,[currPath(inds(end)+1:end),'_relMovDist.csv'])
end

print(f3,'isolProbMovDist.eps','-dpdf','-r300','-painters','-bestfit')
%
%% compare animals in M0 and M90
%for comparing M0 and M90, use only animals aligned with 0 or 90 depending on condition!
if patternRot==0
    dataSpec=specM_probNorm0;
elseif patternRot==90
    dataSpec=specM_probNorm90;
else
    %for freely rotating animals use this:
    dataSpec=specM_probNorm;
end
%
f4=figure('Position',[500 500 1200 400]);
subplot(1,3,1)
shadedErrorBar_anna(fNew+10^-10,nanmean(dataSpec(:,:,1),1)/MpixPerMM,nanstd(dataSpec(:,:,1),1,1)/(sqrt(size(dataSpec,1))-1)*2/MpixPerMM,'m',0.5,true);hold on;
shadedErrorBar_anna(fNew+10^-10,nanmean(dataSpec(:,:,2),1)/MpixPerMM,nanstd(dataSpec(:,:,2),1,1)/(sqrt(size(dataSpec,1))-1)*2/MpixPerMM,'r',0.5,true);
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([fmin fmax]);ylim([0.01 10])
ylabel('amplitude (mm)')
xlabel('frequency (Hz)')
title('isol. proboscis')
myLegend({'parall. bodyaxis', 'perp. bodyaxis'},{'r-','m-'})

subplot(1,3,3)
hold on

%total movement extracted from integral of FFT
binWidth=nanmean(diff(fNew));
binRange=find(fNew>=fSumMin & fNew<=fSumMax);
% calculate movement from fft - smaller movements get relatively more weight
xMov=nansum(dataSpec(:,binRange,1),2)/MpixPerMM*binWidth;
yMov=nansum(dataSpec(:,binRange,2),2)/MpixPerMM*binWidth;
xMov(yMov==0)=[];yMov(yMov==0)=[];

if size(xMov,1)>1
    boxplot([xMov, yMov],'labels',{'parall.','perp.'},'colors','mr')
    for u=1:size(xMov,1)
        plot([1+0.01*randn(1,1),2+0.01*randn(1,1)],[xMov(u),yMov(u)],'k.-','Markersize',10);
    end
else
    if ~isempty(xMov)
        plot([1 2],[xMov, yMov],'ro');xlim([0 3]);
        set(gca,'xtick',[1 2],'xticklabel',{'parall.','perp.'})
    end
end
xlabel('movement direction')
ylabel('total movement (mm)')
ylim([0 10])

if ~isempty(xMov)
    [p1,h,stats1] = signrank(xMov, yMov);
    title(['p = ',num2str(p1)])
end

subplot(1,3,2);hold on;
ratioMean=nanmean(dataSpec(:,:,2)./dataSpec(:,:,1),1);
ratiostd=nanstd(dataSpec(:,:,2)./dataSpec(:,:,1),1,1);
shadedErrorBar_anna(fNew+10^-10,ratioMean,ratiostd/(sqrt(size(dataSpec,1))-1),'k',0.5,true);
plot([1 50],[1 1],'k')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([fmin fmax]);ylim([0.1 10])
ylabel('ratio mov. (parall./perp. bodyaxis)')
xlabel('frequency (Hz)')
title('proboscis')

% print -f4 -dpdf -r300 -painters -bestfit isolProbMov.eps
print(f4,'isolProbMov090.eps','-dpdf','-r300','-painters','-bestfit')

%% rel. movement isolated proboscis in parallel and perpendicular

f51=figure('Position',[500 500 1200 400]);
subplot(1,3,1)
shadedErrorBar_anna(fNew+10^-10,nanmean(specM_probNorm0(:,:,1),1)/MpixPerMM,nanstd(specM_probNorm0(:,:,1),1,1)/(sqrt(size(specM_probNorm0,1))-1)*2/MpixPerMM,'g',0.5,true);hold on;
shadedErrorBar_anna(fNew+10^-10,nanmean(specM_probNorm90(:,:,1),1)/MpixPerMM,nanstd(specM_probNorm90(:,:,1),1,1)/(sqrt(size(specM_probNorm90,1))-1)*2/MpixPerMM,'r',0.5,true);
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([fmin fmax]);
ylim([0.01 50])
ylabel('amplitude (mm)')
xlabel('frequency (Hz)')
title('perp. body axis')
myLegend({'parall. pattern', 'perp. pattern'},{'g-','r-'})

subplot(1,3,2)
shadedErrorBar_anna(fNew+10^-10,nanmean(specM_probNorm0(:,:,2),1)/MpixPerMM,nanstd(specM_probNorm0(:,:,2),1,1)/(sqrt(size(specM_probNorm0,1))-1)*2/MpixPerMM,'g',0.5,true);hold on;
shadedErrorBar_anna(fNew+10^-10,nanmean(specM_probNorm90(:,:,2),1)/MpixPerMM,nanstd(specM_probNorm90(:,:,2),1,1)/(sqrt(size(specM_probNorm90,1))-1)*2/MpixPerMM,'r',0.5,true);
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([fmin fmax]);
ylim([0.01 50])
ylabel('amplitude (mm)')
xlabel('frequency (Hz)')
title('parallel body axis')
myLegend({'parall. pattern', 'perp. pattern'},{'g-','r-'})

print(f51,'isolProbMovHeading_abs.eps','-dpdf','-r300','-painters','-bestfit')


f5=figure('Position',[500 500 1200 400]);
subplot(1,3,1);hold on;
ratioMean0=nanmean(specM_probNorm0(:,:,2)./specM_probNorm0(:,:,1),1);
ratiostd0=nanstd(specM_probNorm0(:,:,2)./specM_probNorm0(:,:,1),1,1);
shadedErrorBar_anna(fNew+10^-10,ratioMean0,ratiostd0/(sqrt(size(specM_probNorm0,1))-1),'k',0.5,true);
plot([1 50],[1 1],'k')
set(gca,'xscale','log')
set(gca,'yscale','log')

xlim([fmin fmax]);ylim([0.1 10])
ylabel('isolated proboscis (parall./perp. bodyaxis)')
xlabel('frequency (Hz)')
title('body orientation 0°')

subplot(1,3,2);hold on;
ratioMean90=nanmean(specM_probNorm90(:,:,2)./specM_probNorm90(:,:,1),1);
ratiostd90=nanstd(specM_probNorm90(:,:,2)./specM_probNorm90(:,:,1),1,1);
shadedErrorBar_anna(fNew+10^-10,ratioMean90,ratiostd90/(sqrt(size(specM_probNorm90,1))-1),'k',0.5,true);
plot([1 50],[1 1],'k')
set(gca,'xscale','log')
set(gca,'yscale','log')

xlim([fmin fmax]);ylim([0.1 10])
ylabel('isolated proboscis (parall./perp. bodyaxis)')
xlabel('frequency (Hz)')
title('body orientation 90°')


%total movement extracted from integral of FFT
subplot(1,3,3)
hold on

binWidth=nanmean(diff(fNew));
% binRange=find(fNew>=fSumMin & fNew<=fSumMax);
binRange=find(fNew>=0.1 & fNew<=3);
binRange=find(fNew>=15 & fNew<=50);
xMov0=nansum(specM_probNorm0(:,binRange,1),2)/MpixPerMM*binWidth;
xMov90=nansum(specM_probNorm90(:,binRange,1),2)/MpixPerMM*binWidth;
yMov0=nansum(specM_probNorm0(:,binRange,2),2)/MpixPerMM*binWidth;
yMov90=nansum(specM_probNorm90(:,binRange,2),2)/MpixPerMM*binWidth;
relMov0=yMov0./xMov0;
relMov90=yMov90./xMov90;

%reorganise table so that data from same animals is in same row of two
%columns file
%also if possible use data from pattern 0 orientation for 0 turn data and
%same for 90 (not animal turned 0 with 90 orientation pattern, as this is
%not equivalent)
uniqueAnimal=unique(animal);
n=length(unique(animal));
data0vs90_sameAnimal=nan(n,2);%first column is 0, second is 90
%make versatile by allowing different inputs
data1=relMov0;
data2=relMov90;
% if the pattern is all the same, data HAS to come from same animals
if length(unique(patterns))==1
    data0vs90_sameAnimal=[data1 data2];
else
    %run through data and sort everything
    for l=1:n
        indA=find(animal==uniqueAnimal(l));
        %check if entries in Mov0 from pattern 0
        ind0=data1(indA).*(patterns(indA)==0)';
        ind0(ind0==0)=nan;%0 values come from "wrong" pattern orientation
        if sum(~isnan(ind0))>0
            data0vs90_sameAnimal(l,1)=nanmean(ind0(~isnan(ind0)));
        end
        %check if entries in Mov0 from pattern 0
        ind90=data2(indA).*(patterns(indA)==90)';
        ind90(ind90==0)=nan;%0 values come from "wrong" pattern orientation
        if sum(isnan(ind90))>0
            data0vs90_sameAnimal(l,2)=nanmean(ind90(~isnan(ind90)));
        end
    end
end

%find indices without entries
% indNan=unique([find(relMov0==0); find(relMov90==0)]);
% relMov0(indNan)=[];relMov90(indNan)=[];
if sum(~isnan(data0vs90_sameAnimal(:)))>0
    if size(data0vs90_sameAnimal,1)>1
        boxplot(data0vs90_sameAnimal,'labels',{'0°','90°'},'colors','kk')
        for u=1:size(data0vs90_sameAnimal,1)
            plot([1+0.01*randn(1,1),2+0.01*randn(1,1)],data0vs90_sameAnimal(u,:),'k.-','Markersize',10);
        end
    else
        plot([1 2],data0vs90_sameAnimal,'ro');xlim([0 3]);
        set(gca,'xtick',[1 2],'xticklabel',{'0°','90°'})
    end
    xlabel('orientation')
    ylabel('rel. movement')
    ylim([0 1.1*max(data0vs90_sameAnimal(:))])
end

%paired statistics
pairedData=data0vs90_sameAnimal;
pairedData(isnan(pairedData(:,1)),:)=[];
pairedData(isnan(pairedData(:,2)),:)=[];
% if ~isempty(pairedData)
% [p1,h,stats1] = signrank(pairedData(:,1), pairedData(:,2));
% %     title(['X,ppaired = ',num2str(p1)])
% end
% non-paired stats
if ~isempty(pairedData(:,1)) && ~isempty(pairedData(:,2))
    [p2,h,stats1] = ranksum(data0vs90_sameAnimal(:,1), data0vs90_sameAnimal(:,2));
    if exist('p2')
        title(['X,p-nopair = ',num2str(p2)])
    end


    print(f5,'isolProbMovHeading.eps','-dpdf','-r300','-painters','-bestfit')

    %print absolute and relative movement and animal ID
    %collect relative proboscis movement in different orientations into table
    % tbl01 = table(animal',patterns',relMov0,relMov90,'VariableNames',{'animalID','patternOri','relMov0','relMov90'});
    minMaxThresh=[fSumMin;fSumMax;nan(length(uniqueAnimal)-2,1)];
    tbl01 = table(uniqueAnimal',data0vs90_sameAnimal(:,1),data0vs90_sameAnimal(:,2),minMaxThresh,'VariableNames',{'animalID','relMov0','relMov90','minmaxThresh'});
    currPath=pwd;inds=strfind(currPath,'\');
    writetable(tbl01,[currPath(inds(end)+1:end),'_relMov0and90.csv'])
end
%


%% plot isolated proboscis position

f9=figure;hold on;
%plot individual pixels
colormap(flipud(bone))
imagesc(nanmean(allIsolProbPos90,3));

% make contour plot with smooth gradient of data and contour lines
colormap((inferno)); hold on
% 
% %should ideally convert these to mm before averaging from different camera
% %views!!!!!!!!!!!!!!
for u=1:size(allIsolProbPos,3)
    allIsolProbPosN(:,:,u)=allIsolProbPos(:,:,u)/quantile(reshape(allIsolProbPos(:,:,u).',1,[]),.9999);
    allIsolProbPos0N(:,:,u)=allIsolProbPos0(:,:,u)/quantile(reshape(allIsolProbPos0(:,:,u).',1,[]),.9999);
    allIsolProbPos90N(:,:,u)=allIsolProbPos90(:,:,u)/quantile(reshape(allIsolProbPos90(:,:,u).',1,[]),.9999);
end

%smooth and normalise contour data
contourData=nanmean(allIsolProbPos90N,3);
contourData=imgaussfilt(contourData,2.5);
contourData=contourData/quantile(contourData(:),.99);%normalise not to peak but max 97.5%
imagesc(contourData);

%generate contour lines
contour(contourData,[0:.1:1],'-w');
caxis([0 1]) % to make axis the same for M0 and M90
maxPix=120;

%adjust axes to mm
maxMMrange=ceil(maxPix/2/MpixPerMM);
%get uneven number for the range distribution
intRange=ceil((maxMMrange*2)/7);
ticklabelMM=(-maxMMrange:intRange:maxMMrange);
ticklabel=round(ticklabelMM*MpixPerMM);
plot([60 60],[0 maxPix],'--w')
% ticklabel=(0:20:maxPix)-maxPix/2;
% ticklabelMM=ticklabel/MpixPerMM;
xlim([0 maxPix]);ylim([0 maxPix]);
set(gca,'xtick',[ticklabel+maxPix/2],'xticklabel',[ticklabelMM])
set(gca,'ytick',[ticklabel+maxPix/2],'yticklabel',[ticklabelMM+maxMMrange])
ylabel('position (mm)');xlabel('position (mm)');
axis square
colorbar

print(f9,'probHisto90.eps','-dpdf','-r300','-painters','-bestfit')

%plot histogram of isolated proboscis data collapsed along the Y-Axis
% to only show the positions along the x axis
xpos=nan(size(allIsolProbPos,2),size(allIsolProbPos,3));
for i=1:size(allIsolProbPos,3)
    xpos(:,i)=nansum(allIsolProbPos(:,:,i),1);
end


f10=figure;hold on;
%normalise xposition to fairly combine data
xpos_norm=xpos./repmat(nansum(xpos,1),size(xpos,1),1);%normalise to sum under the curve
scaleMax=nanmax(nanmean(xpos_norm,2));
shadedErrorBar_anna(1:size(xpos,1),nanmean(xpos_norm,2)/scaleMax, nanstd(xpos_norm,1,2)/(sqrt(size(xpos,2))-1)/scaleMax)
plot([60 60],[0 maxPix],'--')
xlim([0 maxPix]);
ylim([0 max(nanmean(xpos_norm,2))/scaleMax+max(nanstd(xpos_norm,1,2)/(sqrt(size(xpos,2))-1))/scaleMax]);
set(gca,'xtick',[ticklabel+maxPix/2],'xticklabel',[ticklabelMM])
ylabel('rel. number of contacts');xlabel('position (mm)');
print(f10,'isolProb_xHisto.eps','-dpdf','-r300','-painters','-bestfit')
print(f10,'isolProb_xHisto.jpg','-djpeg','-r300','-painters')

%save median isolated proboscis position relative to midline of animal
%need to get the indices, not the values, of the peak (max) of the curve:
% smooth the curves before getting the max
for u=1:size(xpos,2)
    smoothXpos(:,u)=resample(smooth(xpos(:,u),17,'sgolay'),100,1);
end


%reample smoothXpos to get higher resolution and sample back before
%converting to mm -- when resolution to low and smoothing to high, get the
%same peak values from different distributions
%analysing milimeters
[dat indMed] = min(abs(smoothXpos-repmat(nanmax(smoothXpos,[],1),size(smoothXpos,1),1)),[],1);%median does not work because of all the 0 entries
%calculate relative to midline of animal, which is set to 60
prob_pos_midline_tbl=table(animal',((indMed/100)'-maxPix/2)/MpixPerMM,'VariableNames',{'animal','isol_prob_pos_mm'})
writetable(prob_pos_midline_tbl,[currPath(inds(end)+1:end),'_isol_prob_pos.csv'])
% dur_length_data=dur_length_tbl{:,:};save([currPath(inds(end)+1:end),'_dur_length.mat'],'dur_length_data');

%do the same for individual bouts
for u=1:size(xpos_bout_all,1)
    smoothXpos_bout(u,:)=resample(smooth(xpos_bout_all(u,1:end-1),17,'sgolay'),1,1);
end
[dat indMed] = min(abs(smoothXpos_bout-repmat(nanmax(smoothXpos_bout(:,1:end-1),[],2),1,size(smoothXpos_bout,2))),[],2);%median does not work because of all the 0 entries
%calculate relative to midline of animal, which is set to 60
isolPprobPos_bouts=((indMed/2)'-maxPix/2)/MpixPerMM;%need to divide by 2, because made steps 0.5

%to test if analysis of proboscis position per bout works
figure;boxplot(isolPprobPos_bouts,xpos_bout_all(:,end),'plotstyle','compact');
hold on;plot([min(xpos_bout_all(:,end)) max(xpos_bout_all(:,end))],[0,0])
ylim([-15 15])

for i=1:length(unique(xpos_bout_all(:,end)))
    [p(i),h(i),~]=signrank(isolPprobPos_bouts(find(xpos_bout_all(:,end)==i)),0);
    n(i)=length(find(xpos_bout_all(:,end)==i));
end


% left_right_index=[nansum(xpos(1:60,:),1)' nansum(xpos(61:120,:),1)']


%% extract movement in absolute coordinates (x1/y1 to x2/y2)
function [medianM, histM]=extractMov(inputdata,binEdges)
movements=diff(inputdata,1);
movements(find(abs(movements)>50))=nan;%make sure to adjust this value when camera magnification changes!

%save histograms for amplitude and speed
histM=histc(abs(movements),binEdges)/size(movements,1);
medianM=nanmedian(abs(movements));
end

%% fourier analysis of movement
function specM = fftMovement(Fdata,fNew,Fs,boutLength,minLengthFFT)

%find "jumps" in data where proboscis did not touch
checkJumps=diff(isnan(Fdata(:,1)));
starts=find(checkJumps==-1);ends=find(checkJumps==1);

if ~isempty(starts) && ~isempty(ends)
    %if there is no nan in the beginning, but starts right with data,
    %make the first start entry 1
    if isnan(Fdata(1,1))==0
        starts=[1;starts];
    end
    %if the last index is a start, then make end the frame of the data
    if starts(end)>ends(end)
        ends=[ends;size(Fdata,1)];
    end

    lengthJumps=starts(2:end)-ends(1:end-1);
else
    lengthJumps=nan;
end

if ~isnan(lengthJumps)
    %interpolate nan values if mising bits are not longer than 0.25 ms
    for k=1:length(lengthJumps)
        if lengthJumps(k)<5
            startInp=ends(k)-2;
            endInp=starts(k+1)+2;
            %this is new 20230630 to accomodate different sizes of data
            tempData=Fdata(startInp:endInp,:);
            t    = 1:size(tempData,1);
            %run through data and interpolate missing values
            for r=1:size(Fdata,2)
                nanx = isnan(Fdata(startInp:endInp,r));
                tempData(nanx,r) = interp1(t(~nanx), tempData(~nanx,r), t(nanx));
            end

            %this is the old version
            %             tempData=Fdata(startInp:endInp,:);
            %             nanx1 = isnan(Fdata(startInp:endInp,1));
            %             nanx2 = isnan(Fdata(startInp:endInp,2));
            %             t    = 1:size(tempData,1);
            %             tempData(nanx1,1) = interp1(t(~nanx1), tempData(~nanx1,1), t(nanx1));
            %             tempData(nanx2,2) = interp1(t(~nanx2), tempData(~nanx2,2), t(nanx2));

            Fdata(startInp:endInp,:)=tempData;
        end
    end

    %delete all small pieces
    %find "jumps" in data where proboscis did not touch
    delIndices=[];
    checkJumps=diff(isnan(Fdata(:,1)));
    starts=find(checkJumps==-1);ends=find(checkJumps==1);
    %if there is no nan in the beginning, but starts right with data,
    %make the first start entry 1
    if isnan(Fdata(1,1))==0
        starts=[1;starts];
    end
    %if the last index is a start, then make end the frame of the data
    if starts(end)>ends(end)
        ends=[ends;size(Fdata,1)];
    end

    lengthData=ends(1:end)-starts(1:end);

    for k=1:length(lengthData)
        if lengthData(k)<boutLength %bouts need to be at least half a second long
            delIndices=[delIndices starts(k)-1:ends(k)+1];
        end
    end
    %make sure there are no 0 entries
    delIndices(delIndices==0)=[];
    %make sure index is not larger than data length
    delIndices(delIndices>size(Fdata,1))=[];

    if isempty(delIndices)==0
        Fdata(delIndices,:)=[];
    end

    % patch pieces together by ending one at the place where the other
    % starts
    %find "jumps" in data where proboscis did not touch
    checkJumps=diff(isnan(Fdata(:,1)));
    starts=find(checkJumps==-1);ends=find(checkJumps==1);

    if ~isempty(starts) && ~isempty(ends)

        %if there is no nan in the beginning, but starts right with data,
        %make the first start entry 1
        if isnan(Fdata(1,1))==0
            starts=[1;starts];
        end
        %if the last index is a start, then make end the frame of the data
        if starts(end)>ends(end)
            ends=[ends;size(Fdata,1)];
        end
        % lengthData=ends(1:end)-starts(1:end);
        %make every consecutive piece start where previous ended
        for k=2:length(starts)%for second solution start at 1, for first one at 2
            %difference between end previous and start
            diffStartEnd=Fdata(find(~isnan(Fdata(starts(k-1):ends(k-1),1)),1,'last')+starts(k-1)-1,:)-Fdata(find(~isnan(Fdata(starts(k):ends(k),1)),1,'first')+starts(k)-1,:);
            Fdata(starts(k):ends(k),:)=Fdata(starts(k):ends(k),:)+diffStartEnd;
        end
    end
end
%if we remove nan, generate data that doesn't exist but if we keep it, fft
%doesn't work. need to remove mean from all snippets and attach them to
%each other

Fdata(isnan(Fdata(:,1)),:)=[];%remmove nan
%need minimum length for FFT or return nonsensical data
if ~isempty(Fdata) && size(Fdata,1)>minLengthFFT

    %
    T = 1/Fs;             % Sampling period
    L = size(Fdata,1);             % Length of signal
    t = (0:L-1)*T;        % Time vector
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    %     f = Fs/2*linspace(0,1,NFFT/2+1);
    f = Fs*(0:floor(L/2))/L;%this is new from 20230630
    %     %remove moving average 20230630 - gives exactly the same result
    Fdata=Fdata-movmean(Fdata,400);

    Y=fft(Fdata-nanmean(Fdata,1));
    P2 = abs(Y/L);
    %     P1 = 2*P2(1:NFFT/2+1,:);%this is the amplitude
    P1 = 2*P2(1:floor(L/2)+1,:);%this is new from 20230630

    specM=interp1(f,P1,fNew);
else
    specM=nan(length(fNew),size(Fdata,2));
end
end

%% force a legend how i want it
function myLegend(labels,colours)
hold on;
for l=1:length(labels)
    L(l) = plot(nan, nan, colours{l});
end
legend(L, labels);legend boxoff
end