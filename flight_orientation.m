% flight_orientation extracts key descriptors of the animals' flight 
% during the approach, the probing phase and the departure. 
% These are forward and rotational flight speed and predominant rotation direction, 
% as well as body angle relative to the pattern axis. 
% The output is saved as a flightORI_filename.mat file.

clear all
close all

%load data: rotated data
[filenames, path]=uigetfile('allRotData*.mat','multiselect','on');
%if only one selected turn filenames into cell
if iscell(filenames)==0
    temp=filenames;clear filenames;
    filenames=cell(1,1);filenames{1}=temp;
end

Fs = 200;            % Sampling frequency
windowLength=0.25;      % window for flight orientation analysis in s. Duration of min. proboscis contact to analyse -
lengthGap=0.1;    %interpolate nan values if mising bits are not longer than lengthGap - do not consider this an individual approach
            
            
for i=1:length(filenames)
    cd(path)
    load(filenames{i})
    
    if size(allprobR,1)>500 && sum(isnan(allprobR(:,1))==0)>100
        if exist('allHeadR') && exist('allThoraxR')
            
            %% calculate the angle between head and thorax
            %use atan2d so that it goes through all 4 quadrants, i.e. from
            %-180 to 180 and not like with atand just 90 to -90
            alpha=atan2d((allHeadR(:,1)-allThoraxR(:,1)),(allHeadR(:,2)-allThoraxR(:,2)));
            
           
            
            %% calculate rotation speed
            rotSpeed=diff(alpha)*Fs;
            rotSpeed(rotSpeed>5000)=rotSpeed(rotSpeed>5000)-360*Fs;
            rotSpeed(rotSpeed<-5000)=rotSpeed(rotSpeed<-5000)+360*Fs;
            rotSpeed_smooth=sgolayfilt(rotSpeed,3,31);
            
            %             make handiness index
            turnCumSum(1)=nansum(rotSpeed_smooth>0);%cumulative sum of all path pieces that are turning left
            turnCumSum(2)=nansum(rotSpeed_smooth<0);%cumulative sum of all path pieces that are turning right
            
            
            
            
            
            %% sort proboscis touches and generate begin, middle, end indices
            
            %find "jumps" in data where proboscis did not touch
            [starts, ends]=checkDataPieces(allprobR(:,1));
            lengthJumps=starts(2:end)-ends(1:end-1);
            
            %interpolate nan values if mising bits are not longer than lengthGap
            for k=1:length(lengthJumps)
                if lengthJumps(k)<lengthGap*Fs
                    startInp=ends(k)-1;
                    endInp=starts(k+1)+1;if endInp>size(allprobR,1);endInp=size(allprobR,1);end
                    nanx1 = isnan(allprobR(startInp:endInp,1));nanx2 = isnan(allprobR(startInp:endInp,2));
                    tempData=allprobR(startInp:endInp,:);
                    t    = 1:size(tempData,1);
                    tempData(nanx1,1) = interp1(t(~nanx1), tempData(~nanx1,1), t(nanx1));
                    tempData(nanx2,2) = interp1(t(~nanx2), tempData(~nanx2,2), t(nanx2));
                    allprobR(startInp:endInp,:)=tempData;
                end
            end
            
            %delete all small pieces
            %find "jumps" in data where proboscis did not touch
            [starts, ends]=checkDataPieces(allprobR(:,1));
            lengthData=ends-starts;
            
            snipSize=2*windowLength*Fs;%this one needs to match the window length below! %ideally be at least 2x as long
            delIndices=[];
            for k=1:length(lengthData)
                if lengthData(k)<snipSize
                    if starts(k)==1
                        allprobR(starts(k):ends(k)+1,:)=nan;
                    else
                    allprobR(starts(k)-1:ends(k)+1,:)=nan;
                    end
                    
                end
            end
            
            % select the beginning, end and middle speed and rotational
            % angles for all snippets of proboscis contact
            [starts, ends]=checkDataPieces(allprobR(:,1));
            heading_snips=nan(length(starts),6);
            rotation_snips=nan(length(starts),3);
            
            startEndWind=windowLength*Fs;
            for u=1:length(starts)
                %beginning
                if starts(u)-startEndWind>1
                    beginHeading=alpha(starts(u)-startEndWind:starts(u)+startEndWind);
                    heading_snips(u,1)=circ_r(deg2rad(beginHeading(isnan(beginHeading)==0)));
                    heading_snips(u,2)=rad2deg(circ_mean(deg2rad(beginHeading(isnan(beginHeading)==0))));
                    rotation_snips(u,1)=nanmean(abs(rotSpeed_smooth(starts(u)-startEndWind:starts(u)+startEndWind)));
                end
                
                %middle
                middleHeading=alpha(starts(u)+startEndWind:ends(u)-startEndWind);
                heading_snips(u,4)=rad2deg(circ_mean(deg2rad(middleHeading)));
                
                %since length of middle heading can vary, make sure to make
                %a sliding average of vector strength with same length as
                %initial window
                if numel(middleHeading)>=startEndWind
                    for t=1:floor(numel(middleHeading)/startEndWind)
                        tempR(t)=circ_r(deg2rad(middleHeading(1+(t-1)*startEndWind:t*startEndWind))); %vector length running average
                    end
                    heading_snips(u,3)=nanmean(tempR);
                end
                rotation_snips(u,2)=nanmean(abs(rotSpeed_smooth(starts(u)+startEndWind:ends(u)-startEndWind)));
                
                %end
                if ends(u)+startEndWind<size(alpha,1)
                endHeading=alpha(ends(u)-startEndWind:ends(u)+startEndWind);
                heading_snips(u,5)=circ_r(deg2rad(endHeading(isnan(endHeading)==0))); %end
                heading_snips(u,6)=rad2deg(circ_mean(deg2rad(endHeading(isnan(endHeading)==0)))); %end
                rotation_snips(u,3)=nanmean(abs(rotSpeed_smooth(ends(u)-startEndWind:ends(u)+startEndWind))); %end
                end
            end
            
            %% PLOT
            
            f1=figure('Position',[200 200 1400 300]);
            subplot(1,4,[1 2 3]);
            time=1/Fs:1/Fs:numel(alpha)/Fs;
            hold on;plot(time,alpha,'.','MarkerSize',12);hold on;
            plot(time(isnan(allprobR(:,2))==0),alpha(isnan(allprobR(:,2))==0),'c.','MarkerSize',12);
            plot(time([1,end]),zeros(2,1),'k--')
            xlabel('time (s)')
            ylabel('heading (°)')
            
            subplot(1,4,4);
            polarhistogram(deg2rad(alpha(isnan(allprobR(:,2))==0)),deg2rad([0:10:360]))
            set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
            
            f2=figure('Position',[200 200 1400 300]);
            time=1/Fs:1/Fs:numel(rotSpeed_smooth)/Fs;
            subplot(1,4,[1 2 3]);hold on;
            hold on;plot(time,rotSpeed_smooth,'.','MarkerSize',12);hold on;
            plot(time(isnan(allprobR(2:end,2))==0),rotSpeed_smooth(isnan(allprobR(2:end,2))==0),'c.','MarkerSize',12);
            plot(time([1,end]),zeros(2,1),'k--')
            xlabel('time (s)')
            ylabel('angular velocity (°/s)')
            
            subplot(1,4,4);hold on;
            bar([1 2],turnCumSum/turnCumSum(2),'FaceColor','b')
            plot([0.5 2.5],[1 1],'k-.');
            set(gca, 'xtick',[1 2],'xticklabel',{'left','right'})
            xlim([0.5 2.5])
            ylim([0 max(turnCumSum/turnCumSum(2))*1.1])
            set(gca,'Layer','top')
            
            f3=figure;
            if size(heading_snips,1)>1
                boxplot(heading_snips(:,[1 3 5]),'labels',{'start','middle','end'});
                ylabel('mean vector length');
            end
            
            f4=figure;
            if size(heading_snips,1)>1
                
                headingAngle90=abs(heading_snips);
                for t=1:size(headingAngle90,2);headingAngle90(headingAngle90>90)=180-headingAngle90(headingAngle90>90);end
                boxplot(headingAngle90(:,[2 4 6]),'labels',{'start','middle','end'});
                ylabel('mean heading collapsed(0-90°)');
            end
            
            
            f5=figure('position',[300   500    1200    400]);
            subplot(1,3,1)
            polarhistogram(deg2rad(heading_snips(:,2)),deg2rad([0:10:360]))
            set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
            title('start')
            
            subplot(1,3,2)
            polarhistogram(deg2rad(heading_snips(:,4)),deg2rad([0:10:360]))
            set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
            title('middle')
            
            subplot(1,3,3)
            polarhistogram(deg2rad(heading_snips(:,6)),deg2rad([0:10:360]))
            set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
            title('end')
            
            f6=figure;
            if size(heading_snips,1)>1
                
                boxplot(rotation_snips,'labels',{'start','middle','end'});
                ylabel('avg rotational velocity (°/s)');
            end
            
            
            %% save data
            
%             pause
            indStart=strfind(filenames{i},'_T');
            save(['flightORI_',filenames{i}(indStart+1:end)],'alpha','rotSpeed_smooth','heading_snips','rotation_snips','turnCumSum','allprobR','patternR');
            close(f1,f2,f3,f4,f5,f6);
            
        end
        
    end
    
end


%% inbuilt functions

function [startInds, endInds]=checkDataPieces(data)
checkGaps=diff(isnan(data));
startInds=find(checkGaps==-1)+1;endInds=find(checkGaps==1)+1;
%if there is no nan in the beginning, but starts right with data,
%make the first start entry 1
if isnan(data(1))==0
    startInds=[1;startInds];
end
%if there is no nan in the end, make the last end length of
%data
if isnan(data(end))==0
    endInds=[endInds;numel(data)];
end

end