%this script reads the tracking points from DeepLabCut, which are packaged
%in an .h5 file, and extracts the x/y positions of the points, as well as
%their labels, to feed them into a .csv file which we can further process
%with the Hedrick's Lab's dltv tracking software for post processing
clear all
close all

%% load data
[filenames, path]=uigetfile('*.h5','multiselect','on');
%if only one selected turn filenames into cell
if iscell(filenames)==0
    temp=filenames;clear filenames;
    filenames=cell(1,1);filenames{1}=temp;
end


%% read data
for f=1:numel(filenames)
    h5disp(filenames{f});
data = h5read(filenames{f},'/df_with_missing/table');
numData=data.values_block_0'; %this is the numerical values (x/y tracking points plus likelyhood estimates
xypoints=numData;
likelihood=numData(:,3:3:end);
xypoints(:,3:3:end)=[];%remove likelihood to only get xy points

%set all points below likelihood 0.85 to nan
for u=1:size(xypoints,2)/2
  xypoints(likelihood(:,u)<0.9,2*u-1:2*u)=nan;
end

%sometimes the proboscis gets wrongly detected when the animal isn't
%there. and the likelyhood is still very high. for the other labels it
%seems more robustly rejecting wrong tracks. Therefore I would suggest to
%say that if the combined likelyhood for all bodyparts is low, reject the
%proboscis as well. Then likely the animal is not in the frame
sumLikelihood=nanmean(likelihood(:,2:end),2);
xypoints(sumLikelihood<0.75,:)=nan;

%also, sometimes the algorithm defines a point in the flower (worst case)
%or somewhere outside as the "alteranative point" for the proboscis. then
%the tracking keeps jumping there. this is problematic if the likelihood
%there is high! then rather delete those points. need to have high
%repetition rate to not delete "real" tracking points

% for now, can battle this by deleting all points that are not close to the head and thorax point
radCircle=120;
for r=1:size(xypoints,1)
    if isnan(xypoints(r,1))==0
    th = 0:pi/50:2*pi;
    xunitC = radCircle * cos(th) + xypoints(r,5);%use the head as the circle centre
    yunitC = radCircle * sin(th) + xypoints(r,6);
%     plot(xunitC,yunitC,'-');
    
    in = inpolygon(xypoints(r,1),xypoints(r,2),xunitC,yunitC);
    if in==0
    xypoints(r,1:2)=[nan nan];
    end
    end
end

%% !!!ATTENTION!!!!
% somehow one line seems to be missing from the h5 file, and it seems to be the first one.
% Thus, in order to make things synchronous with the video again, add this line
xypoints=[nan(1,size(xypoints,2));xypoints];
%%

% mirror the y-axis. somehow this comes out mirrored from the DLC
% algorithm. 
ysize=input('what is the video height in pixels?');
xypoints(:,2:2:end)=abs(ysize-xypoints(:,2:2:end));



%% get right labels
%this is the labels that we should find, in the order that they were tracked 
% this is not the most elegant solution, but it should work

attval = h5readatt(filenames{f},'/df_with_missing','non_index_axes');
cellAttributes=strsplit(attval);

lookForLabels={'Proboscis','Thorax','Head','AntennaR','AntennaL','Abdomen'};
labels=cell(length(lookForLabels),1);
count=0;
for i=1:length(cellAttributes)
    for j=1:length(lookForLabels)
lookInd(j)=isempty(strfind(cellAttributes{i},['V',lookForLabels{j}]))==0; %i dont know why the algorithm labels it all as V"Label" but lets go with it
    end
if sum(lookInd)~=0
count=count+1;
labels{count}=lookForLabels{lookInd};
end
end

%% save the data as .csv
labels_xy=cell(length(labels)*2,1);
for u=1:length(labels)
    labels_xy{u*2-1}=[labels{u},'_X'];
    labels_xy{u*2}=[labels{u},'_Y'];
end

indEnd=strfind(filenames{f},'DeepCut')+6;%find the end of DeepCut string and end filename there

%check if a file with this name already exists, to make sure not to
%overwrite it (this is especially critical when there are xcel files which
%have already been curated under the same name. because then all the
%curated data will be lost

if exist([filenames{f}(1:indEnd),'_xypts.csv'])
error('A file with this name exists already and will be overwritten when saving this data. Please make sure to resolve the naming conflict');

else
data_table = array2table(xypoints,'VariableNames',labels_xy);
writetable(data_table,[filenames{f}(1:indEnd),'_xypts.csv'],'Delimiter',',');

%need to also save tables for DLTV to work properly with the other
%properties
offset_table=table(zeros(size(xypoints,1),1),'VariableNames',{'cam1_offset'});
writetable(offset_table,[filenames{f}(1:indEnd),'_offsets.csv'],'Delimiter',',');

res_table=table(nan(size(xypoints,1),size(xypoints,2)/2));
writetable(res_table,[filenames{f}(1:indEnd),'_xyzres.csv'],'Delimiter',',');

xzy_table=table(nan(size(xypoints,1),size(xypoints,2)/2*3));
writetable(xzy_table,[filenames{f}(1:indEnd),'_xyzpts.csv'],'Delimiter',',');
end

end

