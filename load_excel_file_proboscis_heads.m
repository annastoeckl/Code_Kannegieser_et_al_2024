function load_excel_file_proboscis_heads(files)
% loads -csv files with the animal data in the following format: X, Y
% coordinate of the following points in this order:
% proboscis, head, thorax, proboscisbase, left antenna, right antenna, outline of the flower, pattern 1, pattern 2
% outline of the flower are 4 points, marking orthogonal points on the
% circular outline of the flower, to reconstruct its circumference
% pattern 1 can be a line or a cross or a circle, pattern 2 is a pattern
% that is to be subtracted from pattern 1
% output: mat file with the extracted data with parameter names to be used
% by append_and_rotate_data.m

if exist('files','var')==0
    %load filename with xypts extension
    files=uigetfile('*xypts.csv');
end

% [num] =  csvread(files,1);
%new way of reading table, allows to read headers too
Tbl=readtable(files);
try %this is the "old" (Matlab 2017) way to read the table
num=Tbl{:,:};
catch %with the new way, some columns that contain NaN are read as text, 
    %use readmatrix instead
    %or use readtable without the spreadsheet bit Tbl=readtable(files,'FileType','spreadsheet');
num=readmatrix(files);   
end
if iscell(num) 
Tbl=readtable(files,'FileType','spreadsheet');
num=Tbl{:,:};
end
labels = Tbl.Properties.VariableNames;

%select outline of the flower (cross through its center on both x and y dimension).
% it will always be 4 x/y points, and the first entry
%of few digits right after the main tracking data
indOutline=find(sum(isnan(num)==0)==4,1,'first');%take the first entry, because pattern might also have 4 digits
outline=num(:,indOutline:indOutline+1);
outline(isnan(outline(:,1)),:)=[];

%if there is a pattern, it will be indicated in the next columns
if size(num,2)>indOutline+1
    pattern=num(:,indOutline+2:indOutline+3);
    indnan=isnan(pattern);
    pattern(indnan(:,1),:)=[];
else
    pattern=[nan nan];
end


%if there is a secondary pattern, it will be indicated in the columns after
%the primary pattern
if size(num,2)>indOutline+3
    pattern2=num(:,indOutline+4:indOutline+5);
    indnan=isnan(pattern2);
    pattern2(indnan(:,1),:)=[];
else
    pattern2=[nan nan];
end

%select moth tracks
%delete columns with orientation and pattern, and then the rest should fall
%in the right order
if sum(isnan(pattern(:)))~=0 %if there is no pattern
    num(:,indOutline:indOutline+1)=[];
    labels(indOutline:indOutline+1)=[];
elseif sum(isnan(pattern(:)))==0 && sum(isnan(pattern2(:)))~=0  %if there is a pattern but no secondary one
    num(:,indOutline:indOutline+3)=[];
    labels(indOutline:indOutline+3)=[];
else %if there is both a primary and secondary pattern
    num(:,indOutline:indOutline+5)=[];
    labels(indOutline:indOutline+5)=[];
end

prob=num(:,1:2);%proboscis track
if size(num,2)>2
    thorax=num(:,3:4);%thorax track
else thorax=[nan nan];
end
if size(num,2)>4
    head=num(:,5:6);%head track
else head=[nan nan];
end
if size(num,2)>6
    probBase=num(:,7:8);%head prob base
else probBase=[nan nan];
end
if size(num,2)>8
    antBaseL=num(:,9:10);%head left antenna base
else antBaseL=[nan nan];
end
if size(num,2)>10
    antBaseR=num(:,11:12);%head right antenna base base
else antBaseR=[nan nan];
end

dataMatrix=num;%save the tracking indices here, to prepare for REWRITE of code later so that all
%indices can be kept in matrix

%swap y-axis in variables, because the DLC variables (prob, head, thorax)
%are already swapped
probBase(:,2)=500-probBase(:,2);
antBaseL(:,2)=500-antBaseL(:,2);
antBaseR(:,2)=500-antBaseR(:,2);

figure;hold on;%plot(num(:,2),'b');
plot(prob(:,1),prob(:,2),'c*');hold on;
plot([pattern(:,1);pattern(1,1)],[pattern(:,2);pattern(1,2)],'k-');
plot([pattern2(:,1);pattern2(1,1)],[pattern2(:,2);pattern2(1,2)],'-','Color',[0.5 0.5 0.5]);
plot(thorax(:,1),thorax(:,2),'bo');hold on;
plot(head(:,1),head(:,2),'go');hold on;
plot(probBase(:,1),probBase(:,2),'r.');hold on;
plot(antBaseL(:,1),antBaseL(:,2),'m.');hold on;
plot(antBaseR(:,1),antBaseL(:,2),'m.');hold on;

%plot flower circle
temp=diff(outline);temp(2,:)=[];
midpoint=[max(abs(temp(:,1)))/2+min(outline(:,1))  max(abs(temp(:,2)))/2+min(outline(:,2))];
lengthCircle=nanmean(nanmax(abs(temp)),2);

th = 0:pi/50:2*pi;
xunitC = lengthCircle/2 * cos(th) + midpoint(1);
yunitC = lengthCircle/2 * sin(th) + midpoint(2);
plot(xunitC,yunitC,'-');


ylabel('px')
xlabel('px')


disp('number of contacts:')
disp(sum(isnan(prob(:,1))==0));

outputname=['data_',files(1:end-4),'.mat'];
save(outputname,'pattern2','pattern','prob','outline','thorax','head','probBase','antBaseL','antBaseR','dataMatrix','labels');


end