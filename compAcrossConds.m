% plots and statistically compares individually saved proboscis data from separate conditions 
% (data generated with plot_heat_map_proboscis_track and analyseProboscisMotion)

clear all
close all


ylabelText='rel. contacts';
animalInd=5;
dataInd=[1 2];
med=0;


animalOverlap=1;%put 1 if same animals in both datasets, 0 if not
med=0;%compare data against this value (i.e. 0 for position around line, 1 for ratios)
animalInd=1;
dataInd=1;
ylim0=0;

%% load data

%first condition
filename1=uigetfile('*.csv');
x=readtable(filename1);

%second condition
filename2=uigetfile('*.csv');
y=readtable(filename2);

%process data to turn into generalised analysis format
animalsX=x{:,animalInd};x(:,animalInd)=[];
animalsY=y{:,animalInd};y(:,animalInd)=[];

%for the half width of pattern contacts, remove animal IDs, these are just
%for general reference but will break logic of analysis here
if ~isempty(strfind(filename1,'_HW'))
x(:,1)=[];y(:,1)=[];
ylabelText='half-width probing area (mm)';
ylim0=1;
end

%check how many animals are the same in both datasets
% if size(x,2)==5 && size(y,2)==5
sameAnimalID = intersect(animalsX,animalsY);
if isempty(sameAnimalID);animalOverlap=0;
else animalOverlap=1;
end
% else
%     animalOverlap=0;
% end
%to force no overlap if we know there is none
% animalOverlap=0;

%compare datasets with shared animals
%detect animals that are the same in both datasets
if animalOverlap==1
comb=nan(length(sameAnimalID),2);

%compare ratio of contacts along pattern and the control
ind=[1 2];
for i=1:length(sameAnimalID)
    indX=animalsX==sameAnimalID(i);indY=animalsY==sameAnimalID(i);
comb(i,:)=[x{indX,ind(1)}./x{indX,ind(2)} y{indY,ind(1)}./y{indY,ind(2)}];
end

% comb=[comb(:,1); comb(:,2)];%reshape data for plotting
% grp=[ones(size(sameAnimalID));2*ones(size(sameAnimalID))];

elseif animalOverlap==0

    %compare ratio of contacts along pattern and the control when there are no
%shared animals in group
if dataInd==2
comb=[x{:,dataInd(1)}./x{:,dataInd(2)};y{:,dataInd(1)}./y{:,dataInd(2)}];
grp=[ones(size(x{:,dataInd(1)}));2*ones(size(y{:,dataInd(1)}))];

else
%non-paired data, no ratio to be calculated
comb=[x{:,dataInd};y{:,dataInd}];
grp=[ones(size(x{:,dataInd}));2*ones(size(y{:,dataInd}))];
end
%remove values that cannot be plotted / statistically analysed
% indnan=find(isnan(comb));comb(indnan)=[];grp(indnan)=[];

end



%% plot against each other
f1=figure('Position',[500 400 775/2 815/2]);hold on;
       
if animalOverlap==1

if size(comb,1)>1
    boxplot(comb);hold on;
    for u=1:size(comb,1)
        plot([1+0.01*randn(1,1),2+0.01*randn(1,1)],[comb(u,1),comb(u,2)],'k.-','Markersize',10);
    end
%     set(gca,'xticklabel','rel. contacts')
else
    plot([1 2],comb,'ro');xlim([0 3]);
    set(gca,'xtick',[1 2]);%,'xticklabel','rel. contacts')
end

[p,h,stats]=signrank(comb(:,1),comb(:,2));
testStatRank=stats.signedrank;
title(sprintf('n=%2.0f, p=%1.3f, rank=%1.3u',sum(isnan(comb(:,1))==0),p,stats.signedrank));
ylabel(ylabelText);
plot([0.5 2.5],[med med],'k--')
ylim([0 1.1*nanmax(comb(~isoutlier(comb(:))))])
xlim([0.5 2.5]);

%perform test against level
[p1,h,stats1]=signrank(comb(:,1),med);
[p2,h,stats2]=signrank(comb(:,2),med);
tbl1 = table({'cond1 vs cond2';['cond1 vs value ',num2str(med)];['cond2 vs value ',num2str(med)]},[p;p1;p2],[stats.zval;stats1.zval;stats2.zval],[stats.signedrank;stats1.signedrank;stats2.signedrank],'VariableNames',{'condition','p-value','z-val','signedrank'});

elseif animalOverlap==0

boxplot(comb,grp);hold on;
plot(1+0.01*randn(1,length(comb(grp==1))),comb(grp==1),'k.','Markersize',10);
plot(2+0.01*randn(1,length(comb(grp==2))),comb(grp==2),'k.','Markersize',10);
plot([0.5 2.5],[med med],'k--')
ylim([1.1*nanmin(comb(~isoutlier(comb(:)))) 1.1*nanmax(comb(~isoutlier(comb(:))))])
if ylim0==1
    ylim([0 1.1*nanmax(comb(~isoutlier(comb(:))))])
end
ylabel(ylabelText);
xlim([0.5 2.5]);

[p,h,stats]=ranksum(comb(grp==1),comb(grp==2))
testStatRank=stats.ranksum;

%perform test against level
%perform test against level
[p1,h,stats1]=signrank(comb(grp==1),med);
[p2,h,stats2]=signrank(comb(grp==2),med);
if ~isfield(stats1,'zval');stats1.zval=nan;end;if ~isfield(stats2,'zval');stats2.zval=nan;end;

tbl1 = table({'cond1 vs cond2';['cond1 vs value ',num2str(med)];['cond2 vs value ',num2str(med)]},[p;p1;p2],[stats.zval;stats1.zval;stats2.zval],[testStatRank;stats1.signedrank;stats2.signedrank],'VariableNames',{'condition','p-value','z-val','signedrank'});


end

figName=input('specify figure name');
print(f1,[figName,'.eps'],'-dpdf','-r300','-painters','-bestfit')

disp(figName)
disp(['signrank test against ',num2str(med)])
disp(tbl1)
writetable(tbl1, [figName,'_stats.txt'])

%% correlate x and y
%i.e. proboscis position rel. to body vs turn index, or vs prob pos rel to
%pattern

%load data

%first condition
filename1=uigetfile('*.csv');
x=readtable(filename1);

%second condition
filename2=uigetfile('*.csv');
y=readtable(filename2);


f1=figure('Position',[500 400 775/2 815/2]);hold on;
plot(x{:,2},x{:,3},'.b','MarkerSize',18)
plot(y{:,2},y{:,3},'.c','MarkerSize',18)
%for turn index
% set(gca,'YScale','log');
% ylim([0.3 3]);xlim([-6 6])
% plot([0 0],[0.3 3],'-k')
% plot([-6 6],[1 1],'-k')

% ylabel('turn index')
%for pattern position
% set(gca,'YScale','log');
ylim([-5 5]);xlim([-6 6])
ylabel('position (pattern is 0 mm)')
plot([0 0],[-5 5],'-k')
plot([-6 6],[0 0],'-k')

axis square
xlabel('prob position')


[rho,pval] = corr([x{:,2};y{:,2}],[x{:,3};y{:,3}]);
title(sprintf('n=%2.0f, p=%1.3f, rho=%1.3u',sum(isnan([x{:,2};y{:,2}])==0),pval,rho));

figure;histogram([x{:,2};y{:,2}],15);xlim([-6 6])
figure;histogram([x{:,3};y{:,3}],12);xlim([-5 5]);

figure;histogram(log([x{:,3};y{:,3}]),12);xlim(log([0.3 3]));