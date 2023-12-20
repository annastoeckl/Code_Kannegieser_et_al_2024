% plots and statistically compares individually saved data of halfwidth of proboscis contact distributions from separate conditions 
% (data generated with plot_heat_map_proboscis_track)

clear all
close all

ylabelText='half-width (mm)';
animalInd=1;
dataInd=2;
med=0;
numelConds=3;
ylim0=0;

%% load data
comb=[];grp=[];

for i=1:numelConds
%first condition
filename=uigetfile('*.csv');
x=readtable(filename);
legend{i}=filename;
% 
% %second condition
% filename2=uigetfile('*.csv');
% y=readtable(filename2);
% 
% %third condition
% filename2=uigetfile('*.csv');
% y=readtable(filename2);



%generate data structure for stats
comb=[comb;x{:,dataInd}];
grp=[grp;i*ones(size(x{:,dataInd}))];
end
% %check how many animals are the same in both datasets
% % if size(x,2)==5 && size(y,2)==5
% sameAnimalID = intersect(animalsX,animalsY);
% if isempty(sameAnimalID);animalOverlap=0;
% else animalOverlap=1;
% end

%% plot against each other if paired data
f1=figure('Position',[500 400 775/2 815/2]);hold on;
       
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


%% plot unpaired data
f1=figure('Position',[500 400 775/2 815/2]);hold on;

boxplot(comb,grp);hold on;
for i=1:numelConds
plot(i+0.01*randn(1,length(comb(grp==i))),comb(grp==i),'k.','Markersize',10);
end
plot([0.5 numelConds+.5],[med med],'k--')
ylim([1.1*nanmin(comb(~isoutlier(comb(:)))) 1.1*nanmax(comb(~isoutlier(comb(:))))])
if ylim0==1
    ylim([0 1.1*nanmax(comb(~isoutlier(comb(:))))])
end
ylabel(ylabelText);
xlim([0.5 numelConds+.5]);
ylim([0 20])

%stats
[p,h,stats]=signrank(comb(grp==1),comb(grp==3))

[p,tbl,stats] = anova1(comb,grp)

[p,h,stats]=ranksum(comb(grp==1),comb(grp==3))
testStatRank=stats.ranksum;

%perform test against level
%perform test against level
[p1,h,stats1]=signrank(comb(grp==1),med);
[p2,h,stats2]=signrank(comb(grp==2),med);

if ~isfield(stats1,'zval');stats1.zval=nan;end;if ~isfield(stats2,'zval');stats2.zval=nan;end;

tbl1 = table({'cond1 vs cond2';['cond1 vs value ',num2str(med)];['cond2 vs value ',num2str(med)]},[p;p1;p2],[stats.zval;stats1.zval;stats2.zval],[testStatRank;stats1.signedrank;stats2.signedrank],'VariableNames',{'condition','p-value','z-val','signedrank'});




figName=input('specify figure name');
print(f1,[figName,'.eps'],'-dpdf','-r300','-painters','-bestfit')

disp(figName)
disp(['signrank test against ',num2str(med)])
disp(tbl1)
writetable(tbl1, [figName,'_stats.txt'])

