clc;
clear;
% load normalized_Leukemia
% load found_Genes_Leukemia
% load normalized_Prostate
% load found_Genes_Prostate
% load normalized_Leukemia_ATL
% load found_Genes_ATL
% load normalized_Liver
% load found_Genes_Liver
load normalized_Colorectal
load found_Genes_Colorectal
n = size(d,2)-1;
c = d(:,n+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(0,0,'dk')
hold on
x = d(:,genes(1));
p = 1-abs(corr(c,x,'type','Spearman'));
x = (1 - 2*rand(1,1))*p;
y = ((p^2-x.^2).^(1/2)).*randsrc(1,1);
plot(x,y,'or')
hold on;
x = d(:,1);
p = 1-abs(corr(c,x,'type','Spearman'));
x = (1 - 2*rand(1,1))*p;
y = ((p^2-x.^2).^(1/2)).*randsrc(1,1);
plot(x,y,'sg')
hold on;
legend('Disease variable','Discovered genes','Other genes');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
genes(1) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n
    x = d(:,i);
    p = 1-abs(corr(c,x,'type','Spearman'));
    x = (1 - 2*rand(1,1))*p;
    y = ((p^2-x.^2).^(1/2)).*randsrc(1,1);
    if isempty(intersect(i,genes))
        plot(x,y,'sg')
        hold on;
    end
end
for i = 1:n
    x = d(:,i);
    p = 1-abs(corr(c,x,'type','Spearman'));
    x = (1 - 2*rand(1,1))*p;
    y = ((p^2-x.^2).^(1/2)).*randsrc(1,1);
    if ~isempty(intersect(i,genes))
        plot(x,y,'or')
        hold on;
    end
end

xlabel('1-|\rho_{xy}|');
ylabel('1-|\rho_{xy}|');

legend('Disease variable','Discovered genes','Other genes');

set(gca,'xtick',[-1:1:1]);
set(gca,'xTickLabel',{0,1,2});
set(gca,'ytick',[-1:1:1]);
set(gca,'yTickLabel',{0,1,2});