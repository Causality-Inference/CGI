clear;
clc;
load found_Genes_Leukemia.mat
load normalized_Leukemia.mat
% genes = [3426,6753,4324,4153]; %V-strucure Cai
len = length(genes);
d = d(:,[genes,7130]);
d = sortrows(d,len+1);
M = zeros(2,len);
for i = 1:len
    M(1,i) = sum(d(1:25,i));
    M(2,i) = sum(d(26:72,i));
end

[a,b] = sortrows(M',-1);
c = d;
c(:,len+1) = [];
c = c(:,b);
c = [c,d(:,len+1)];
d = c;

for i = 1:72
    for j = 1:len %(len+1)
        if d(i,j) > 0 && d(i,j) <= 0.5
            plot(i,j,'s','MarkerFace',[255,150,255]/256,'MarkerSize',10)
            hold on
        end
        if d(i,j) > 0.5 && d(i,j) <= 1
            plot(i,j,'s','MarkerFace',[255,100,200]/256,'MarkerSize',10)
            hold on
        end
        if d(i,j) > 1 && d(i,j) <= 1.5
            plot(i,j,'s','MarkerFace',[255,75,150]/256,'MarkerSize',10)
            hold on
        end
        if d(i,j) > 1.5 && d(i,j) <= 2
            plot(i,j,'s','MarkerFace',[255,50,100]/256,'MarkerSize',10)
            hold on
        end
        if d(i,j) > 2 && d(i,j) <= 2.5
            plot(i,j,'s','MarkerFace',[220,25,50]/256,'MarkerSize',10)
            hold on
        end
        if d(i,j) > 2.5 
            plot(i,j,'s','MarkerFace',[255,0,0]/256,'MarkerSize',10)
            hold on
        end
%-------------------------------------------------
        if d(i,j) <= 0 && d(i,j) > -0.5
            plot(i,j,'s','MarkerFace',[176,196,222]/256,'MarkerSize',10)
            hold on
        end
        if d(i,j) <= -0.5 && d(i,j) > -1
            plot(i,j,'s','MarkerFace',[100,149,237]/256,'MarkerSize',10)
            hold on
        end
        if d(i,j) <= -1 && d(i,j) > -1.5
            plot(i,j,'s','MarkerFace',[65,105,225]/256,'MarkerSize',10)
            hold on
        end
        if d(i,j) <= -1.5 && d(i,j) > -2
            plot(i,j,'s','MarkerFace',[0,0,128]/256,'MarkerSize',10)
            hold on
        end
        if d(i,j) <= -2 && d(i,j) > -2.5
            plot(i,j,'s','MarkerFace',[0,0,139]/256,'MarkerSize',10)
            hold on
        end
        if d(i,j) <= -2.5 
            plot(i,j,'s','MarkerFace',[25,25,112]/256,'MarkerSize',10)
            hold on
        end
    end
end
plot([47.40,47.40],[0,40],':g','LineWidth',2)
hold on
plot([0,80],[5.40,5.40],':g','LineWidth',2)
hold on
xlabel('Sample size');
ylabel('Features');