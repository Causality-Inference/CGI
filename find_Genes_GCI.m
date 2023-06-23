clear;
clc;
addpath(genpath('./algo'))
% load normalized_Prostate.mat
load normalized_Leukemia.mat
data = d;
n = size(data,2)-1;
x = data(:,n+1);
data = data(:,1:n);
non = []; %non-adjacent
nontry1 = []; %non-adjacent
nontry2 = []; %non-adjacent
cov='covSEiso';
Ncg=100;
hyp=[4;log(4);log(sqrt(0.01))];
tempRes = []; % save res x-Z
tempZ = []; % save Z
%---------------------------- 0-order CI tests-----------------------------
fprintf('--------------- 0-order CI tests \n');
for i = 1:n
    ind1 = PaCoTest(x,data(:,i),[],0.05);
    if ind1 
        ind2 = KCIT(x,data(:,i),[],[]);
        if ind2
            non = [non,i];
        end
    end
end
%---------------------------- 1-order CI tests-----------------------------
fprintf('--------------- 1-order CI tests \n');
idx1 = setdiff(1:n,non);
len1 = length(idx1);
for j = 1: len1
    [j,len1]
    for k = 1: len1
        if j~=k && isempty(intersect(idx1(k),non))
            y = data(:,idx1(j));
            z = data(:,idx1(k));
            ind1 = PaCoTest(x,y,z,0.05);
            if ind1
                try
                    xf = fit_gpr(z,x,cov,hyp,Ncg);
                    res1 = xf-x;
                    yf = fit_gpr(z,y,cov,hyp,Ncg);
                    res2 = yf-y;
                    ind2 = KCIT(res1, res2,[],[]);
                    if ind2
                        tempRes = [tempRes,res1];
                        tempZ = [tempZ;idx1(k),0];
                        non = [non,idx1(j)];
                        break;
                    end
                catch
                    non = [non,idx1(j)];
                    break;
                end
            end
        end
    end
end
%--------------------------- 2-order CI tests------------------------------
fprintf('--------------- 2-order CI tests \n');
idx2 = setdiff(1:n,non);
len2 = length(idx2);
M = nchoosek(1:len2,2);
m = 1:size(M,1);
for p = 1: len2
    j = idx2(p);
    for k = 1:m
        y = data(:,j);
        z = data(:,idx2(M(k,:)));
        if idx2(M(k,1))~= j && idx2(M(k,2))~= j
            ind = PaCoTest(data(:,i), data(:,j),z,0.05);
            if ind
                try
                    xf = fit_gpr(z,x,cov,hyp,Ncg);
                    res1 = xf-x;
                    yf = fit_gpr(z,y,cov,hyp,Ncg);
                    res2 = yf-y;
                    ind2 = KCIT(res1, res2,[],[]);
                    if ind2
                        tempRes = [tempRes,res1];
                        tempZ = [tempZ;idx2(M(k,:))];
                        non = [non,j];
                    end
                catch
                    non = [non,j];
                    break;
                end
            end
        end
    end
end
%--------------- find genes by regression 1st ---------------------------------
fprintf('--------------- find genes by regression 1st \n');
pa = setdiff(1:n,non);
l = length(pa);
found_Genes_1st = [];
[~,idz1] = intersect(tempZ(:,1),pa);
[~,idz2] = intersect(tempZ(:,2),pa);
for k = 1:length(idz1)
    if KCIT(data(:,tempZ(idz1(k),1)), tempRes(:,idz1(k)),[],[])
        found_Genes_1st = [found_Genes_1st,idz1(k)];
    end
end
for k = 1:length(idz2)
    if KCIT(data(:,tempZ(idz2(k),1)), tempRes(:,idz2(k)),[],[])
        found_Genes_1st = [found_Genes_1st,idz2(k)];
    end
end
%--------------- find genes by regression 2nd--------------------------------
fprintf('--------------- find genes by regression 2nd \n');
found_Genes_2nd = [];
M = nchoosek(1:l,2);
lenM = size(M,1);
for p = 1: lenM
    [p,lenM]
    j = pa(M(p,:));
    z = data(:,j);
    try
        xf = fit_gpr(z,x,cov,hyp,Ncg);
        res1 = xf-x;
        ind3 = KCIT(res1, z(:,1),[],[]);
        ind4 = KCIT(res1, z(:,2),[],[]);
        if ind3 
            found_Genes_2nd = [found_Genes_2nd,j(1)];
        end
        if ind3
            found_Genes_2nd = [found_Genes_2nd,j(2)];
        end
    catch
    end
end
%--------------- results--------------------------------
found_Genes = unique([found_Genes_1st,found_Genes_2nd])
