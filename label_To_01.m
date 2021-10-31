clear;
clc;
load Lung.mat
l = length(labels);
c = ones(1,l);

for i = 1:l
    q =  strfind(char(labels(i)),'ATL');
    if ~isempty(q)
        c(i) = 0;
    end
end

d = [data;c];
d = d';

for i = 1:(size(d,2)-1)
    a = d(:,i);
    a = (a-mean(a))/std(a);
    d(:,i) = a;
end