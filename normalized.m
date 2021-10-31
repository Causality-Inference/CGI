clc;
clear;
load normalized_Colorectal
% d = ColorectalGSE25070;
for i = 1:size(d,2)-1
    a = d(:,i);
    a = (a-mean(a))/std(a);
    d(:,i) = a;
end
