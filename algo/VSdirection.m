function[Dskeleton] = VSdirection(Cskeleton,AM,BM)
Dskeleton = zeros(size(Cskeleton,1),size(Cskeleton,1));
for p = 1:size(AM,2)  %%%%% AM = zeros(2,Lmn);  存在很多0的情况
    x = AM(1,p);
    y = AM(2,p);
    
    if x == 0    
        break;
    end

    pax = find(Cskeleton(:,x)==1)';
    chx = find(Cskeleton(x,:)==1);
    pcx = unique([pax,chx]);
    pay = find(Cskeleton(:,y)==1)';
    chy = find(Cskeleton(y,:)==1);
    pcy = unique([pay,chy]);
    interPc = intersect(pcx,pcy);
    if ~isempty(interPc)
        if  BM(1,p) == 0
            Dskeleton(x,interPc)=1;
            Dskeleton(y,interPc)=1;
            Dskeleton(interPc,x)=0;
            Dskeleton(interPc,y)=0;
        end
        if  BM(1,p) > 0
            diffPc = setdiff(interPc,BM(:,p)');
            if ~isempty(diffPc)
                Dskeleton(x,diffPc)=1;
                Dskeleton(y,diffPc)=1;
                Dskeleton(diffPc,x)=0;
                Dskeleton(diffPc,y)=0;
            end
        end
    end
end

end