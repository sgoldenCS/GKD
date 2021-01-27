function [thickright,yold] = updateV(v,index,minRS,vrold)

y1 = v(:,index(1:minRS));
numOld = size(vrold,2);
yold =[];
if numOld ~= 0
    yold = cgs(y1, [vrold; zeros(size(y1,1)-size(vrold,1),numOld)]);
    thickright = [y1 yold];
else
    thickright = y1;
end