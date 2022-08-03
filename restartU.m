function [Qtilde,R] = restartU(u,s,v,index,minRS,yold)

index1 = index(1:minRS);
index2 = index(minRS+1:end);
numOld = size(yold,2);
y2 = v(:,index2);
w1 = u(:,index1); w2 = u(:,index2);
s1 = s(index1); s2 = s(index2);

if numOld ~= 0
    [q,r] = qr(diag(s2)*(y2'*yold),0);
    Qtilde = [w1 w2*q];
    R(1:minRS+numOld,1:minRS+numOld) = blkdiag(diag(s1),r);
else
    R(1:minRS,1:minRS) = diag(s1);
    Qtilde = w1;
end
end