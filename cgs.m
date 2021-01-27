function [x,r] = cgs(B,x)
k = size(B,2);
nv = size(x,2);
r = zeros(k+nv,nv);

if isempty(B)
    [x,r] = qr(x,0);
    return;
end

for i = 1:nv
    iter = 0;
    normx = 0;
    prev = 1;
    while normx <= 0.75*prev && iter < 4
        iter = iter + 1;
        prev = norm(x(:,i));
        coefB = B'*x(:,i);
        coefX = x(:,1:i-1)'*x(:,i);
        x(:,i) = x(:,i) - B*coefB - x(:,1:i-1)*coefX;
        r(1:k+i-1,i) = r(1:k+i-1,i) + [coefB; coefX];
        normx = norm(x(:,i));
    end
    if iter > 3
        warning('Vector not orthogonalized!');
    end
    r(k+i,i) = normx;
    x(:,i) = x(:,i)/normx;
end

end