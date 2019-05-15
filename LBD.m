function [u,B,v,atu] = LBD(A,v0,iter,tr)
if tr == 1
    transp = 'notransp';
    notransp = 'transp';
else
    transp = 'transp';
    notransp = 'notransp';
end
v(:,1) = v0/norm(v0);
u(:,1) = A(v(:,1),notransp);
a(1) = norm(u(:,1));
v(end,iter) = 0;
u(:,1) = u(:,1)/a(1);
u(end,iter) = 0;
atu = zeros(size(v));
for i = 2:iter
    atu(:,i-1) = A(u(:,i-1),transp);
    v(:,i) = atu(:,i-1) - a(i-1)*v(:,i-1);
    [v(:,i),x] = cgs(v(:,1:i-1),v(:,i));
    b(i-1) = x(end);
    
    u(:,i) = A(v(:,i),notransp) - b(i-1)*u(:,i-1);
    [u(:,i),y] = cgs(u(:,1:i-1),u(:,i));
    a(i) = y(end);
end
atu(:,i) = A(u(:,i),transp);
B = diag(a)+[zeros(iter-1,1), diag(b); zeros(1,iter)];

end