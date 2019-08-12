function [u,B,v] = LBD(A,v0,iter,tr)
if tr == 1
    transp = 'notransp';
    notransp = 'transp';
else
    transp = 'transp';
    notransp = 'notransp';
end
BlS = size(v0,2);
[v(:,1:BlS),~] = qr(v0,0);
u(:,1:BlS) = A(v(:,1:BlS),notransp);
[u(:,1:BlS), a{1}] = qr(u(:,1:BlS),0);
v(end,iter*BlS) = 0;
u(end,iter*BlS) = 0;
for i = 1:iter-1
    v(:,BlS*i+1:(i+1)*BlS) = A(u(:,(i-1)*BlS+1:BlS*i),transp) - v(:,(i-1)*BlS+1:BlS*i)*a{i};
    [v(:,BlS*i+1:(i+1)*BlS),x] = cgs(v(:,1:BlS*i),v(:,BlS*i+1:(i+1)*BlS));
    b{i} = x(end-BlS+1:end,end-BlS+1:end)';
    
    u(:,BlS*i+1:(i+1)*BlS) = A(v(:,BlS*i+1:(i+1)*BlS),notransp) - u(:,(i-1)*BlS+1:BlS*i)*b{i};
    [u(:,BlS*i+1:(i+1)*BlS),y] = cgs(u(:,1:BlS*i),u(:,BlS*i+1:(i+1)*BlS));
    a{i+1} = y(end-BlS+1:end,end-BlS+1:end);
end
B = blkdiag(a{:})+[zeros(BlS*(iter-1),BlS), blkdiag(b{:}); zeros(BlS,iter*BlS)];

end