function [reset,newmvs] = checkReset(rc,rcf,run,A,u,s,v,notransp)
%Check for reset
reset = 0;
newmvs = 0;
if any(rc*rcf > run)
    rvn = norm(A(v,notransp)-u*diag(s)); newmvs = size(u,2);
    rcf = rvn/(rc) * 1.22; %Modify rc to match the actual rvn
    if (rc*rcf > run)
        reset = 1;
    end
end