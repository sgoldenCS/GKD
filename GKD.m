function [varargout] = GKD(A,varargin)
%% GKD   Find a few singular values and vectors for large sparse matrices
%  
% Outputs:
% 
% [s] = GKD(...)
% [s,r] = GKD(...)
% [u,s,v] = GKD(...)
% [u,s,v,r] = GKD(...)
% [u,s,v,r,stats] = GKD(...)
% [u,s,v,r,stats,hist] = GKD(...)
%
% u             Left Singular Vectors
% s             Singular Values
% v             Right Singular Vectors
% r             norm(A'*u - s*v) for each singular triplet
% stats         Matvecs, Time and estimated norm(A)
% hist          Convergence History
%
% Inputs:
%
% [...] = GKD(A,numValues,target)
% [...] = GKD(A,numValues,target,opts)
% [...] = GKD(A,numValues,target,opts,Pata)
%
% A             m by n matrix
% numValues     number of singular values to return
% target        seek singular value nearest target
% opts          structure containing extra solver options
% Pata          preconditioner for A'*A
%
% Options for GKD
%   opts.tol                residual norm tolerance
%   opts.maxBasis           max number of basis vectors in V, U
%   opts.aNorm              norm(A) estimate
%   opts.maxMV              maximum number of matrix-vector multiplications
%   opts.v0                 initial vector for V
%   opts.disp               options for printing history to console
%   opts.minRS              number of vectors to maintain after restart
%   opts.numPk              +k criteria
%   opts.maxII              maximum number of inner solver iterations
%   opts.seed               Sets seed for the random number generator
%   opts.locking            Turns hard locking on if set
%   opts.isdouble           0 = Single Precision, 1 = Double Precision
%   opts.LBD                1 = Start with LBD basis up to maxBasis-1
%   opts.AtQ                1 = Keep extra storage for AtQ. Reduces Matvecs
%                                   with locking = 0

inputs = varargin;
[m,n,transp,notransp,Transpose,numValues,target,opts,Pata] = varginParse(A,inputs);


[ERROR,LBDflag,xtrStorage,normA,tol,maxBasis,maxMVs,v0,display,minRS,numOld,maxII,HLock] = optsParse(A,m,n,opts);
if ERROR ~= 0
    for i = nargout:-1:1
        varargout{i} = [];
    end
    return;
end

if HLock == 0 && minRS < numValues+2
    error('minRestartSize too small for soft-locking\n\tIncrease minRS > %d or locking = 1',numValues+2);
end

if isfield(opts,'isdouble') && opts.isdouble == 0
    if ~isa(A,'function_handle')
        A = @(x,tr) singleA(A,x,tr);
    end
    macheps = eps('single');
    V = single(zeros(n,maxBasis));
    Q = single(zeros(m,maxBasis));
    R = single(zeros(maxBasis,maxBasis));
    if xtrStorage
        AtQ = single(zeros(n,maxBasis));
    end
else
    if ~isa(A,'function_handle')
        A = @(x,tr) doubleA(A,x,tr);
    end
    V = zeros(n,maxBasis);
    Q = zeros(m,maxBasis);
    R = zeros(maxBasis,maxBasis);
    if xtrStorage
        AtQ = zeros(n,maxBasis);
    end
    macheps = eps;
end

%If tol is too small, it may not be possible to converge
if tol < 5*normA*macheps
    error("Tolerance too small. opts.tol must be larger than %e",string(5*macheps));
end

reset = 0;
restarts = 0;
mvs = 0;     %MV Counter
k = 1;       %Current basis size
outerits = 0;
t = 1;
touch = 1;
sigma_prev = 0;
sigma = 0;
rcf = 1; %reset criteria factor
hist = [];
found = 0;

%Storage for returned SVD and residuals
LV = [];
LQ = [];
LSig = [];
LR = [];

tic
if LBDflag
    mb = maxBasis-1;
    [Q(:,1:mb), R(1:mb,1:mb), V(:,1:mb), AtQ(:,1:mb)] = LBD(A,v0,mb,Transpose);
    k = mb;
    mvs = mvs + 2*k;
else
    v0=v0/norm(v0);
    V(:,1:k) = v0;
    Q(:,1:k) = A(v0,notransp); mvs = mvs + 1;
    R(1,1) = norm(Q(:,1)); Q(:,1) = Q(:,1)/R(1,1);   % QR of AV
end

%% Main iteration
while (maxMVs <= 0 || mvs<maxMVs) && found < numValues
    tempQ = Q(:,1:k);
    tempV = V(:,1:k);
    if xtrStorage
        tempAtQ = AtQ(:,1:k);
    end
    cflag = 0; %Convergence flag
    
    % Extraction process
    [ur,sr,vr]=svd(R(1:k,1:k));
    [~,index] = sort(abs(target - diag(sr)));
    sigma_prev = sigma;
    sigma = sr(index(t),index(t));
    u = tempQ*ur(:,index(t));
    v = tempV*vr(:,index(t));
    vrold = vr(:,index(t:(min(size(vr,2),numOld))+t-1));
    
    %Calculate left residual
    %ru = A(u,transp) - sigma*v; mvs = mvs + 1;
    if xtrStorage
        ru = tempAtQ*ur(:,index(t)) - sigma*v;
    else
        ru = A(u,transp) - sigma*v; mvs = mvs+1;
    end
    run = norm(ru);
    
    outerits = outerits+1;
    
    %%% DISPLAY %%%
    if display == 3 %Requires an extra matvec to compute rvn
        rv = A(v,notransp)-sigma*u;
        rvn = norm(rv);
        resetcriteria = 4*normA*macheps*sqrt(restarts)*rcf;
        orthv = norm(tempV'*tempV - eye(k));
        orthq = norm(tempQ'*tempQ - eye(k));
        fprintf('iter %4d mv %4d sig %7.3e Ru %7.3e Rv %7.3e conv %3d RCrit %7.3e ',...
            [outerits mvs sigma run rvn found resetcriteria]);
        fprintf('OrthV %7.3e OrthU %7.3e Restarts %4d\n', [orthv orthq restarts]);
        if nargout >= 6
            hist = [hist; [outerits mvs sigma run rvn found ...
                resetcriteria orthv orthq restarts]];
        end
    end
    if display == 2
        fprintf('iter %5d mv %6d sig %7.6e Ru %7.5e conv %3d\n',...
            [outerits mvs sigma run found]);
    end
    if display <= 2 && nargout >= 6
        hist = [hist; [outerits mvs sigma run found]];
    end
    
    if run <= tol
        %%%Locking procedure%%%
        cflag = 0;
        iter = 1;
        while run < tol %Assuming rvn = 0;
            rcf = 1;
            cflag = 1;
            found = found + 1;
            
            if HLock
                LV = [LV v];
                LQ = [LQ u];
                LSig = [LSig;sigma];
                LR = [LR run];
                iter = iter + 1;
                t = t+1;
            else
                t = t+1;
            end
            
           
            %Calculate next left residual
            sigma = sr(index(t),index(t));
            u = tempQ*ur(:,index(t));
            v = tempV*vr(:,index(t));
            vrold = vr(:,index(t:(min(size(vr,2),numOld))+t-1));
            if xtrStorage
                ru = tempAtQ*ur(:,index(t)) - sigma*v;
            else
                ru = A(u,transp) - sigma*v; mvs = mvs+1;
            end
            run = norm(ru);
           
        end
        
        if HLock
            t = 1;
            if (iter ~= size(index,1)+1 && cflag)
                V = tempV*vr(:,index(iter:end));
                Q = tempQ*ur(:,index(iter:end));
                tempQ = Q;
                tempV = V;
                R(1:k-iter+1,1:k-iter+1) = sr(index(iter:end),index(iter:end));
                k = k - (iter - 1);
                restarts = restarts + 1;
            end
        else
            if t-1 >= numValues
                endQ = tempQ*ur(:,index(1:t-1));
                endV = tempV*vr(:,index(1:t-1));
                endS = sr(index(1:t-1),index(1:t-1));
                for endi = 1:t-1
                    endRv = norm(A(endV(:,endi),notransp) - endS(endi,endi)*endQ(:,endi));
                    if xtrStorage
                        endRu = norm(tempAtQ*ur(:,index(endi)) - endS(endi,endi)*endV(:,endi));
                    else
                        endRu = norm(A(endQ(:,endi),transp) - endS(endi,endi)*endV(:,endi));
                    end
                    endR(endi) = sqrt(endRv^2 + endRu^2);
                    if endR(endi) > tol
                        t = endi;
                        found = endi -1;
                        run = endRu;
                        if endRu < tol
                            reset = 1;
                            restarts = -1;
                            rcf = 1;
                        end
                        break;
                    end
                end
                if t-1 >= numValues
                    LQ = endQ;
                    LV = endV;
                    LSig = endS;
                    LR = endR;
                    break;
                end
            end
        end
    end
    
    %%%Inner solver%%%
    shift = min(sigma*run,sigma^2-sigma_prev^2);
    g = @(x) x - v*(v'*x);
    f = @(x,~) g(A((A(x,notransp)),transp)-(sigma^2-shift)*x);
    [prv,iters,touch,~] = qmrs(f,sigma*ru,target,tol,...
        maxII,Pata,sigma^2,sigma^2-shift,touch,[]);
    mvs = mvs+2*iters;
    
    %%%Basis Expansion%%%
    if HLock
        X = [LV tempV];
        [V(:,k+1),~] = cgs(X,prv);
    else
        [V(:,k+1),~] = cgs(tempV,prv);
    end
    
    Q(:,k+1) = A(V(:,k+1),notransp); mvs = mvs + 1;
    
    if HLock
        X = [LQ tempQ];
        [Q(:,k+1),tempR] = cgs(X,Q(:,k+1));
        R(1:k+1,k+1) = tempR(size(LQ,2)+1:end);
    else
        [Q(:,k+1),R(1:k+1,k+1)] = cgs(tempQ,Q(:,k+1));
    end
    
    if xtrStorage
        AtQ(:,k+1) = A(Q(:,k+1),transp); mvs = mvs + 1;
    end
    
    % Increase current basis size
    k = k + 1;
    
    %%%Restart/Reset procedure%%%
    if k >= maxBasis %&& ~cflag
        restarts = restarts + 1;
        rc = 4*normA*macheps*sqrt(restarts); %reset criteria
        %Check for reset
        if rc*rcf > run
            rvn = norm(A(v,notransp)-sigma*u); mvs = mvs + 1;
            rcf = rvn/(rc) * 1.22; %Modify rc to match the actual rvn
            if (rc*rcf > run)
                reset = 1;
                rcf = 1; %Reset convergence factor to default conservative value
                restarts = 0;
                if ~HLock
                    t = 1;
                    found = 0;
                end
            end
        end
        
        [ur,sr,vr]=svd(R(1:k,1:k));
        [~,index] = sort(abs(target - diag(sr)));
        y1 = vr(:,index(1:minRS));
        
        %{
        for itr = 1:t-1
            if sr(index(itr),index(itr)) - stored_sigma(itr) > tol
                t = itr;
                disp('Avoiding misconvergence')
            end
        end
        %}
        
        if numOld ~= 0
            yold = cgs(y1, [vrold(:,1:numOld);zeros(1,numOld)]);
            thickright = [y1 yold];
        else
            thickright = y1;
        end
        
        oldk = k;
        k = size(thickright,2);
        V(:,1:k) = V(:,1:oldk)*thickright;
        
        if (reset)
            V(:,1:k) = cgs(LV,V(:,1:k));
            Q(:,1:k) = A(V(:,1:k),notransp); mvs = mvs+k;
            [Q(:,1:k),R(1:k,1:k)] = qr(Q(:,1:k),0);
            if xtrStorage
                AtQ(:,1:k) = A(Q(:,1:k),transp); mvs = mvs+k;
            end
            reset = 0;
        else
            y2 = vr(:,index(minRS+1:end));
            w1 = ur(:,index(1:minRS));
            w2 = ur(:,index(minRS+1:end));
            s1 = sr(index(1:minRS),index(1:minRS));
            s2 = sr(index(minRS+1:end),index(minRS+1:end));
            
            if numOld ~= 0
                x = s2*(y2'*yold);
                [q,r] = qr(x,0);
                Qtilde = [w1 w2*q];
                rsize = size(s1,2)+size(r,2);
                R(1:rsize,1:rsize) = blkdiag(s1,r);
            else
                R(1:size(s1,2),1:size(s1,2)) = s1;
                Qtilde = w1;
            end
            
            Q(:,1:k) = Q(:,1:oldk)*Qtilde;
            if xtrStorage
                AtQ(:,1:k) = AtQ(:,1:oldk)*Qtilde;
            end
        end
        
    end %restart/reset procedure
    
end %while

if isempty(LR)
    LR = [LR min(hist(:,4))];
end

if Transpose
    Temp = LQ;
    LQ = LV;
    LV = Temp;
end


switch nargout
    case 1
        varargout{1} = LSig;
    case 2
        varargout{1} = LSig;
        varargout{2} = LR;
    otherwise
        varargout{1} = LQ;
        varargout{2} = LSig;
        varargout{3} = LV;
        if nargout >= 4
            varargout{4} = LR;
        end
        if nargout >= 5
            stats.numMatvecs = mvs;
            stats.elapsedTime = toc;
            stats.aNorm = normA;
            varargout{5} = stats;
        end
        if nargout >= 6
            varargout{6} = hist;
        end
        if nargout > 6
            warning('Too many output arguments');
        end
end

end %function

function result = singleA(A,x,tr)
if strcmp(tr,'transp')
    result = single(A'*double(x));
else
    result = single(A*double(x));
end
end


function result = doubleA(A,x,tr)
if strcmp(tr,'transp')
    result = A'*x;
else
    result = A*x;
end
end