function [varargout] = GKD(A,varargin)
%% GKD   Find a few singular values and vectors for large sparse matrices
% 
% S = GKD(A,k,target) returns the k singular values closest to target
% 
% S = GKD(A,k,target,opts) configures additional options within the 
% opts structure.
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
%   opts.LBD                1 = Start with LBD basis up to maxBasis-1
%   opts.BlockSize          Number of vectors in each block
%   opts.noCheck            Removes double-checking of soft-locked vectors
%
%
% S = GKD(A,k,target,opts,P) adds a preconditioner (for A'A)
%
% [U,S,V] = GKD(A,...) computes the singular vectors as well. 
% If A is M-by-N and K singular values are computed, then U is M-by-K
% with orthonormal columns, S is K-by-K diagonal, and V is N-by-K with
% orthonormal columns.
% 
% [U,S,V,H] = GKD(A,...) also returns convergence history

inputs = varargin;
[m,n,transp,notransp,Transpose,numValues,target,opts,Pata] = varginParse(A,inputs);


[ERROR,LBDflag,BlS,normA,tol,maxBasis,maxMVs,v0,display,minRS,numOld,maxII,noCheck] = optsParse(A,m,n,numValues,opts);
if ERROR ~= 0
    for i = nargout:-1:1
        varargout{i} = [];
    end
    return;
end

if minRS < numValues+BlS
    error('minRestartSize too small for soft-locking\n\tIncrease minRS > %d',numValues+BlS);
end

if ~isa(A,'function_handle')
    A = @(x,tr) doubleA(A,x,tr);
end
V = zeros(n,maxBasis);
Q = zeros(m,maxBasis);
R = zeros(maxBasis,maxBasis);
macheps = eps;


%If tol is too small, it may not be possible to converge
if tol < 5*normA*macheps
    error("Tolerance too small. opts.tol must be larger than %e",string(5*macheps));
end

reset = 0;
restarts = 0;
mvs = 0;     %MV Counter
k = BlS;       %Current basis size
outerits = 0;
touch = 1;
sigma_prev = 0;
sigma = 0;
rcf = 1; %reset criteria factor
hist = {};
found = 0;
converged = zeros(numValues+BlS,1);

%Storage for returned SVD and residuals
LV = [];
LQ = [];
LSig = [];
LR = [];

tic
if LBDflag
    mb = maxBasis-BlS;
    [Q(:,1:mb), R(1:mb,1:mb), V(:,1:mb)] = LBD(A,v0,floor(mb/BlS),Transpose);
    k = mb;
    mvs = mvs + 2*k;
else
    [v0,~]=qr(v0,0);
    V(:,1:k) = v0;
    [Q(:,1:k),R(1:k,1:k)] = qr(A(v0,notransp),0); mvs = mvs + 1;
end

%% Main iteration
while (maxMVs <= 0 || mvs<maxMVs) && found < numValues
    t = find(converged == 0, BlS, 'first');
    t = t(t <= k);
    t_size = size(t,1);
    tempQ = Q(:,1:k);
    tempV = V(:,1:k);
    cflag = 0; %Convergence flag
    
    % Extraction process
    [ur,sr,vr]=svd(R(1:k,1:k));
    if target==inf
        index = 1:k;
    elseif target == 0
        index = k:-1:1;
    else
        [~,index] = sort(abs(target - diag(sr)));
    end
    sigma_prev = sigma;
    sigma = sr(index(t),index(t));
    u = tempQ*ur(:,index(t));
    v = tempV*vr(:,index(t));
    vrold = vr(:,index(t:(min(size(vr,2),numOld))+t-1));
    
    %Calculate left residual
    ru = A(u,transp) - v*sigma; mvs = mvs + 1;
    run = vecnorm(ru);
    
    outerits = outerits+1;
    
    hist{end+1} = {toc, outerits, mvs , sigma, run, sum(converged)};
    fprintf('Time: %7.3f Iter: %4d Matvecs: %4d Converged: %4d Restarts: %4d\n',...
        toc, outerits, mvs, sum(converged), restarts);
    
    indices = find(run <= tol);

    if any(indices)
        %%%Locking procedure%%%
        cflag = 1;
        converged(t(indices)) = 1;
        if sum(converged) >= numValues
            LQ = tempQ*ur(:,index(1:numValues));
            LV = tempV*vr(:,index(1:numValues));
            LSig = sr(index(1:numValues),index(1:numValues));
            if ~noCheck
                for endi = 1:numValues
                    endRv = norm(A(LV(:,endi),notransp) - LSig(endi,endi)*LQ(:,endi));
                    endRu = norm(A(LQ(:,endi),transp) - LSig(endi,endi)*LV(:,endi));
                    endR(endi) = sqrt(endRv^2 + endRu^2);
                end
                indices = find(endR > tol);
                if any(indices)
                    converged(indices) = 0;
                    reset = 1;
                    restarts = 0;
                    rcf = 1;
                else
                    LR = endR;
                    break;
                end
            else
                break;
            end
        end
    end

    
    %%%Inner solver%%%
    prv = [];
    for innerIter = 1:size(ru,2)
            sigmaInner = sigma(innerIter,innerIter);
            shift = sigmaInner*run(innerIter); %use sigmaInner^2-sigma_prev^2 if smaller...
            g = @(x) x - v*(v'*x);
            f = @(x,~) g(A((A(x,notransp)),transp)-(sigmaInner^2-shift)*x);
            [prv(:,innerIter),iters,touch,~] = qmrs(f,sigmaInner*ru(:,innerIter),target,tol,...
                maxII,Pata,sigmaInner^2,sigmaInner^2-shift,touch,[]);
            mvs = mvs+2*iters;
    end
    
    [V(:,k+1:k+size(t,1)),~] = cgs(tempV,prv);
    
    Q(:,k+1:k+t_size) = A(V(:,k+1:k+t_size),notransp); mvs = mvs + 1;

    [Q(:,k+1:k+t_size),R(1:k+t_size,k+1:k+t_size)] = cgs(tempQ,Q(:,k+1:k+t_size));
    
    % Increase current basis size
    k = k + t_size;
    
    %%%Restart/Reset procedure%%%
    if k >= maxBasis %&& ~cflag
        
        restarts = restarts + 1;
        rc = 4*normA*macheps*sqrt(restarts); %reset criteria
        %Check for reset
        if any(rc*rcf > run)
            rvn = norm(A(v,notransp)-u*sigma); mvs = mvs + size(u,2);
            rcf = rvn/(rc) * 1.22; %Modify rc to match the actual rvn
            if (rc*rcf > run)
                reset = 1;
                rcf = 1; %Reset convergence factor to default conservative value
                restarts = 0;
            end
        end
        
        [ur,sr,vr]=svd(R(1:k,1:k));
        if target == inf
            index = 1:k;
        elseif target == 0
            index = k:-1:1;
        else
            [~,index] = sort(abs(target - diag(sr)));
        end
        
        index1 = index(1:minRS);
        index2 = index(minRS+1:end);
        y1 = vr(:,index1);
        
        if numOld ~= 0
            yold = cgs(y1, [vrold(:,1:numOld);zeros(BlS,numOld)]);
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
            reset = 0;
        else
            y2 = vr(:,index2);
            w1 = ur(:,index1);
            w2 = ur(:,index2);
            s1 = sr(index1,index1);
            s2 = sr(index2,index2);
            
            if numOld ~= 0
                x = s2*(y2'*yold);
                [q,r] = qr(x,0);
                Qtilde = [w1 w2*q];
                rsize = size(s1,2)+size(r,2);
                R(1:rsize,1:rsize) = blkdiag(s1,r);
            else
                R(1:size(s1,1),1:size(s1,2)) = s1;
                Qtilde = w1;
            end
            
            Q(:,1:k) = Q(:,1:oldk)*Qtilde;
        end
        
    end %restart/reset procedure
    
end %while

if isempty(LR)
    LR = hist{end}{5};
end

hist{end+1} = {toc, outerits, mvs , sigma, LR, size(LQ,2)};

if Transpose
    Temp = LQ;
    LQ = LV;
    LV = Temp;
end


switch nargout
    case 1
        varargout{1} = LS;
    case 3
        varargout{1} = LQ;
        varargout{2} = LS;
        varargout{3} = LV;
    otherwise
        varargout{1} = LQ;
        varargout{2} = LS;
        varargout{3} = LV;
        varargout{4} = hist;
end

end %function

function result = doubleA(A,x,tr)
if strcmp(tr,'transp')
    result = A'*x;
else
    result = A*x;
end
end