function [LU,LS,LV,history] = GKD(A,numVals,varargin)
%% GKD   Find a few singular values and vectors for large sparse matrices
%
% [U,S,V,H] = GKD(A,numVals,target) returns numVals singular triplets
%       closest to target with convergence history
%
% [U,S,V,H] = GKD(A,numVals,target,...) configures additional options
%
% Options for GKD (Name, Value pairs)
% tol           residual norm tolerance
% maxMV         maximum number of matrix-vector multiplications
% normA         norm(A) estimate
% showHist      Prints partial history to console if set
% v0            Initial vector for V
% b             Block size
% maxBasis      max number of basis vectors in V,U
% minRS         number of vectors to maintain after restart
% numOld        number of +k vectors to keep during restart
% maxII         max number of inner solver iterations
% seed          random seed
% m             number of rows in A (if A is a function_handle)
% n             number of cols in A (if A is a function_handle)
% P             preconditioner for AtA

% If A is M-by-N and K singular values are computed, then U is M-by-K
% with orthonormal columns, S is K-by-K diagonal, and V is N-by-K with
% orthonormal columns.


%% Input Parsing %%

p = inputParser();
addRequired(p,'A');
addRequired(p,'numVals');
addOptional(p,'mode','resid');
addParameter(p,'tol',1e-6);
addParameter(p,'sort_order','L');
addParameter(p,'target_func',[]);
addParameter(p,'stop_func',[]);
addParameter(p,'maxMV',inf);
addParameter(p,'normA',-1);
addParameter(p,'showHist',0);
addParameter(p,'v0',[]);
addParameter(p,'b',1); %block size
addParameter(p,'maxBasis',15);
addParameter(p,'minRS',7);
addParameter(p,'numOld',1);
addParameter(p,'maxII',0);
addParameter(p,'seed','shuffle')
addParameter(p,'m',-1);
addParameter(p,'n',-1);
addParameter(p,'P',1);

parse(p,A,numVals,varargin{:});

A = p.Results.A;
numVals = p.Results.numVals;
mode = p.Results.mode;
tol = p.Results.tol;
sort_order = p.Results.sort_order;
target_func = p.Results.target_func;
stop_func = p.Results.stop_func;
maxMV = p.Results.maxMV;
normA = p.Results.normA;
showHist = p.Results.showHist;
v0 = p.Results.v0;
b = p.Results.b;
maxBasis = p.Results.maxBasis;
minRS = p.Results.minRS;
numOld = p.Results.numOld;
maxII = p.Results.maxII;
seed = p.Results.seed;
m = p.Results.m;
n = p.Results.n;
P = p.Results.P;

rng(seed);
if m == -1 || n == -1, [m,n] = size(A); end
transp = 'transp';
notransp = 'notransp';

%Reduce size of A'*A for rectangular matrices
if m < n
    temp = m; m = n; n = temp;
    transp = 'notransp'; notransp = 'transp';
end

if ~isa(A,'function_handle')
    A = @(x,tr) doubleA(A,x,tr);
end

if strcmp(mode,'resid') || isempty(target_func) || isempty(stop_func)
    target_func = @(normA,k,sr,f_vecs,f_resid) resid_target(tol,normA,k,sr,f_vecs,f_resid);
    stop_func = @(normA,k,sr,f_vecs,f_resid) resid_stop(tol,numVals,normA,k,sr,f_vecs,f_resid);
end

%Set min restart to at least one block larger than desired # of values
if minRS < numVals
    minRS = numVals+b;
end

if maxBasis < max(minRS+3*b,floor(1.3*minRS))
    maxBasis = max(minRS+3*b,floor(1.3*minRS));
end

%Set initial vectors if not given
if isempty(v0), v0 = randn(n,b); end

%% Set up variables
reset = 0;
restarts = 0;
mvs = 0;     %MV Counter
outerits = 0;
touch = 1;
rcf = 1; %reset criteria factor
mipv = 1; %minimum iterations per value
history = struct(); histindex = 1;
conv = zeros(maxBasis + 2*b,1);

%Storage for returned SVD and residuals
LV = [];
LU = [];
LS = [];

%Preallocate memory for bases, and R
V = zeros(n,maxBasis);
Q = zeros(m,maxBasis);
R = zeros(maxBasis,maxBasis);

%Allocate memory for a block of vectors/residuals
u = zeros(m,b);
v = zeros(n,b);
ru = zeros(n,b);

%Allocate memory for a running storage of recent residual norms
allrun = inf(maxBasis,1);

%% Set up method
starttime = tic;
[v0,~]=qr(v0,0);
k = size(v0,2);
V(:,1:k) = v0;
[Q(:,1:k),R(1:k,1:k)] = qr(A(v0,notransp),0); mvs = mvs + k;

%% Main Iteration
while (maxMV <= 0 || mvs<maxMV)
    One2K = 1:k;
    tempQ = Q(:,One2K);
    tempV = V(:,One2K);
    
    [ur,sr,vr]=svd(R(One2K,One2K));  %sr goes from largest to smallest
    sr = diag(sr);
    if max(sr) > normA
        normA = max(sr);
    end
    
    %Reorder SVD results
    order = f_order(k,sort_order,sr);
    ur = ur(:,order);
    sr = sr(order);
    vr = vr(:,order);
    
    f_resid = @(indices,compute) getResidNorms(indices,compute);
    f_vecs = @(indices) getVecs(indices);
    f_A = @(x,tr) A(x,tr);
    
    index = target_func(normA,k,sr,f_vecs,f_resid);
    [done,numVals] = stop_func(normA,k,sr,f_vecs,f_resid);

    index = index(1:min(length(index),b));
    
    if showHist
        fprintf('Time: %7.3f Iter: %4d Matvecs: %4d Restarts: %4d\n',...
            toc(starttime), outerits, mvs, restarts);
    end
    if done
        LU = tempQ*ur(:,1:numVals);
        LV = tempV*vr(:,1:numVals);
        LS = sr(1:numVals);
        break;
    end
    
    t_size = length(index);
    u = tempQ*ur(:,index);
    v = tempV*vr(:,index);
    s = sr(index);
    ru = A(u,'transp') - v*diag(s);
    run = vecnorm(ru);
    allrun(index) = run;
    
    prv = zeros(size(ru));
    for j = 1:size(ru,2)
        si = s(j);
        shift = si^2 - si*run(j);
        g = @(x) x - v*(v'*x);
        f = @(x,~) g(A((A(x,notransp)),transp)-shift*x);
        [prv(:,j),iters,touch] = qmrs(f,si*ru(:,j),normA*tol,maxII,P,si^2,shift,touch);
        mvs = mvs+2*iters;
    end
    
    %% Basis Expansion
    [V(:,k+1:k+t_size),~] = cgs(tempV,prv);
    Q(:,k+1:k+t_size) = A(V(:,k+1:k+t_size),notransp); mvs = mvs + t_size;
    [Q(:,k+1:k+t_size),R(1:k+t_size,k+1:k+t_size)] = cgs(tempQ,Q(:,k+1:k+t_size));
    k = k + t_size;
    
    outerits = outerits+1;
    
    %% Restart/Reset procedure %%
    if k >= maxBasis
        vrold = vr(:,index(1:numOld));  %Used for restarting with GD+k
        restarts = restarts + 1;
        rc = 4*normA*eps*sqrt(restarts); %reset criteria
        if ~reset, [reset,newmvs] = checkReset(rc,rcf,run,A,u,s,v,notransp); end
        mvs = mvs + newmvs;
        [ur,sr,vr] = svd(R(1:k,1:k));
        sr = diag(sr);
        
        order = f_order(k,sort_order,sr);
        
        [Vtilde,yold] = updateV(vr,order,minRS,vrold);
        oldk = k; k = size(Vtilde,2);
        V(:,1:k) = V(:,1:oldk)*Vtilde;
        if reset
            reset = 0; restarts = 0; rcf = 1;
            allrun = inf(maxBasis,1);
            V(:,1:k) = qr(V(:,1:k),0);
            Q(:,1:k) = A(V(:,1:k),notransp);
            [Q(:,1:k),R(1:k,1:k)] = qr(Q(:,1:k),0);
            mvs = mvs + k;
        else
            [Qtilde,R] = restartQ(ur,sr,vr,order,minRS,yold);
            Q(:,1:k) = Q(:,1:oldk)*Qtilde;
        end
    end
    
end

if strcmp(transp,'notransp')
    Temp = LU;
    LU = LV;
    LV = Temp;
end
    
    %% Nested functions for getting Residuals and Vectors
    function [residnorm] = getResidNorms(i,compute)
        if compute
            residnorm = vecnorm(A(tempQ*ur(:,i),'transp') - tempV*(vr(:,i)*diag(sr(i))))';
            allrun(i) = residnorm(i);
        else
            residnorm = allrun(i);
        end
    end

    function [u,v] = getVecs(i)
        u = tempQ*ur(:,i);
        v = tempV*vr(:,i);
    end

end

function [index] = resid_target(tol,normA,k,~,~,f_resid)
    r = f_resid(1:k,0);
    index = find(r > normA*tol); %1e-6 is the tolerance
end

function [done,numVals] = resid_stop(tol,numVals,normA,k,~,~,f_resid)
    if k > numVals
        r = f_resid(1:numVals,0); %6 is the chosen number of values
        if all(r < normA*tol)
            done = 1;
        else
            done = 0;
        end
    else
        done = 0;
    end
end

function order = f_order(k,target,sr)
if strcmpi(target,'L')
    order = 1:k;
elseif strcmpi(target,'S')
    order = k:-1:1;
elseif isa(target,'double')
    [~,order] = sort(abs(target-sr));
else
    error('sort_order must be ''L'', ''S'' or a double')
end
end


function result = doubleA(A,x,tr)
if strcmp(tr,'transp')
    result = A'*x;
else
    result = A*x;
end
end
