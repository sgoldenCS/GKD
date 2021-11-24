function [U,S,V,HIST,UD] = GKD(A,numVals,varargin)
%% GKD   Find a few singular values and vectors for large sparse matrices
%
% [U,S,V,H,UD] = GKD(A,numVals) returns numVals singular triplets
%       closest to target with convergence history
%
% If A is M-by-N and K singular values are computed, then U is M-by-K
% with orthonormal columns, S is K-by-K diagonal, and V is N-by-K with
% orthonormal columns.
%
% [U,S,V,H,UD] = GKD(A,numVals,...) configures additional options
%
% Options for GKD (Name, Value pairs)
% tol           residual norm tolerance
% SIGMA         Sets which values are desired ('L' or 'S')
% target_fn  Function used for expanding the basis
%    [index,userdata] = target_fn(solverdata,userdata)
% stop_fn     Function used for stopping the solver
%    [done,numVals,userdata] = stop_fn(numVals,solverdata,userdata)
% maxMV         maximum number of matrix-vector multiplications
% normA         norm(A) estimate
% display      Prints partial history to console if set
% v0            Initial vector for V
% b             Block size
% minRestart         number of vectors to maintain after restart
% maxBasis      max number of basis vectors in V,U
% numOld        number of +k vectors to keep during restart
% maxII         max number of inner solver iterations
% seed          random seed
% m             number of rows in A (if A is a function_handle)
% n             number of cols in A (if A is a function_handle)
% P             preconditioner for AtA
%
% Default Options Settings
%
% tol          1e-6
% SIGMA        'L' (Largest)
% target_fn    Targeting based on first values with residuals above 'tol'
% stop_fn      Stopping based on residual tolerance ('tol')
% maxMV        -1 (No stopping based on matvecs)
% normA        Largest value seen (Accurate when SIGMA = 'L')
% display      0 (Off)
% v0           Gaussian random vectors
% b            1
% minRestart        numVals+max(b,15)
% maxBasis     max(minRestart+2*b,floor(1.3*minRestart))
% numOld       1
% maxII        0
% seed         'shuffle' (sets rng based on current time)
% P            1 (Identity matrix)


%% Input Parsing %%
p = inputParser();
addRequired(p,'A');
addRequired(p,'numVals');
addParameter(p,'tol',1e-6);
addParameter(p,'SIGMA','L');
addParameter(p,'target_fn','resid');
addParameter(p,'stop_fn','resid');
addParameter(p,'maxMV',inf);
addParameter(p,'maxTime',inf);
addParameter(p,'normA',-1);
addParameter(p,'display',0);
addParameter(p,'v0',[]);
addParameter(p,'b',1); %block size
addParameter(p,'maxBasis',-1);
addParameter(p,'minRestart',-1);
addParameter(p,'numOld',1);
addParameter(p,'maxQMR',0);
addParameter(p,'seed','shuffle');
addParameter(p,'m',-1);
addParameter(p,'n',-1);
addParameter(p,'P',1);
addParameter(p,'user_data',[]);

parse(p,A,numVals,varargin{:});

p = p.Results;
m = p.m; n = p.n;
normA = p.normA;

rng(p.seed);
if m == -1 || n == -1, [m,n] = size(p.A); end
transp = 'transp';
notransp = 'notransp';

%Reduce size of A'*A for rectangular matrices
if m < n
    temp = m; m = n; n = temp;
    transp = 'notransp'; notransp = 'transp';
end

if ~isa(p.A,'function_handle')
    p.A = @(x,tr) Afun(p.A,x,tr);
end

%Default Parameters for minRestart and maxBasis
if p.minRestart < p.numVals
    p.minRestart = max(7,p.numVals+p.b);
end

if p.maxBasis < p.minRestart
    p.maxBasis = max(15,p.numVals+3*p.b);
end

if isa(p.target_fn,'char')
    if strcmpi(p.target_fn,'prog_tol')
        p.target_fn = @(sd,ud) ...
            prog_tol_target(sd,ud);
    elseif strcmpi(p.target_fn,'resid')
        p.target_fn = @(sd,ud) ...
            resid_target(sd,ud);
    else
        error(string(p.target_fn)+' is not a known target function');
    end
end

if isa(p.stop_fn,'char')
    if strcmpi(p.stop_fn,'resid')
        p.stop_fn = @(numSV,sd,ud) resid_stop(numSV,sd,ud);
    else
        error(string(p.stop_fn)+' is not a known stopping function');
    end
end

%Preallocate memory for bases, and R
V = zeros(n,p.maxBasis);
U = zeros(m,p.maxBasis);
R = zeros(p.maxBasis,p.maxBasis);

% Set up variables
reset = 0;          %Reset flag
restarts = 0;       %Restart Counter
mvp = 0;            %MV Counter
outerits = 0;       %Iteration Counter
touch = 1;          %Used for QMR Convergence Criterion
rcf = 1;            %Reset criteria factor
HIST = [];          %Convergence History 

allrun = inf(p.maxBasis,1);   %Storage for Residual Norms

starttime = tic;
if isempty(p.v0)
    k = p.b;
    [V(:,1:k),~] = qr(randn(n,p.b),0);
else
    k = size(p.v0,2);
    [V(:,1:k),~] = qr(p.v0,0);
end
[U(:,1:k),R(1:k,1:k)] = qr(p.A(V(:,1:k),notransp),0); mvp = mvp + k;

while mvp < p.maxMV && toc(starttime) < p.maxTime
    [ur,sr,vr] = svd(R(1:k,1:k));
    sr = diag(sr);
    
    if max(sr) > normA
        normA = max(sr);
    end
    
    %%%% It would be nice not to do this if possible %%%%
    %Reorder SVD results
    order = f_order(k,p.SIGMA,sr);
    ur = ur(:,order);
    sr = sr(order);
    vr = vr(:,order);
    
    solver_data = struct('mvp',mvp,'iters',outerits,'time',toc(starttime), ...
        'numVals',p.numVals,'normA',normA,'k',k,'s',sr, ...
        'b',p.b,'tol',p.tol,'rn',allrun);
    
    [index,p.user_data] = p.target_fn(solver_data,p.user_data);
    
    u = U(:,1:k)*ur(:,index);
    v = V(:,1:k)*vr(:,index);
    s = sr(index);
    ru = p.A(u,transp) - v*diag(s);
    mvp = mvp + length(index);
    allrun(index) = vecnorm(ru);
    
    if p.display
        fprintf('Time: %7.3f Iter: %4d Matvecs: %4d Restarts: %4d Num Conv: %4d Min Unconverged Resid: %7.3e\n',...
            toc(starttime), outerits, mvp, restarts, length(find(allrun < normA*p.tol)), min(allrun(allrun > normA*p.tol)));
    end
    HIST = [HIST; toc(starttime), outerits, mvp, restarts, length(find(allrun < normA*p.tol)), min(allrun(allrun > normA*p.tol))];
    solver_data.rn = allrun;
    [done,p.numVals,p.user_data] = p.stop_fn(p.numVals,solver_data,p.user_data);
    
    if isempty(index) || done
        U = U(:,1:k)*ur(:,1:min(p.numVals,k));
        V = V(:,1:k)*vr(:,1:min(p.numVals,k));
        S = sr(1:min(p.numVals,k));
        break;
    end
    
    
    cb_size = size(ru,2);
    if p.maxQMR > 0
        for j = 1:cb_size
            si = s(j);
            shift = si^2 - si*run(j);
            g = @(x) x - v*(v'*x);
            f = @(x,~) g(p.A((p.A(x,notransp)),transp)-shift*x);
            [ru(:,j),iters,touch] = qmrs(f,si*ru(:,j),normA*p.tol,p.maxII,p.P,si^2,shift,touch);
            mvp = mvp+2*iters;
        end
    end
    
    %% Basis Expansion
    [V(:,k+1:k+cb_size),~] = cgs(V(:,1:k),ru); %Ortho on V
    U(:,k+1:k+cb_size) = p.A(V(:,k+1:k+cb_size),notransp); mvp = mvp + cb_size;
    [U(:,k+1:k+cb_size),R(1:k+cb_size,k+1:k+cb_size)] = cgs(U(:,1:k),U(:,k+1:k+cb_size)); %Ortho on U
    
    k = k + cb_size;
    outerits = outerits+1;
    
    %% Restart/Reset procedure %%
    if k >= p.maxBasis
        vrold = vr(:,index(1:p.numOld));  %Used for restarting with GD+k
        restarts = restarts + 1;
        rc = 4*normA*eps*sqrt(restarts); %reset criteria
        if ~reset, [reset,newmvp] = checkReset(rc,rcf,allrun(index),p.A,u,s,v,notransp); end
        mvp = mvp + newmvp;
        [ur,sr,vr] = svd(R(1:k,1:k));
        sr = diag(sr);
        
        order = f_order(k,p.SIGMA,sr);
        
        [Vtilde,yold] = updateV(vr,order,p.minRestart,vrold);
        oldk = k; k = size(Vtilde,2);
        V(:,1:k) = V(:,1:oldk)*Vtilde;
        if reset
            reset = 0; restarts = 0; rcf = 1;
            allrun = inf(p.maxBasis,1);
            [V(:,1:k),~] = qr(V(:,1:k),0);
            U(:,1:k) = p.A(V(:,1:k),notransp);
            [U(:,1:k),R(1:k,1:k)] = qr(U(:,1:k),0);
            mvp = mvp + k;
        else
            [Utilde,R] = restartU(ur,sr,vr,order,p.minRestart,yold);
            U(:,1:k) = U(:,1:oldk)*Utilde;
        end
        allrun(p.minRestart+1:end) = inf;
    end
end

if strcmp(transp,'notransp')
    Temp = U;
    U = V;
    V = Temp;
end

UD = p.user_data;

end

%% Standard Residual Stopping Criteria
function [done,numVals,ud] = resid_stop(numVals,sd,ud)

done = 0;
if sd.k > numVals
    r = sd.rn(1:numVals);
    if all(r < sd.normA*sd.tol)
        done = 1;
    end
end
end

%% Target based on residual tolerance
% Returns at most one block size of indices
% Index will be empty if all ||r|| are below ||A||*tol and k >= numVals
function [index,ud] = resid_target(sd,ud)
r = sd.rn(1:sd.k);
index = find(r > sd.normA*sd.tol);
if sd.k < sd.numVals && isempty(index)
    [~,index] = sort(r,'descend');
end
index = index(1:min(sd.b,length(index)));
end

%% Target based on a progressively more restrictive residual tolerance
function [index,ud] = prog_tol_target(sd,ud)
r = sd.rn(1:min(sd.numVals,sd.k));
r(r > sd.normA) = sd.normA;
% If tolerance has been met for all sd.numVals, return empty index
% to stop the solver
if length(r) > sd.numVals && all(r(1:sd.numVals) < sd.tol*sd.normA)
    index = [];
    return
end
r_order = abs(ceil(log10(r./sd.normA)./2));
[~,index] = sort(r_order);
%Max length of index should be p.b
if length(index) < sd.b
    index = [index; [sd.numVals+1:sd.k]'];
    index = index(1:min(sd.b,length(index)));
else
    index = index(1:sd.b);
end

end


%% Function if given exact matrix A
function result = Afun(A,x,tr)
if strcmp(tr,'transp')
    result = A'*x;
else
    result = A*x;
end
end

%% Function to order singular triplets based on target
function order = f_order(k,target,sr)
if strcmpi(target,'L')
    order = 1:k;
elseif strcmpi(target,'S')
    order = k:-1:1;
elseif isa(target,'numeric')
    [~,order] = sort(abs(target-sr));
else
    error('SIGMA must be ''L'', ''S'' or numeric')
end
end
