function [U,S,V,HIST,UD] = GKD(A,numVals,varargin)
%% GKD   Find a few singular triplets for large sparse matrices
%
% [U,S,V,H,UD] = GKD(A,numVals) returns numVals singular triplets
%       closest to target with convergence history
%
% If A is M-by-N and K singular values are computed, then U is M-by-K
% with orthonormal columns, S is K-by-K diagonal, and V is N-by-K with
% orthonormal columns.
%
% HIST returns the total time, number of iterations, matrix vector products,
% number of restarts, number of converged values (residual < normA*tol), and
% the smallest unconverged residual.
%
% UD returns any userdata needed or generated by stopping functions.
% The standard stopping criteria for GKD do not utilize this structure.
%
% If A is a function, A(x,'notransp') must return the matrix vector
% product A*x and A(x,'transp') must return A'*x. Additionally, the
% optional name value pairs 'm' and 'n' must be configured with the
% number of rows and columns of A respectively.
%
% [U,S,V,H,UD] = GKD(A,numVals,...) configures additional options
%
% Options for GKD (Name, Value pairs)
% tol           residual norm tolerance
% SIGMA         Sets which values are desired ('L' or 'S')
% target_fn     Function used for expanding the basis ('resid' or 'prog_tol')
% stop_fn       Function used for stopping the solver
%                   [done,numVals,userdata] = stop_fn(numVals,solverdata,userdata)
% maxMV         maximum number of matrix-vector multiplications
% maxTime       maximum solver time
% normA         norm(A) estimate
% display       Prints partial history to console if set
% v0            initial vectors with min(size(A)) rows
% b             block size
% minRestart    number of vectors to maintain after restart
% maxBasis      max number of basis vectors in V,U
% numOld        number of +k vectors to keep during restart (must be <= b)
% maxQMR        max number of QMR iterations
% seed          random seed
% m             number of rows in A (if A is a function_handle)
% n             number of cols in A (if A is a function_handle)
% P             preconditioner for AtA
% userdata      external data needed for custom stopping functions
%
% Default Options Settings
%
% tol          1e-6
% SIGMA        'L' (Largest)
% target_fn    'resid'
% stop_fn      'resid' (Stops when the first numVals residuals are < tol)
% maxMV        inf (No stopping based on matvecs)
% maxTime      inf (No time based stopping)
% normA        Largest value seen (Accurate when SIGMA = 'L')
% display      0 (Off)
% v0           Gaussian random vectors
% b            1
% minRestart   max(7,p.numVals+5)
% maxBasis     max([15,p.minRestart+4*p.b,floor(1.3*p.minRestart)])
% numOld       -1 (Uses current block size for +k)
% maxQMR       0
% seed         'shuffle' (sets rng based on current time)
% P            []
% userdata     []
%
% Using the stop_fn option
%
% Requires a funtion_handle with the following form:
%       [done,numVals,userdata] = stop_fn(numVals,solverdata,userdata)
% Inputs:
% numVals      The number of singular values desired. This is given as an
%              input so the user function does not need to assign the
%              numVals output variable unless desired.
% solverdata   Includes the following information:
%                   mvp: number of matrix vector products,
%                   outerits: number of outer iterations,
%                   time: total time elapsed,
%                   normA: ||A||,
%                   k: current basis size,
%                   s: a (k x 1) vector of ritz values,
%                   b: current block size,
%                   tol: minimum residual tolerance,
%                   rn: a (maxBasis x 1) vector of the most recent residual norms
% userdata     Inputs the userdata structure given as an optional parameter.
%              This can be updated inside of the stopping function as long
%              as the updates are returned in the third return argument
%
% Outputs:
% done         Convergence flag. If set, the solver will stop and return
%              the most recent numVals Ritz vectors and values.
% numVals      Determines the number of values to return when the 'done'
%              flag is set. Can be set separately, but should always be
%              less than or equal to the numVals input. Doing this without
%              setting the 'done' flag may change targeting to improve
%              convergence.
% userdata     Allows for updates to the userdata structure to be stored
%              for later iterations.



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
addParameter(p,'numOld',-1);
addParameter(p,'maxQMR',0);
addParameter(p,'seed','shuffle');
addParameter(p,'m',-1);
addParameter(p,'n',-1);
addParameter(p,'P',[]);
addParameter(p,'userdata',[]);

parse(p,A,numVals,varargin{:});

p = p.Results;
m = p.m; n = p.n;
normA = p.normA;

rng(p.seed);
if m == -1 || n == -1
    if isa(p.A,'function_handle')
        error('Using a Function Handle for A requires dimensions [m,n] as parameters');
    else
        [m,n] = size(p.A);
    end
end
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
    p.minRestart = max(7,p.numVals+5);
end

if p.maxBasis < p.minRestart
    p.maxBasis = max([15,p.minRestart+4*p.b,floor(1.3*p.minRestart)]);
end

%Restrict +k restarting to +b
if p.numOld > p.b
    error('+'+string(p.numOld)+' restarting is not available with block size '+string(p.b));
end

if isa(p.target_fn,'char')
    if strcmpi(p.target_fn,'prog_tol')
        p.target_fn = @(sd,ud,A,U,V) ...
            prog_tol_target(sd,ud,A,U,V);
    elseif strcmpi(p.target_fn,'resid')
        p.target_fn = @(sd,ud,A,U,V) ...
            resid_target(sd,ud,A,U,V);
    elseif strcmpi(p.target_fn,'large')
        p.target_fn = @(sd,ud,A,U,V) ...
            large_target(sd,ud,A,U,V);
    else
        error(string(p.target_fn)+' is not a known target function');
    end
end

if isa(p.stop_fn,'char')
    if strcmpi(p.stop_fn,'resid')
        p.stop_fn = @(numSV,sd,ud,A,U,V) resid_stop(numSV,sd,ud,A,U,V);
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
outerits = 1;       %Iteration Counter
touch = 1;          %Used for QMR Convergence Criterion
rcf = 1;            %Reset criteria factor
HIST = [];          %Convergence History

allrun = inf(p.maxBasis,1);   %Storage for Residual Norms
resid_est = inf(p.maxBasis,1);

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
    
    %Reorder SVD results
    order = f_order(k,p.SIGMA,sr);
    ur = ur(:,order);
    sr = sr(order);
    vr = vr(:,order);
    
    solver_data = struct('mvp',mvp,'iters',outerits,'time',toc(starttime), ...
        'numVals',p.numVals,'normA',normA,'k',k,'s',sr,'transp',transp,'notransp',notransp,...
        'b',p.b,'tol',p.tol,'rn',allrun,'ur',ur,'vr',vr,'resid_est',resid_est);
    
    [index,p.userdata] = p.target_fn(solver_data,p.userdata,p.A,U,V);
    
    numVecs = min(length(index),p.b);
    if numVecs ~= 0
        u = U(:,1:k)*ur(:,index(1:numVecs));
        v = V(:,1:k)*vr(:,index(1:numVecs));
        s = sr(index(1:numVecs));
        ru = p.A(u,transp) - v*diag(s);
        mvp = mvp + numVecs;
        run = vecnorm(ru);
        allrun(index(1:numVecs)) = run;
    end
    
    if p.display
        fprintf('Time: %7.3f Iter: %4d Matvecs: %4d Restarts: %4d Num Conv: %4d Min Resid: %7.3e\n',...
            toc(starttime), outerits, mvp, restarts, length(find(allrun < normA*p.tol)), ...
            min(allrun(allrun > normA*p.tol))/normA);
    end
    HIST = [HIST; toc(starttime), outerits, mvp, restarts, ...
        length(find(allrun < normA*p.tol)), min(allrun(allrun > normA*p.tol))/normA];
    
    solver_data.rn = allrun;
    [done,p.numVals,p.userdata] = p.stop_fn(p.numVals,solver_data,p.userdata,p.A,U,V);
    
    if numVecs == 0 || all(allrun(1:p.numVals) < normA*p.tol) || done
        U = U(:,1:k)*ur(:,1:min(p.numVals,k));
        V = V(:,1:k)*vr(:,1:min(p.numVals,k));
        S = diag(sr(1:min(p.numVals,k)));
        break;
    end
    
    cb_size = size(ru,2);
    if p.maxQMR > 0 || ~isempty(p.P)
        for j = 1:cb_size
            si = s(j);
            shift = si^2 - si*run(j);
            g = @(x) x - v*(v'*x);
            f = @(x,~) g(p.A(p.A(x,notransp),transp)-shift*x);
            [ru(:,j),iters,touch] = qmrs(f,si*ru(:,j),p.SIGMA,normA*p.tol,p.maxQMR,p.P,si^2,shift,touch);
            mvp = mvp+2*iters;
        end
    end
    
    %% Basis Expansion
    [V(:,k+1:k+cb_size),~] = cgs(V(:,1:k),ru); %Ortho on V
    U(:,k+1:k+cb_size) = p.A(V(:,k+1:k+cb_size),notransp); mvp = mvp + cb_size;
    [U(:,k+1:k+cb_size),R(1:k+cb_size,k+1:k+cb_size)] = cgs(U(:,1:k),U(:,k+1:k+cb_size)); %Ortho on U
    
    resid_est = vecnorm(R(1:k,k+1:k+cb_size)'*ur,2,1);
    resid_est(resid_est < normA*eps) = normA*eps;

    k = k + cb_size;
    outerits = outerits+1;
    
    %% Restart/Reset procedure %%
    if k > p.maxBasis - p.b
        %Get +k vectors
        if p.numOld == -1
            vrold = vr(:,index(1:numVecs));
        elseif p.numOld > numVecs
            warning('Not enough +k restart vectors. Restarting with +'+string(numVecs));
            vrold = vr(:,index(1:numVecs));
        else
            vrold = vr(:,index(1:p.numOld));
        end
        
        restarts = restarts + 1;
        rc = 4*normA*eps*sqrt(restarts); %reset criteria
        if ~reset, [reset,newmvp] = checkReset(rc,rcf,allrun(index),p.A,u,s,v,notransp); end
        mvp = mvp + newmvp;
        
        [ur,sr,vr] = svd(R(1:k,1:k));
        sr = diag(sr);
        
        order = f_order(k,p.SIGMA,sr);
        
        [Vtilde,yold] = updateV(vr,order,p.minRestart,vrold);
        oldk = k; k = size(Vtilde,2);
        resid_est = resid_est(1:min(k,length(resid_est)));
        V(:,1:k) = V(:,1:oldk)*Vtilde;
        if reset
            reset = 0; restarts = 0; rcf = 1;
            allrun = inf(p.maxBasis,1);
            [V(:,1:k),~] = qr(V(:,1:k),0);
            U(:,1:k) = p.A(V(:,1:k),notransp);
            [U(:,1:k),R(1:k,1:k)] = qr(U(:,1:k),0);
            mvp = mvp + k;
        else
            %Has issues if numOld > oldk - k
            [Utilde,R] = restartU(ur,sr,vr,order,p.minRestart,yold);
            U(:,1:k) = U(:,1:oldk)*Utilde;
        end
        allrun(p.minRestart+1:end) = inf;
    end
end

if mvp >= p.maxMV || toc(starttime) >= p.maxTime
    [ur,sr,vr] = svd(R(1:k,1:k));
    U = U(:,1:k)*ur;
    V = V(:,1:k)*vr;
    S = sr;
end


if strcmp(transp,'notransp')
    Temp = U;
    U = V;
    V = Temp;
end

UD = p.userdata;

end

%% Standard Residual Stopping Criteria
function [done,numVals,ud] = resid_stop(numVals,sd,ud,A,U,V)
done = 0;
if sd.k > numVals
    r = sd.rn(1:numVals);
    if all(r < sd.normA*sd.tol)
        done = 1;
    end
end
end

%% Targeting Functions:
%       GKD will expand with as many vectors as elements in returned index
%       In general, index should be no larger than the block size

%% Target based on residual tolerance
%   Returns at most one block size of indices
%   Index will be empty if all ||r|| are below ||A||*tol and k >= numVals
function [index,ud] = resid_target(sd,ud,A,U,V)
r = sd.rn(1:sd.k);
index = find(r > sd.normA*sd.tol);
end

%% Target based on a progressively more restrictive residual tolerance
function [index,ud] = prog_tol_target(sd,ud,A,U,V)
if sd.iters == 1
    index = 1:sd.k;
    return
end

if sd.numVals <= sd.b
    % Always target the first unconverged indices
    index = [1:sd.k]';
    conv = find(sd.rn(1:sd.k) < sd.normA*sd.tol);
    index = setdiff(index,conv,'stable');
else
    % Find converged indices 
    r = sd.rn(1:sd.k);
    conv = find(r < sd.normA*sd.tol);
    
    % Sort 1:nv residuals by order of magnitude (groups of 2)
    r = r(1:min(sd.k,sd.numVals));
    r(r > sd.normA) = sd.normA;
    r_order = abs(ceil(log10(r./sd.normA)./2));
    [~,index] = sort(r_order);
    
    % Fill index with remaining values
    index = [index; [sd.numVals+1:sd.k]'];
    
    % Remove converged indices
    index = setdiff(index,conv,'stable');
end

end

function [index,ud] = large_target(sd,ud,A,U,V)
if sd.iters == 1
    index = 1:sd.k;
    return;
end

r = sd.rn(1:sd.k);

conv = find(r < sd.normA*sd.tol);

r = r(1:min(sd.k,sd.numVals));
[~,index] = sort(r,'descend');

index = [index; [sd.numVals+1:sd.k]'];

index = setdiff(index,conv,'stable');

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
