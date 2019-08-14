function [ERROR,LBDflag,BlS,normA,tol,maxBasis,maxMVs,v0,display,minRS,numOld,maxII,HLock,noCheck] = optsParse(A,m,n,opts)
ERROR = 0;
optionnames = {'tol','maxBasis','aNorm','maxMV','v0','disp','minRS','numPk',...
    'maxII','seed','locking','LBD','BlockSize','noCheck'};
if isempty(opts)
    names = [];
else
    names = fieldnames(opts);
end
for i = 1:size(names,1)
    x = strcmp(names{i},optionnames);
    if ~x
        disp(strcat("Cannot find option named: ",names{i}));
        ERROR = -1;
    end
end

%% SETUP defaults %%
if isfield(opts,'aNorm')
    normA = opts.aNorm;
else
    if isa(A,'function_handle')
        normA = 1;
        disp('Reverting to non-relative tolerance');
    else
        normA = normest(A);
    end
end

if isfield(opts,'seed')
    rng(opts.seed);
end

if isfield(opts,'LBD')
    LBDflag = opts.LBD;
else
    LBDflag = 0;
end

if isfield(opts,'BlockSize')
    BlS = opts.BlockSize;
else
    BlS = 1;
end

if isfield(opts,'tol')
    tol = opts.tol*normA;
else
    tol = normA*1e-13;
end

if isfield(opts,'maxBasis')
    maxBasis = opts.maxBasis;
else
    maxBasis = 35;
end
maxBasis = floor(maxBasis/BlS)*BlS;

if isfield(opts,'maxMV')
    maxMVs = opts.maxMV;
else
    maxMVs = -1; %Will not stop based on MVs
end

if isfield(opts,'v0')
    v0 = opts.v0;
else
    v0 = randn(n,BlS);
end

if isfield(opts,'disp')
    display = opts.disp;
else
    display = 0;
end

if isfield(opts,'minRS')
    minRS = opts.minRS;
else
    minRS = floor(maxBasis*0.4 + 1);
end

if isfield(opts,'numPk')
    numOld = opts.numPk;
else
    numOld = 1;
end

if isfield(opts,'maxII')
    maxII = opts.maxII;
else
    maxII = 0;
end

if isfield(opts,'locking')
    HLock = opts.locking;
else
    HLock = 0;
end

if isfield(opts,'noCheck')
    noCheck = opts.noCheck;
else
    if tol > 1e-8
        noCheck = 1;
    else
        noCheck = 0;
    end
end
end