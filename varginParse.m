function [m,n,transp,notransp,Transpose,numValues,target,opts,Pata] = varginParse(A,inputs)
nextArg = 1;
numIN = size(inputs,2)+1;
transp = 'transp';
notransp = 'notransp';
Transpose = 0;
if isa(A,'function_handle')
    m = inputs{1};
    n = inputs{2};
    nextArg = 3;
    if n > m
        temp = m;
        m = n;
        n = temp;
        transp = 'notransp';
        notransp = 'transp';
        Transpose = 1;
    end
else
    [m,n] = size(A);
    if n > m
        temp = m;
        m = n;
        n = temp;
        transp = 'notransp';
        notransp = 'transp';
        Transpose = 1;
    end
end

if numIN-1 >= nextArg
    numValues = inputs{nextArg};
else
    numValues = 1;
end
nextArg = nextArg + 1;

if numIN-1 >= nextArg
    target = inputs{nextArg};
else
    target = inf;
end
if ischar(target)
    if strcmp(target,'L')
        target = inf;
    elseif strcmp(target,'S')
        target = 0;
    end
end

nextArg = nextArg + 1;

if numIN-1 >= nextArg
    opts = inputs{nextArg};
else
    opts = [];
end
nextArg = nextArg + 1;

if numIN-1 >= nextArg
    Pata = inputs{nextArg};
else
    Pata = 1;
end

end