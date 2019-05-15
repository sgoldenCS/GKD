# GKD


A Golub-Kahan Davidson Method for Accurately Computing a Few Singular Triplets


## Usage:
 
#### Outputs:
 
    [s] = GKD(...)
    [s,r] = GKD(...)
    [u,s,v] = GKD(...)
    [u,s,v,r] = GKD(...)
    [u,s,v,r,stats] = GKD(...)
    [u,s,v,r,stats,hist] = GKD(...)


    u             Left Singular Vectors
    s             Singular Values
    v             Right Singular Vectors
    r             norm(A'*u - s*v) for each singular triplet
    stats         Matvecs, Time and estimated norm(A)
    hist          Convergence History

#### Inputs:

    [...] = GKD(A,numValues,target)
    [...] = GKD(A,numValues,target,opts)
    [...] = GKD(A,numValues,target,opts,Pata)

    A             m by n matrix
    numValues     number of singular values to return
    target        seek singular value nearest target
                      if target is 'L', GKD computes the largest singular values
                      if target is 'S', GKD computes the smallest singular values
    opts          structure containing extra solver options
    Pata          preconditioner for A'*A

#### Options for GKD
  
    opts.tol                residual norm tolerance  
    opts.maxBasis           max number of basis vectors in V, U  
    opts.aNorm              norm(A) estimate 
    opts.maxMV              maximum number of matrix-vector multiplications  
    opts.v0                 initial vector for V 
    opts.disp               options for printing history to console  
    opts.minRS              number of vectors to maintain after restart  
    opts.numPk              +k criteria  
    opts.maxII              maximum number of inner solver iterations  
    opts.seed               Sets seed for the random number generator  
    opts.locking            Turns hard locking on if set 
    opts.isdouble           0 = Single Precision, 1 = Double Precision 
    opts.LBD                1 = Start with LBD basis up to maxBasis-1  
    opts.AtQ                1 = Keep extra storage for AtQ. Reduces Matvecs  
                              with locking = 0  
