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
In order to give the user better control of the algorithm, the following
options can be provided through the opts structure: 

| Field Name      | Description                                             | Default       |
| -------------   | -------------------------------------                   | ------------- |
| tol             |      residual norm tolerance                            | \|\|A\|\|*1e-13 |
| maxBasis        |      max number of basis vectors in V, U                |  35         |
|    aNorm        |      norm(A) estimate                                   |     1    |
|    maxMV        |      maximum number of matrix-vector multiplications    |     Inf    |
|    v0           |      initial vector for V                               | randn(size(A,2),1) |
|    disp         |      options for printing history to console            | 0         |
|    minRS        |      number of vectors to maintain after restart        | floor(0.4*maxBasis + 1)         |
|    numPk        |      +k criteria                                        |    1      |
|    maxII        |      maximum number of inner solver iterations          |     0     |
|    seed         |      Sets seed for the random number generator          |     None     |
|    locking      |      Turns hard locking on if set                       |     0     |
|    isdouble     |      0 = Single Precision, 1 = Double Precision         |     1     |
|    LBD          |      1 = Start with LBD basis up to maxBasis-1          |     0     |
