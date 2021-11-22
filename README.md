# GKD


A Golub-Kahan Davidson Method for Accurately Computing a Few Singular Triplets


## Usage:
 
#### Outputs:
 
    [U,S,V,H,UD] = GKD(...);

    U             Left Singular Vectors
    S             Singular Values
    V             Right Singular Vectors
    H             Convergence History
    UD            User Data

#### Inputs:

    [U,S,V,H,UD] = GKD(A,numVals,...)
    
    A             m by n matrix
    numValues     number of singular values to return
    
  Additional Options for GKD (Name, Value pairs)
    
    tol           residual norm tolerance
    SIGMA         Sets which values are desired ('L' or 'S')
    target_fn     Function used for expanding the basis
                       [index,userdata] = target_fn(solverdata,userdata)
    stop_fn       Function used for stopping the solver
                       [done,numVals,userdata] = stop_fn(numVals,solverdata,userdata)
    maxMV         maximum number of matrix-vector multiplications
    maxTime       maximum allowed computed time
    normA         norm(A) estimate
    display       Prints partial history to console if set
    v0            Initial vector for V
    b             Block size
    minRestart    Number of vectors to maintain after restart
    maxBasis      max number of basis vectors in V,U
    numOld        number of +k vectors to keep after restart
    maxII         max number of inner solver iterations
    seed          random seed
    m             number of rows in A (if A is a function_handle)
    n             number of cols in A (if A is a function_handle)
    P             preconditioner for AtA

  Default Options Settings

    tol          1e-6
    SIGMA        'L' (Largest)
    target_fn    Targeting based on first values with residuals above 'tol'
    stop_fn      Stopping based on residual tolerance ('tol')
    maxMV        inf (No stopping based on matvecs)
    maxTime      inf (No stopping based on time)
    normA        Largest value seen (Accurate when SIGMA = 'L')
    display      0 (Off)
    v0           Gaussian random vectors
    b            1
    minRestart        numVals+max(b,15)
    maxBasis     max(minRestart+2*b,floor(1.3*minRestart))
    numOld       1
    maxII        0
    seed         'shuffle' (sets rng based on current time)
    P            1 (Identity matrix)
