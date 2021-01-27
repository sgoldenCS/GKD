# GKD


A Golub-Kahan Davidson Method for Accurately Computing a Few Singular Triplets


## Usage:
 
#### Outputs:

    [u,s,v,hist] = GKD(...)

    u             Left Singular Vectors
    s             Singular Values
    v             Right Singular Vectors
    hist          Convergence History

#### Inputs:

    [...] = GKD(A,numValues)
    [...] = GKD(A,numValues,...)

    A             m by n matrix or a matrix function
    numValues     maximum number of singular values to return
    
    If A is a function pointer, the user must provide name/value 
        pairs for 'm' and 'n' (the size of the matrix)

#### Name/Value Pair Options for GKD
  
    tol                residual norm tolerance for standard stopping criteria
    maxBasis           max number of basis vectors in V, U  
    normA              norm(A) estimate 
    maxMV              maximum number of matrix-vector multiplications  
    v0                 initial vector for V 
    showHist           options for printing history to console  
    minRS              number of vectors to maintain after restart  
    numOld             +k criteria  
    maxII              maximum number of inner solver iterations (GKJD)
    seed               Sets seed for the random number generator  
    sort_order         Determines whether the largest ('L'), smallest ('S') or values 
                          closest to a target (double) will be returned.                            
    b                  Block Size

    
#### Advanced Name/Value Pair Options
    P                  Preconditioner for AtA (either a matrix or function pointer)
    target_func        Function pointer to determine which values are targeted 
                          for improvement at each step
    stop_func          Function to determine when GKD will terminate

    Both target_func and stop_func must take 5 input arguments (but do not need to use them)
        normA                    Current estimate of the ||A||_2
        k                        Current number of vectors in the basis
        s                        Current estimates of the singular values (sorted by sort_order)
        [u,v] = f_vecs(i)        Function handle to return the i'th singular vectors (u,v)
        r = f_resid(i,compute)   Function handle to return the i'th residual norms (r)
                                 The compute flag (0/1) forces direct computation when set
                                 Otherwise, the most recently computed norms are returned.
                                 These may be out-of-date or 'inf'
    
    [indices] = target_func(...) must return an array of indices for the algorithm to target. 
        If more than b (block size) indices are returned, only the first b indices are targeted
        If fewer than b indices are returned, the block will not be padded
        If no indices are returned, the algorithm terminates and returns 
            min(k,numVals) singular triplets
        
    [done,numVals] = stop_func(...) must return:
        'done'      A flag (0/1) to indicate when to terminate the algorithm 
        'numVals'   The number of values the user would like returned 
    
    Note: The algorithm will not terminate when k < numVals even if the 'done' flag is set. 
