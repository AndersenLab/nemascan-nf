import numpy as np
from scipy import linalg

#############################################################################
###
###  Supporting Functions
###
#############################################################################

def mat_vec(the_vec, nr):
    """Turn a vector into a matrix by repeating it row-wise"""
    return np.tile(the_vec, (nr, 1))

def col_center(the_mat):
    """Center a matrix by column means"""
    n = the_mat.shape[0]
    return the_mat - np.mean(the_mat, axis=0, keepdims=True)

def col_sd(the_mat):
    """Calculate the standard deviation for columns in a matrix"""
    return np.std(the_mat, axis=0, ddof=1)

def col_norm(the_mat):
    """Normalize columns of a matrix"""
    return col_center(the_mat) / col_sd(the_mat).reshape(1, -1)

#############################################################################
###
###  medTest
###
#############################################################################

def medTest(E, M, Y, Z=None, useWeightsZ=True, nperm=100, w=1):
    """
    Mediation test function
    
    INPUT:
    Y:              Outcome (n x 1 array)
    E:              Exposure (n x 1 array)
    M:              Mediator (n x p array, where p is the number of mediators)
    Z:              Additional covariates (n x c array)
    useWeightsZ:    Whether to use weights for Z
    nperm:          Number of permutations for estimating p-value
    w:              Weight assigned to each subject (scalar or array of length n)
    
    OUTPUT:
    Sp:             Statistic and p-value in a p x 2 matrix
    """
    
    # Convert inputs to numpy arrays
    E = np.array(E).reshape(-1, 1)
    Y = np.array(Y).reshape(-1, 1)
    M = np.array(M)
    
    # If M is 1D, reshape to 2D
    if M.ndim == 1:
        M = M.reshape(-1, 1)
    
    n = E.shape[0]
    p = M.shape[1]
    
    # Handle weights
    if np.isscalar(w):
        w = np.full(n, w)
    else:
        w = np.array(w)
        if len(w) != n:
            raise ValueError("The length of w must be either 1 or the length of E.")
    
    # Standardize weights
    w = w / np.sum(w)
    
    # If Z is not null, take residuals
    if Z is not None:
        Z = np.array(Z)
        if Z.ndim == 1:
            Z = Z.reshape(-1, 1)
        
        # Add intercept
        Z1 = np.column_stack([np.ones(n), Z])
        
        # Get residuals for Y
        Y = Y - Z1 @ np.linalg.lstsq(Z1, Y, rcond=None)[0]
        
        if useWeightsZ:
            # Weighted least squares for M and E
            for m in range(p):
                W_sqrt = np.sqrt(w)
                Z1_weighted = Z1 * W_sqrt[:, np.newaxis]
                M_weighted = M[:, m] * W_sqrt
                M[:, m] = M[:, m] - Z1 @ np.linalg.lstsq(Z1_weighted, M_weighted, rcond=None)[0]
            
            E_weighted = E.flatten() * W_sqrt
            E = E - Z1 @ np.linalg.lstsq(Z1_weighted, E_weighted, rcond=None)[0]
        else:
            # Ordinary least squares for M and E
            for m in range(p):
                M[:, m] = M[:, m] - Z1 @ np.linalg.lstsq(Z1, M[:, m], rcond=None)[0]
            E = E - Z1 @ np.linalg.lstsq(Z1, E, rcond=None)[0]
        
        # Ensure E and Y are 2D
        E = E.reshape(-1, 1)
        Y = Y.reshape(-1, 1)
    
    # Normalization
    En = E - np.nanmean(E)
    Mn = col_center(M)
    Yn = Y - np.nanmean(Y)
    sdE = np.nanstd(En, ddof=1)
    sdY = np.nanstd(Yn, ddof=1)
    
    # Getting the residuals
    tEn = En.T
    invCrossEn = 1 / np.nansum(En**2)
    
    B1_obs = invCrossEn * np.nansum(En * Yn)
    B2_obs = invCrossEn * tEn @ Mn
    rY = Yn - B1_obs * En
    rM = Mn - En @ B2_obs
    
    # Calculate these only once
    col_norm_Mn = col_norm(Mn)
    col_norm_rM = col_norm(rM)
    
    # Getting the observed value of the statistic
    # Use weighted E-M correlation
    En_weight = np.sqrt(w).reshape(-1, 1) * (E - np.nansum(w.reshape(-1, 1) * E))
    Mn_weight = np.sqrt(w).reshape(-1, 1) * (M - mat_vec(np.nansum(w.reshape(-1, 1) * M, axis=0), nr=n))
    sdE_weight = np.sqrt(np.nansum(En_weight**2))
    sdMn_weight = mat_vec(np.sqrt(np.nansum(Mn_weight**2, axis=0)), nr=n)
    
    cEM = (En_weight / sdE_weight).T @ (Mn_weight / sdMn_weight)
    cYM = col_norm(rY).T @ col_norm_rM / (n-1)
    
    S = np.abs(cEM * cYM)
    nmed = rM.shape[1]
    
    # Getting p-values
    # Identify the metabolites in Group A
    groupA = np.abs(cEM) >= np.abs(cYM)
    groupB = ~groupA
    
    # Max(S) from nperm permutations
    max_S_mat = np.zeros((nperm, 1))
    
    # Only calculate this value once
    EEiE = invCrossEn * En
    
    for the_perm in range(nperm):
        # Randomize Y
        if the_perm == nperm - 1:  # For the last permutation, use observed data
            rYt = rY
        else:
            rYt = rY[np.random.permutation(n), :]
        
        B1 = np.sum(EEiE * rYt)
        rY2 = rYt - B1 * En
        cYM_A = col_norm(rY2).T @ col_norm_rM / (n-1)
        
        # Randomize E
        if the_perm == nperm - 1:  # For the last permutation, use observed data
            Ent = En
            Et = E
        else:
            permInd = np.random.permutation(n)
            Ent = En[permInd]
            Et = E[permInd, :]
        
        tEnt = Ent.T
        
        B2 = invCrossEn * tEnt @ Mn[:, groupB.flatten()]
        rM3 = rM.copy()
        rM3[:, groupB.flatten()] = Mn[:, groupB.flatten()] - Ent @ B2
        B3 = invCrossEn * np.sum(Ent * Y)
        rY3 = rY - B3 * Ent
        
        # Use weighted E-M correlation
        Ent_weight = np.sqrt(w).reshape(-1, 1) * (Et - np.sum(w.reshape(-1, 1) * Et))
        sdEt_weight = np.sqrt(np.sum(Ent_weight**2))
        
        cEM_B = (Ent_weight / sdEt_weight).T @ (Mn_weight / sdMn_weight)
        cYM_B = cYM
        
        max_S_mat[the_perm] = max(
            np.max(np.abs(cEM_B * cYM_B)[groupB]) if np.any(groupB) else -np.inf,
            np.max(np.abs(cEM * cYM_A)[groupA]) if np.any(groupA) else -np.inf
        )
    
    pval = np.zeros((nmed, 1))
    for i in range(nmed):
        pval[i] = np.sum(max_S_mat > S[0, i]) / nperm
    
    Sp = np.column_stack([S.flatten(), pval.flatten()])
    
    return Sp