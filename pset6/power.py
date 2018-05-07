import numpy as np
def normalize(x):
    """Used to normalize vectors."""
    try:
        x = x/np.linalg.norm(x)
    except: 
        x = 0
    return x

def checktol(x,y,tol):
    """Check absolute difference between two values and compare to 
    a defined tolerance. Return boolean. """
    err = abs(x-y)/x
    if err<tol:
        return True
    if err>tol:
        return False

def checkrms(x,y,tol):
    """Calculate RMS error and compare to a defined tolerance. Return 
    boolean. """
    n = len(x)
    rms = ((1/n)*np.sum(((x-y)/y)**2))**(1/2)
    if rms<tol:
        return True
    if rms>tol:
        return False

def power_iteration(hmat,fmat,phi_guess,k_guess,ncells):
    """Power iteration solver for two group diffusion. 
    Solves diffusion problem: H*phi=F*phi

    Paramaters
    ----------
    hmat: numpy.ndarray
        Matrix with dimensions [2*ncells,2*ncells]
    fmat: numpy.ndarray
        Matrix with dimensions [2*ncells,2*ncells]
    phi_guess : numpy.ndarray
        Vector with length 2*ncells, first half is group1, second is group2
    k_guess: float
        Scalar guess of dominant eigenvalue
    ncells: int
        Number of cells in diffusion problem

    Returns
    -------
    numpy.ndarray
        Vector with length 2*ncells with dominant eigenfunction
    float
        Dominant eigenvalue [k_eff]
    int
        Number of iterations to converge
    numpy.ndarray
        Vector with length 2*nells with fission source (1/k)*(F*phi)
    """


    phi_guess = normalize(phi_guess)
    b = (1/k_guess)*np.dot(fmat,phi_guess)
    kold = k_guess
    phiold = phi_guess
    phioldsum = phiold[:ncells] + phiold[ncells:]
    hinv = np.linalg.inv(hmat)
    counter = 0
    while True:
        counter += 1
        # Calculate new flux by solving H*phi_1=b_0
        phinew = np.dot(hinv,b)
        phinewsum = phinew[:ncells]+phinew[ncells:]
        knew = np.sum(np.dot(fmat,phinew))/np.sum(np.dot(fmat,phiold))*kold
        # knew = np.sum(phinew)/np.sum(phiold)*kold
        if checktol(kold,knew,tol=1e-7):
            phi = normalize(phinew)
            b = normalize(b)
            return phi,knew,counter,b
        if checkrms(phioldsum,phinewsum,tol=1e-5):
            phi = normalize(phinew)
            b = normalize(b)
            return phi,knew,counter,b
        kold = knew
        phiold = normalize(phinew)
        phioldsum = phinewsum
        b = (1/kold)*np.dot(fmat,phiold)
