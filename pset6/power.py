import numpy as np

def checktol(x,y,tol):
    err = abs(x-y)/x
    print(err)
    if err<tol:
        return True
    if err>tol:
        return False

def checkrms(x,y,tol):
    n = len(x)
    rms = ((1/n)*np.sum(((x-y)/y)**2))**(1/2)
    print(rms)
    if rms<tol:
        return True
    if rms>tol:
        return False

def power_iteration(hmat,fmat,phi_guess,k_guess,ncells):
    b = (1/k_guess)*np.dot(fmat,phi_guess)
    kold = k_guess
    phiold = phi_guess
    phioldsum = phiold[:ncells] + phiold[ncells:]
    hinv = np.linalg.inv(hmat)
    counter = 0
    while True:
        counter += 1
        phinew = np.dot(hinv,b)
        phinewsum = phinew[:ncells]+phinew[ncells:]
        knew = np.sum(np.dot(fmat,phinew))/np.sum(np.dot(fmat,phiold))*kold
        if checktol(kold,knew,tol=1e-7):
            return phinew,knew,counter
        if checkrms(phioldsum,phinewsum,tol=1e-5):
            return phinew,knew,counter
        kold = knew
        phiold = phinew/sum(phinew)
        phioldsum = phinewsum
        b = (1/kold)*np.dot(fmat,phiold)
