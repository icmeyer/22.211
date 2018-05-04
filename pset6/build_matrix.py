import numpy as np
from materials import mat_properties

def dtilde(d1,d2,delta):
    return 2*d1*d2/(delta*(d1+d2))

def build_matrix(problem):
    """
    Parameters
    ----------
    problem : Dictionary
        Dictionary with members that are each a zone as a list in the form
        [Width [cm],Spacing,Material ID]

    Returns
    -------
    numpy.ndarray
        H matrix for performing power iteration
    numpy.ndarray
        F matrix for performing power iteration
    """
    mesh = {'size': [],'materials': [],'spacing': []}
    zones = 0
    mesh['size'] = 0
    for zone in problem:
        width = problem[zone][0]
        spacing = problem[zone][1]
        material = problem[zone][2]
        ncells = int(width/spacing)
        mesh['materials'] = mesh['materials']+[material]*ncells
        mesh['spacing'] = np.concatenate([mesh['spacing'],
                                         np.ones([ncells])*spacing])
        mesh['size'] += ncells
        zones += 1

    ncells = mesh['size']
    matsize = mesh['size'] * 2
    hmat = np.zeros([matsize,matsize])
    fmat = np.zeros([matsize,matsize])

    #Building the first and last line for each group for H matrix
    #First Line
    delta = mesh['spacing'][0]
    nf1 =  mat_properties[mesh['materials'][0]]['NF1']
    nf2 =  mat_properties[mesh['materials'][0]]['NF2']
    sigma_s12 = mat_properties[mesh['materials'][0]]['S12']
    sigma_r11 = (mat_properties[mesh['materials'][0]]['A1'] 
               + mat_properties[mesh['materials'][0]]['S12'])
    sigma_a2 =  mat_properties[mesh['materials'][0]]['A2']
    #d11 is the material property
    d11 = mat_properties[mesh['materials'][0]]['D1']
    d12 = mat_properties[mesh['materials'][1]]['D1']
    D112 = dtilde(d11,d12,delta)
    d21 = mat_properties[mesh['materials'][0]]['D2']
    d22 = mat_properties[mesh['materials'][1]]['D2']
    D212 = dtilde(d21,d22,delta)
    #D11 is the D-tilde used in the matrix
    D11 = 2*d11/delta*(1/(1+4*d11/delta))
    D21 = 2*d21/delta*(1/(1+4*d21/delta))
    hmat[0,0] = sigma_r11*delta+D11+D112
    hmat[0,1] = -D112
    hmat[ncells,ncells] = sigma_a2*delta+D21+D212
    hmat[ncells,ncells+1] = -D212
    hmat[ncells,0] = -sigma_s12*delta
    fmat[0,0] = nf1*delta
    fmat[0,ncells] = nf2*delta

    #Last Line
    delta = mesh['spacing'][-1]
    nf1 =  mat_properties[mesh['materials'][-1]]['NF1']
    nf2 =  mat_properties[mesh['materials'][-1]]['NF2']
    sigma_s12 = mat_properties[mesh['materials'][-1]]['S12']
    sigma_r1N = (mat_properties[mesh['materials'][-1]]['A1'] 
               + mat_properties[mesh['materials'][-1]]['S12'])
    sigma_a2 =  mat_properties[mesh['materials'][-1]]['A2']
    #d1n is the material property
    d1n = mat_properties[mesh['materials'][-2]]['D1']
    d1N = mat_properties[mesh['materials'][-1]]['D1']
    D1nN = dtilde(d1n,d1N,delta)
    d2n = mat_properties[mesh['materials'][-2]]['D2']
    d2N = mat_properties[mesh['materials'][-1]]['D2']
    D2nN = dtilde(d2n,d2N,delta)
    #D1N is the D-tilde used in the matrix
    D1N = 2*d1N/delta*(1/(1+4*d1N/delta))
    D2N = 2*d2N/delta*(1/(1+4*d2N/delta))
    hmat[ncells-1,ncells-2] = -D1nN
    hmat[ncells-1,ncells-1] = sigma_r1N*delta+D1N+D1nN
    hmat[2*ncells-1,2*ncells-2] = -D2nN
    hmat[2*ncells-1,2*ncells-1] = sigma_a2*delta+D2N+D2nN
    hmat[2*ncells-1,ncells-1] = -sigma_s12*delta
    fmat[ncells-1,ncells-1] = nf1*delta
    fmat[ncells-1,2*ncells-1] = nf2*delta

    # Building upper half of matrix for energy group 1
    for i in range(ncells-2):
        cell = i+1
        index = i+1
        delta = mesh['spacing'][cell]
        sigma_s12 =  mat_properties[mesh['materials'][cell]]['S12']
        sigma_a1  =  mat_properties[mesh['materials'][cell]]['A1'] 
        sigma_a2  =  mat_properties[mesh['materials'][cell]]['A2'] 
        sigma_r1n = sigma_a1+sigma_s12
        # Defining d1l d1c d1r as left, center, right
        d1l = mat_properties[mesh['materials'][cell-1]]['D1']
        d1c = mat_properties[mesh['materials'][cell]]['D1']
        d1r = mat_properties[mesh['materials'][cell+1]]['D1']
        D1lc = dtilde(d1l,d1c,delta)
        D1cr = dtilde(d1c,d1r,delta)
        nf1 =  mat_properties[mesh['materials'][cell]]['NF1']
        nf2 =  mat_properties[mesh['materials'][cell]]['NF2']
        hmat[index,index-1] = -D1lc 
        hmat[index,index]   = sigma_r1n*delta + D1lc + D1cr
        hmat[index,index+1] = -D1cr
        fmat[index,index] = nf1*delta
        fmat[index,index+ncells] = nf2*delta

    # Building lower half of matrix for energy group 2
    for i in range(ncells-2):
        cell = i+1
        index = ncells+i+1
        delta = mesh['spacing'][cell]
        sigma_s12 =  mat_properties[mesh['materials'][cell]]['S12']
        sigma_a1  =  mat_properties[mesh['materials'][cell]]['A1'] 
        sigma_a2  =  mat_properties[mesh['materials'][cell]]['A2'] 
        sigma_r1n = sigma_a1+sigma_s12
        # Defining d1l d1c d1r as left, center, right
        d2l = mat_properties[mesh['materials'][cell-1]]['D2']
        d2c = mat_properties[mesh['materials'][cell]]['D2']
        d2r = mat_properties[mesh['materials'][cell+1]]['D2']
        D2lc = dtilde(d2l,d2c,delta)
        D2cr = dtilde(d2c,d2r,delta)
        nf1 =  mat_properties[mesh['materials'][cell]]['NF1']
        nf2 =  mat_properties[mesh['materials'][cell]]['NF2']
        hmat[index,index-1] = -D2lc 
        hmat[index,index]   = sigma_a2*delta + D2lc + D2cr
        hmat[index,cell] = -sigma_s12*delta
        hmat[index,index+1] = -D2cr

    return hmat, fmat, ncells

        

        

