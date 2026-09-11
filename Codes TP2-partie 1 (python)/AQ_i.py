# -*- coding: utf-8 -*-
"""
Created on Tue Jan  6 10:54:13 2026

@author: apoll
"""

import numpy as np

def AQ_i(dh, Ne, m): # A1 la sous matrice, AQ la matrice finale (creee en sommant des sous blocs de A1)
# Ne le nombre d'elements   
     
    n = 3;
    N = 2*Ne+1;
    
    # Definition de la matrice raideur K et de celle de masse elementaire M
    K = np.array([[7,-8,1], [-8,16,-8], [1,-8,7]]) * 1/(3*dh);
    M = np.array([[4,2,-1], [2,16,2], [-1,2,4]]) * dh/30;
        
    # Definition de la sous-matrice A1 ============
    A1 = np.zeros((n,n));  
    A1 = K+M*m**2;
    
    # =============================================

    A = np.zeros((N,N));
    B = np.zeros((N,N)); # matrice intermediaire
    i = 0;
    
    while(i<(N+1-n)): 
        #redef de B====
        j = i+n;
        ki=i;
        B = np.zeros((N,N));
        while(ki<j):
            kj=i;
            while(kj<j):
                B[ki,kj] = A1[ki-i,kj-i]; 
                kj = kj+1; 
            ki=ki+1;
        #==============
        A = A + B;
        i = i+2;
  
    #conditions limites
    A[0,:] = 0;
    A[0,0] = 1;
    return (A);