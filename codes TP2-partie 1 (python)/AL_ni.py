# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 13:05:35 2026

@author: apoll
"""

import numpy as np

def AL_ni(dh, Ne, m, hp, l): # A1 la sous matrice, A la matrice finale (creee en sommant des sous blocs de A1)
# Ne le nombre d'elements   
     
    n = 2;
    
    # Definition de la sous-matrice A1 ============
    K = np.array([[1,-1], [-1, 1]]) * 1/dh;
    M = np.array([[2,1], [1,2]]) * dh/6;
    
    A1 = np.zeros((n,n));
    A1 = K+M*m**2;
    # =============================================

    A = np.zeros((Ne,Ne));
    B = np.zeros((Ne,Ne)); # matrice intermediaire
    i = 0;
    
    while(i<(Ne+1-n)): #i<4
        #redef de B====
        j = i+n;
        ki=i;
        B = np.zeros((Ne,Ne));
        while(ki<j):
            kj=i;
            while(kj<j):
                B[ki,kj] = A1[ki-i,kj-i]; 
                kj = kj+1; 
            ki=ki+1;
        #==============
        A = A + B;
        i = i+1;
  
    #conditions limites
    A[0,1] = 0;
    A[0,0] = 1;
    A[Ne-1, Ne-1] = A[Ne-1, Ne-1] + hp/l;
    return (A);