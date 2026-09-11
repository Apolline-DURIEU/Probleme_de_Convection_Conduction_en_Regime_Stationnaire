# -*- coding: utf-8 -*-
"""
Created on Wed Dec 31 10:34:10 2025

@author: apoll
"""

# TP1 - C) 

import numpy as np

def LU(A, B):
    
    n = np.shape(A)[0]; # taille de la matrice A (taille nxn) et autres vecteurs
    
    e = 1e-12;
    if(abs(np.linalg.det(A))<e):
        print("erreur la matrice n'est pas inversible");
        return np.zeros(n);
    
    Y = np.zeros(n);
    
    X = np.zeros(n);
    
    # Initialisation de U & L
    U = np.zeros((n,n));
    L = np.eye(n);
    
    # Factorisation
    
    j = 0;
    S=0;
    
    while(j<n): 
        i = 0;
        while(i<(j+1)): # (U)
            k=0;
            S=0;
            while(k<i): # calcul de la somme
                S = S + L[i,k]*U[k,j];
                k=k+1;
            U[i,j]=A[i,j]-S;
            i=i+1;
        while(i<n): # L
            k=0;
            S=0;
            while(k<j): # calcul de la somme
                S = S + L[i,k]*U[k,j];
                k=k+1;
            L[i,j] = 1/U[j,j]*(A[i,j]-S);
            i = i+1;
        j=j+1;
    
    # Resolution
    i=0; #remise a zero des indices
    
    while(i<n): #descente (Y)
        k = 0;
        S = 0;
        while(k<i):
            S = S+L[i,k]*Y[k];
            k=k+1;
        Y[i] = B[i]-S;
        i=i+1;
    
    
    i = n-1;
    while(i>-1): #remonte (X)
        k=i+1;
        S=0;
        while(k<n):
            S = S + X[k]*U[i,k];
            k=k+1;
        X[i] = 1/(U[i,i])*(Y[i]-S);
        i=i-1;

    return X;