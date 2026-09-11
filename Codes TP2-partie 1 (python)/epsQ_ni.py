# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 21:45:24 2026

@author: apoll
"""

import numpy as np
from numpy.linalg import norm

from AQ_ni import AQ_ni
from LU import LU
from T_ex2 import T_ex2


def epsQ_ni(H, Ne, m, hp, l):
    
    T_0 = 100;
    T_inf = 20;
   
    N = 2*Ne+1; # taiile des vecteurs et matrices
    dh = H/Ne;
    # Definition du système numerique =====
    A1 = AQ_ni(dh, Ne, m, hp, l); # Definition matrice A
    
    #Définition du vecteur B
    B = np.zeros(N);
    B[0]=T_0-T_inf;
    
    X = LU(A1,B); # Vecteur resultat de theta (diff de temperature)
    
    i=0; # indice
    while(i<N):
        X[i]=X[i]+20; # vecteur resultat de la temperature
        i=i+1;
    #======================================
    
    # Definition du système exacte=========
    T = np.zeros(N);
    k=0;
    dhp = H/(N-1);
    while(k<N):
        T[k] = T_ex2(dhp*k);
        k=k+1;
    # =====================================
      
    # Calcul de l'erreur ==================
    S = 0;
    s = 0;
    i = 0;
    while(i<N):
        S = S + (X[i]-T[i])**2;
        s = s + T[i]**2;
        i = i+1;
    
    eps = np.sqrt(S)/np.sqrt(s);
    #======================================
    
    return eps;