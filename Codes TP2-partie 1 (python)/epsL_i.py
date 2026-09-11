# -*- coding: utf-8 -*-
"""
Created on Fri Jan  2 09:43:37 2026

@author: apoll
"""

import numpy as np
from numpy.linalg import norm

from AL_i import AL_i
from LU import LU
from T_ex1 import T_ex1


def epsL_i(dh, Ne, m):
    
    T_0 = 100;
    T_inf = 20;
   
    # Definition du système numerique =====
    A1 = AL_i(dh, Ne, m); # Definition matrice A
    
    #Définition du vecteur B
    B = np.zeros(Ne);
    B[0]=T_0-T_inf;
    
    X = LU(A1,B); # Vecteur resultat de theta (diff de temperature)
    
    i=0; # indice
    while(i<Ne):
        X[i]=X[i]+20; # vecteur resultat de la temperature
        i=i+1;
    #======================================
    
    # Definition du système exacte=========
    T = np.zeros(Ne);
    k=0;
    while(k<Ne):
        T[k] = T_ex1(dh*k);
        k=k+1;
    # =====================================
      
    # Calcul de l'erreur ==================
    S = 0;
    s = 0;
    i = 0;
    while(i<Ne):
        S = S + (X[i]-T[i])**2;
        s = s + T[i]**2;
        i = i+1;
    
    eps = np.sqrt(S)/np.sqrt(s);
    #======================================
    
    return eps;