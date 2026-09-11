# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 13:17:18 2026

@author: apoll
"""

import numpy as np
from numpy.linalg import norm

from AL_ni import AL_ni
from LU import LU
from T_ex2 import T_ex2


def epsL_ni(dh, Ne, m, hp, l):
    
    T_0 = 100;
    T_inf = 20;
   
    # Definition du système numerique =====
    A1 = AL_ni(dh, Ne, m, hp, l); # Definition matrice A
    
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
        T[k] = T_ex2(dh*k);
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