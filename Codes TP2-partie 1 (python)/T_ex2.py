# -*- coding: utf-8 -*-
"""
Created on Fri Jan  2 09:53:14 2026

@author: apoll
"""

import numpy as np

def T_ex2(x):
# Fonction de calcul de la solution exacte de la temperature dans le barreau
    
    # Donnees numeriques =======
    T_inf = 20;
    T_0 = 100;
    H = 5*1e-3;
    r = 0.5*1e-3;
    h = 100
    lamb = 40;
    hp = 80;
    # ==========================
    
    m = np.sqrt((h*2)/(lamb*r)); # m = sqtr(hP/(lambda*S)) # Constante 
    
    T_ex = T_inf + (T_0-T_inf)*(np.cosh(m*(H-x))+hp/(lamb*m)*np.sinh(m*(H-x)))/(np.cosh(m*H)+hp/(lamb*m)*np.sinh(m*H)); # Solution Exacte
    
    return T_ex;