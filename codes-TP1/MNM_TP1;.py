# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import table
import trace
import maillage
from LU import LU

"""
Created on Wed Dec  3 14:10:11 2025

@author: apoll
"""

# A. CREATION DE MATRICES PARTICULIERES ================
# A.1)
def A1(n): # fonction matrice pentadiagonale
    a1 = 4*np.ones(n*n);
    b1 = -1*np.ones(n*n-1);
    c1 = 2*np.ones(n*(n-1));
    A1 = np.diag(a1) + np.diag(b1,1) + np.diag(b1,-1) + np.diag(c1,n) + np.diag(c1,-n);
    return A1;

# A.2)
def A2(n): # n est le nombe de lignes et colonnes dont on décale la matrice A avant de la sommer 
    a2 = np.ones(3);
    A = np.ones((3,3))+np.diag(a2);
    B = np.zeros((9,9));
    A2 = np.zeros((9,9));
    i = 0;
    while(i<6):
        B[i:i+3,i:i+3]=A;
        A2 = A2 + B;
        B = np.zeros((9,9));
        i = i+n;
    return A2;

#B. TRACE DE GRAPHIQUES ===============================
def Nuage(n):
    X = np.random.random(n);
    Y = np.random.random(n);
    
    plt.plot(X,Y, 'o')

    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend("Mes points")
    plt.title('Nuage de point')
    plt.grid()
    plt.show()
    return;
    
def NuageLogLog(n):
    X = np.random.random(n);
    Y = np.random.random(n);
    
    plt.loglog(X,Y, 'o')

    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend("Mes points")
    plt.title('Nuage de point en échelle logarithmique')
    plt.grid()
    plt.show()
    return;

#C. PROGRAMMATION DE LA METHODE INVERSION DE SYSTEME LINEAIRE PAR FACTORISATION LU ==
# voir fonction LU exporté (LU.py)
  
# D. PRINCIPE DE LA STRUCTURE DES DONNEES POR ELEMENTS FINIS ========================
# D.1) table de connexion :
    
Connex,Coord = table.tables(10,20,3,5,'carre'); #test
# table de connexion & table de coordonnées pour maillage structuré avec nœuds disposés en quadrillage tel que (Nx = 3) et (Ny = 5)

C = np.array([[2,3],[4,1],[4,7]]); #test
trace.trace(C); #test

# D.2) Tracé du maillage :
# cef maillage.py
maillage.maillage(); #test

