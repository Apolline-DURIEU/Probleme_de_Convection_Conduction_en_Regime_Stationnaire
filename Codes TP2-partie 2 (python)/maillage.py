# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 12:36:15 2026

@author: apoll
"""

import matplotlib.pyplot as plt
import numpy as np
import table
import trace

def maillage(Lx=10,Ly=20,Nx=3,Ny=5,e='triangle'):
    """Parameters
    ----------
    Lx : TYPE, optional
        DESCRIPTION. The default is 10.
    Ly : TYPE, optional
        DESCRIPTION. The default is 20.
    Nx : TYPE, optional
        DESCRIPTION. The default is 3.
    Ny : TYPE, optional
        DESCRIPTION. The default is 5.
    e : TYPE, optional
        DESCRIPTION. The default is 'triangle'.
    """
    """
        Fonction de tracé d'un maillage triangulaire, se basant sur un tableau de coordonnées et un tableau de points. 
        Tableaux créés par la fonction table, tracé de chaques trianges individuellement tracés par la fonction trace
    """
    
    
    if(e=='carre'):
        print("erreur : le maillage tracé ne peut pas être quadratique");
    else :
        Tbc,Noeud = table.tables(Lx,Ly,Nx,Ny,e); # création des tableaux de points et de coordonnées
        nT = np.shape(Tbc)[0]; # nombre de triangle
        CT = np.zeros((3,2)); # Tableau intermédiaire des coordonnées de chaque triangle
        
        plt.figure()
        
        k=0; #indice
        while(k<nT): # pour chaque triangle existant
            i=0;
            while(i<3): # pour chacun des trois points du triangle
                p = Tbc[k,i]; # défini quel est le point
                CT[i,:] = Noeud[p,:]; #met les coordonnés du point dans le tableaux de coordonnées du triangle
                i = i+1;
            trace.trace(CT); # rajoute le triangle au tracé du maillage
            k=k+1;
            
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title('Tracé du maillage')
        plt.grid()
        plt.axis('equal')
        plt.show()
        
    return ;