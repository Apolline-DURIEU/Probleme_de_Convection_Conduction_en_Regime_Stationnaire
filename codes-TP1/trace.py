# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 12:04:25 2026

@author: apoll
"""
import numpy as np
import matplotlib.pyplot as plt

def trace(C): 
    """Fonction de traçage d'un triangle qui prend comme argument le tableau C 
        de coordonnées x,y de trois points 
    """
    # C le tableau de coordonnées du triangle
    if(np.shape(C)[0]!=3 or np.shape(C)[1]!=2):
        print("erreur : le tableau de coordonnees n'est pas aux bonnes dimensions")
    else :
        #Séparation des X et des Y
        X = np.array([C[0,0], C[1,0], C[2,0], C[0,0]]);
        Y = np.array([C[0,1], C[1,1], C[2,1], C[0,1]]);
        
        #trace
        plt.plot(X, Y, '-', color='red') 
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title('Tracé du triangle')
        plt.grid()
        plt.axis('equal')
    
    return;