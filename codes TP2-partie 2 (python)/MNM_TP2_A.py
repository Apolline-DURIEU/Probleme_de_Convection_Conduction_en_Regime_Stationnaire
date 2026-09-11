# -*- coding: utf-8 -*-
"""
Created on Fri Jan 23 12:23:25 2026

@author: apoll
"""

import matplotlib.pyplot as plt
import numpy as np
import table
from LU import LU

# Definition des constantes
H = 5*1e-3; # hauteur (selon y) (m)
r = 0.5*1e-3; # largeur = 2r (selon x) (m)

T_0 = 100; # température du système à y=0 (°C)
T_inf = 20; # température du système à l'infini (°C)

h = 100; # coefficient d'échange par convection à x=r (W/(m^2.°C))
hp = 80 # coefficient d'échange par convection à y=H (W/(m^2.°C))
lambd = 40; # conductivité thermique dans le barreau (W/(m.°C))

def profil_triangle(Nx='3',Ny='5',graph='false'):

    '''Description :
        Fonction de résolution et traçage du profil de température d'une coupe d'ailette via un maillage triangulaire
        Renvoie le tableau des températures du maillage et trace le profil
        H : longueur selon y (m)
        r : largueur selon x (m)
        isolée en x=0
        maintenue à T_0 en y=0 (°C)
        h : coefficient de convection en x=r (W/(m^2.°C))
        hp : coefficient de convection en y=H (W/(m^2.°C))
        lambd : conductivité thermique dans le barreau (W/(m.°C))
        T_inf = température de l'air ambiant
        
        Nx : nombre de point de discrétisation selon x
        Ny : nombre de point de discrétisation selon y
    '''
    
    # Définition des constantes ===============================================
    
    dhx = r/(Nx-1); # pas de discrétisation selon x
    dhy = H/(Ny-1); # pas de discrétisation selon y
    
    
    Ne = Nx*Ny; # nombre de noeuds
    Nt = 2*(Nx-1)*(Ny-1);
    
    # Initialisation des matrices (qui sont remplie après)
    
    KeX = np.zeros((3,3)); # Matrice de raideur élémentaire K_x
    KeY = np.zeros((3,3)); # Matrice de raideur élémentaire K_y
    Ke = np.zeros((3,3)); # Matrice de raideur élémentaire globale
    K_glob = np.zeros((Ne,Ne)); # Matrice de raideur globale
    Ap = np.zeros((3,3)); # Matrice de coordonné des points de chaque éléments
    
    Mry = np.zeros((2,2)); # Matrice de masse elementaire pour (x=r)
    MHx = np.zeros((2,2)); # Matrice de masse elementaire pour (y=H)
    M = np.zeros((Ne,Ne)); # Matrice de masse globale
    
    A = np.zeros((Ne,Ne)); # Matrice Globale
    B = np.zeros(Ne); # Vecteur Global 
    X = np.zeros(Ne); # Vecteur des diff de tempréture (l'inconnu)
    T = np.zeros((Ny,Nx)); # Vecteur des températures
    Tp = np.zeros((Ny,Nx));
    # =========================================================================
    
    # Definition du maillage
    Point, Coord = table.tables(r,H,Nx,Ny,'triangle'); # Definition des tableaux de Points et de Coordonnées
    #maillage.maillage(r,H,Nx,Ny,'triangle'); # Tracé du maillage
    
    # Definition des matrices elementaires KeX et KeY
    
    # Definition de la matrice de raideur K_glob ==============================
    k=0; # indice de ligne
    i=0; # indice de colonne de I
    j=0; # indice de colonne de J
    
    while(k<Nt):
        # calcul Ke
        # _remplissage de la matrice intermediaire de coordonnées
        
        Ap = np.array([[1,Coord[Point[k,0],0],Coord[Point[k,0],1]],[1,Coord[Point[k,1],0],Coord[Point[k,1],1]],[1,Coord[Point[k,2],0],Coord[Point[k,2],1]]]);
    
        y12 = Ap[0,2] - Ap[1,2];
        y23 = Ap[1,2] - Ap[2,2];
        y31 = Ap[2,2] - Ap[0,2];
        x12 = Ap[0,1] - Ap[1,1];
        x23 = Ap[1,1] - Ap[2,1];
        x31 = Ap[2,1] - Ap[0,1];
        
        KeX = np.array([[y23**2, y23*y31, y23*y12],[y23*y31, y31**2, y31*y12], [y23*y12, y31*y12, y12**2]]);
        KeY = np.array([[x23**2, x23*x31, x23*x12],[x23*x31, x31**2, x31*x12], [x23*x12, x31*x12, x12**2]]);
    
        Ke = (1/(2*np.linalg.det(Ap)))*(KeX+KeY); # matrice de raideur élémentaire
    
        i = 0;
        while(i<3):
            j=0;
            while(j<3):
                Ig = Point[k,i];
                Jg = Point[k,j];
                K_glob[Ig,Jg] = K_glob[Ig,Jg] + Ke[i,j]; # matrice de raideur globale
                
                j=j+1;
            i=i+1;
        k=k+1;
    # =========================================================================
    
    # Definition des matrices de masse ========================================
    i=0; # indice de colonne
    j=0; # indice de ligne
    
    # ___definition des matrices de masse élémentaires
    MHx = (dhx/6)*np.array([[2,1],[1,2]]);
    Mry = (dhy/6)*np.array([[2,1],[1,2]]);
    
    m=0; # indice du nombre de matrices élémentaires à mettre dans la matrice gloable
    
    while(m<(Nx-1)): # sommage de la matrice selon x
        i=0;
        while(i<2):
            j=0;
            while(j<2):
                M[(Ne-2-m)+i,(Ne-2-m)+j] = M[(Ne-2-m)+i,(Ne-2-m)+j]+(hp/lambd)*MHx[i,j];
                j=j+1;
            i=i+1;    
        m=m+1;
    
    m=0;
    
    while(m<(Ny-1)): # sommage de la matrice selon y
        i=0;
        while(i<2):
            j=0;
            while(j<2):
                M[(Ne-1)+(i-1-m)*Nx,(Ne-1)+(j-1-m)*Nx] = M[(Ne-1)+(i-1-m)*Nx,(Ne-1)+(j-1-m)*Nx]+(h/lambd)*Mry[i,j];
                j=j+1;
            i=i+1;
        m=m+1;
    # =========================================================================
    
    A = K_glob+M; # Matrice Globale
    
    # Conditions initiales ====================================================
    i=0; # indice de colonne
    
    while(i<Nx): # Condition de Dirichlet
        A[i,:]=0;
        A[i,i]=1;
        B[i]=T_0-T_inf;
        i=i+1;
    # =========================================================================
    
    X = LU(A,B); #RESOLUTION !!
    
    # Création du tableau des températures ====================================
    i=0;
    j=0;
    while(j<Ny):
        i=0;
        while(i<Nx):
            T[j,i]=X[i+j*Nx]+T_inf;
            i=i+1;
        j=j+1;
    
    # Création du tableau miroir de T
    i=0;
    while(i<Nx):
        Tp[:,i]=T[:,Nx-1-i];
        i = i+1;
    # =========================================================================
    
    # Création du maillage des coordonnées ====================================
    x = np.linspace(0,r,Nx)
    y = np.linspace(0,H,Ny)
    x,y = np.meshgrid(x,y)
    
    xp = np.linspace(-r,0,Nx)
    yp = np.linspace(0,H,Ny)
    xp,yp = np.meshgrid(xp,yp)
    # =========================================================================
    
    # Tracé de la résolution ==================================================
    if(graph=='true'):
        plt.figure(2);
        plt.contourf(x,y,T,levels=30,cmap="RdBu_r")
        plt.contourf(xp,yp,Tp,levels=30,cmap="RdBu_r")
        plt.axis('equal')
        plt.colorbar()
        
        plt.xlabel("x")
        plt.ylabel("y")
        #plt.title("Repartition de la chaleur dans une coupe de l'ailette")
    # =========================================================================
    
    return T;


