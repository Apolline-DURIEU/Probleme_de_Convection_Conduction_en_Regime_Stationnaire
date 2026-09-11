# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 11:14:07 2026

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

def profil_quadr(Nx='3',Ny='5',graph='false'):
    
    '''Description :
        Fonction de résolution et traçage du profil de température d'une coupe d'ailette via un maillage quadrangulaire
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

    dhx = r/(Nx-1); # pas de discrétisation selon x
    dhy = H/(Ny-1); # pas de discrétisation selon y
    
    Ne = Nx*Ny; # nombre de noeuds
    Nt = (Nx-1)*(Ny-1);
    
    # Initialisation des matrices (certaines sont remplie après) ==============

    KeX = (dhy/(6*dhx))*np.array([[2,-2,-1,1],[-2,2,1,-1],[-1,1,2,-2],[1,-1,-2,2]]);
    KeY = (dhx/(6*dhy))*np.array([[2,1,-1,-2],[1,2,-2,-1],[-1,-2,2,1],[-2,-1,1,2]]);
    Ke = KeX+KeY; # matrice de raideur élémentaire
    K_glob = np.zeros((Ne,Ne)); # Matrice de raideur globale
    
    Mry = (dhy/6)*np.array([[2,1],[1,2]]); # Matrice de masse elementaire pour (x=r)
    MHx = (dhx/6)*np.array([[2,1],[1,2]]); # Matrice de masse elementaire pour (y=H)
    M = np.zeros((Ne,Ne)); # Matrice de masse globale
    
    A = np.zeros((Ne,Ne)); # Matrice Globale
    B = np.zeros(Ne); # Vecteur Global 
    X = np.zeros(Ne); # Vecteur des diff de tempréture (l'inconnu)
    T = np.zeros((Ny,Nx)); # Matrice des températures
    Tp = np.zeros((Ny,Nx)); # Matrice miroir des températures
    # =========================================================================
    
    # Definition du maillage
    Point, Coord = table.tables(r,H,Nx,Ny,'carre'); # Definition des tableaux de Points et de Coordonnées
    
    # Definition de la matrice K_glob =========================================
    k=0; # indice de ligne
    i=0; # indice de colonne de I
    j=0; # indice de colonne de J
    
    while(k<Nt): 
        i = 0;
        while(i<4):
            j=0;
            while(j<4):
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
    
    # Définition du tableau des températures ==================================
    i=0;
    j=0;
    while(j<Ny):
        i=0;
        while(i<Nx):
            T[j,i]=X[i+j*Nx]+T_inf;
            i=i+1;
        j=j+1;
    
    # Définition de la matrice miroir de T
    i=0;
    while(i<Nx):
        Tp[:,i]=T[:,Nx-1-i];
        i = i+1;
    # =========================================================================
    
    # Définition du maillage des coordonées ===================================
    x = np.linspace(0,r,Nx)
    y = np.linspace(0,H,Ny)
    x,y = np.meshgrid(x,y)
    
    xp = np.linspace(-r,0,Nx)
    yp = np.linspace(0,H,Ny)
    xp,yp = np.meshgrid(xp,yp)
    # =========================================================================
    
    # Tracé du profil de température ==========================================
    if(graph=='true'):
        plt.figure(1);
        plt.contourf(x,y,T,levels=30,cmap="RdBu_r")
        plt.contourf(xp,yp,Tp,levels=30,cmap="RdBu_r")
        plt.axis('equal')
        plt.colorbar()
        
        plt.xlabel("x")
        plt.ylabel("y")
        #plt.title("Repartition de la chaleur dans une coupe de l'ailette")
    # =========================================================================
    
    return T;

def evol_quadr(n):
    '''
    Description :
    fonction de traçage de l'évolution de la température en (0,H) en fonction du pas de discrétisation
    n : nombre de points tracés.

    return : H1 vecteur des pas de discrétisation selon x
            TH vecteur des températures selon H1
    -------
    None.

    '''
    p = 0; # indice
    TH = np.zeros(n); # matrice des résultats
    H1 = np.zeros(n);
    while(p<n):
        print(p);
        Nx = 3+p; # nombre de point de discrétisation selon x
        Ny = 2*Nx-1; # nombre de point de discrétisation selon y
    
        T = profil_quadr(Nx, Ny, 'false');
    
        TH[p] = T[Ny-1,0];
        H1[p] =  r/(Nx-1); # = dhx
          
        p=p+1;
    
    # Tracé de la température =================================================
    plt.figure(2);
    plt.plot(H1,TH,'+',color='r')
    plt.plot(H1,TH, '..', color='k')
    
    plt.xlabel("h1")
    plt.ylabel("T(H,0)")
    # =========================================================================
    
    return TH, H1;
        