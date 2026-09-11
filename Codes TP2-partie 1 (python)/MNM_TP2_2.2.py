# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 22:11:37 2026

@author: apoll
"""

import numpy as np
from matplotlib import pyplot as plt
from LU import LU
from AQ_ni import AQ_ni
from epsQ_ni import epsQ_ni
from T_ex2 import T_ex2


# 2. Probleme de convection/conduction en régime stationnaire
# 2.2 ÉLÉMENTS FINIS UNIDIMENSIONNELS QUADRATIQUES

#   On cherche à obtenir la répartition de la température dans un barreau métallique (type ailette de refroidissement d’un processeur). Le barreau a une longueur H = 5 mm et le rayon de la section est r = 0.5 mm.

# La température sur l’extrémité du barreau en x = 0 est imposée et égale à T0 = 100°C. A l’autre extrémité, on a placé un matériau isolant (paroi adiabatique, flux thermique nul). L’ailette est refroidie par convection avec le milieu ambiant (h = 100 W.m-2.°C-1, T∞ = 20°C). 

H=5*1e-3; #hauteur du barreau (m)
r=0.5*1e-3; #rayon du bareau (mm)
T_0 = 100; #temperature a x=0 (°C)
h=100; # (W.m^(-2)/°C)
hp = 80;
T_inf = 20; #temperature a l'infini
lambd = 40; #conductivite thermique (W.m^(-1)/°C)

P = np.pi*2*r; # perimetre du barreau (mm^2)
S = np.pi*r**2; # surface de la section du barreau 
m = np.sqrt(h*P/(lambd*S));

Ne = 5; # nombre d'éléments
N = 2*Ne+1; # taille des matrices et vecteurs

dh = H/(Ne);


#   6. Résolution du problème matriciel obtenue à la fin de la question 5
#============================================
# Definition de la matrice AQ
A = AQ_ni(dh, Ne, m, hp, lambd);

print("det(A) = ",np.linalg.det(A));

#Définition du vecteur B
B = np.zeros(N);
B[0]=T_0-T_inf;


print("A = ",A);
print("B = ",B);
#============================================

# Resolution grace a la fonction LU 
X = LU(A,B); #vecteur diff de temperature
print("X = ",X);

#Vecteur temperature
T = np.zeros(N);
i=0; # indice
while(i<N):
    T[i]=X[i]+20; # vecteur resultat de la temperature
    i=i+1;
print("T = ",T);

DH = np.linspace(0,H,N);

x = np.linspace(0,H,50);
dhp = H/(50-1);
Tex = np.zeros(50);
i=0; # indice
while(i<50):
    Tex[i]=T_ex2(i*dhp);
    i=i+1;

plt.figure(1);
plt.plot(x, Tex, color='lightpink');
plt.plot(DH,T, '+', markersize=6, color='r');

plt.legend(("Solution Exacte","Solution numérique"));
plt.xlabel("x (m)");
plt.ylabel("Température (°C)");
plt.title("Méthode Quadratique", fontsize=9);
plt.suptitle("Profil de température dans le barreau (non-isolé en x=H)")
plt.grid();
plt.show;

# ===========================================
# Evaluation de l'erreur en fonction de dh 

t = 50; # tailles des matrices suivantes

E = np.zeros(t); # tableau de collecte de l'erreur en fonction de Ne
DH = np.zeros(t); # tableau de collecte des valeurs de dh
i = 0; #indice

while(Ne<(t+1)):
    print("i= ",i);
    dh = H/(Ne);
    N = 2*Ne+1;
    DH[i] = dh;
    E[i] = epsQ_ni(H, Ne, m, hp, lambd);
    Ne = Ne + 1;
    i = i+1;

print("E = ", E);

# Trace de la figure ========================
plt.figure(2);
plt.loglog(DH,E, '+', color='r', markersize=5);
plt.loglog(DH,E, color='k', linewidth=0.7);

plt.xlabel("dh(m)");
plt.ylabel("erreur relative");
plt.title("Méthode Quadratique", fontsize=9);
plt.suptitle("Erreur Relative en fonction de dh");
plt.grid();
plt.show;
# ===========================================
# ===========================================

#============================================================


