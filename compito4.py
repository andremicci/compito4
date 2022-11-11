'''
 Schema per lo sviluppo del metodo di Numerov
'''


import matplotlib.pyplot as plt 
import numpy as np
import math

#scelgo potenziale elastico

def V(xi):
  return xi*xi


def b(eps, xi):

  return (h**2)*(1/12)*(2*eps-V(xi))





def numerov(n1, n2, eps):
  psi = np.array(xi)*0  # copio xi in psi e lo azzero
  j = np.sign(n2-n1)
  
  psi[n1] = 0
  psi[n1+j] = 1e-06
  for i in range(n1+2*j, n2+j, j):
    psi[i] = ( 2*psi[i-j]*(1-5*b(eps,xi[i-j]))-psi[i-2*j]*(1+b(eps,xi[i-2*j])))/(1+b(eps,xi[i]))
  return psi


def evalDerivative(eps):
    global psir,psil
    psil = numerov(0,nmatch+1,eps)
    psir = numerov(n-1,nmatch-1,eps)
    k=psil[nmatch]/psir[nmatch]
    psir=k*psir
    diff=(psil[nmatch+1]-psil[nmatch-1]-psir[nmatch+1]+psir[nmatch-1])/2*h
   
    return diff

''' 
  Metodo di bisezione per trovare l'energia in cui la funzione evalDerivative
  e' nulla. Cioe' l'enegia per cui  la derivata sinistra e destra coincidono 
'''



def findE(emin,emax,tol):
    while (emax-emin>tol):
        emed = (emin+emax)/2
        if evalDerivative(emin)*evalDerivative(emed)<0:
            emax = emed
        else:
            emin = emed
    return (emin+emax)/2


#compito 4: faccio una funzione per calcolare le richieste per una data funzione d'onda

def dati_psi(psi,e):
  
  mod = psi*psi
  k=np.trapz(mod,xi)
  
  fig, (ax1, ax2) = plt.subplots(2)
  ax1.plot(xi,psi/np.sqrt(k),label=f'$\\varepsilon = {e:.3}$')
  ax1.set_title("Funzione d'onda")
  
  ax2.set_title("Densità di probabilità")
  ax2.plot(xi,mod/k,'tab:red')
  ax1.legend()
  fig.tight_layout()
  
  #calcolo numero di zeri di psi
  n_zeros=0
  for i in range(1,len(psi),1):
    if(psi[i-1]*psi[i]<0):
      n_zeros+=1
  
  #stampo sui grafici
  ax1.text(-7.45,0.62,f"Numero di zeri={n_zeros}",fontsize = 10, 
         bbox = dict(facecolor = 'red', alpha = 0.5))
  
  
  


 
''' 
  Codice principale:
'''
n       = 14000
nmatch  = 10000
xi      = np.linspace(-7.,7,n)
h       = xi[1]-xi[0]



Emax=10  
precisione=0.00001
N=300
E=np.linspace(0,Emax,N)
diff=np.ones(N)

for i in range(0,N,1):
 diff[i]=evalDerivative(E[i])

energie=[]
stati=[]
#trova tutti gli stati fino a Emax e crea un vettore di stati

for i in range(1,len(diff),1):
  if(diff[i]*diff[i-1]<0):
   e=findE(E[i-1],E[i],precisione)
   psil[nmatch+1:n]=psir[nmatch+1:n]
   stati.append(psil)
   energie.append(e)
   
   

m=0
for x in stati:
 dati_psi(x,energie[m])
 m+=1

plt.show() 
