# 1-d Wave Equation Solver
import numpy as np
import matplotlib.pyplot as plt

L=1.0;c=.5;

def f(x):
    return x


A_0=L

def A_n(n,L):
    return (2*L/(n*np.pi)**2)*(-1+(-1)**n)


def B_n(n,L):
    return 0

def omega_n(n,L):
    return n*np.pi/L


# a_list=np.array([an(i,P) for i in n_list])
N_l=200

n_list=np.arange(N_l)+1



def u(x,t,N=N_l):
    uu=np.array([(A_n(n,L)*np.cos(omega_n(n,L)*c*t)+B_n(n,L)*np.sin(omega_n(n,L)*c*t))*np.cos(omega_n(n,L)*x) for n in np.arange(N)+1])
    return A_0/2+np.sum(uu,axis=0)




x_pl=np.linspace(0,L,1000)

fig = plt.figure('Heat Equation')
plt.clf()

ax = fig.gca()


ax.plot(x_pl,f(x_pl))

ax.plot(x_pl,u(x_pl,1.0))

ax.plot


frames=300;
T_0=0.0;T_f=7.0;

tt=np.linspace(T_0,T_f,frames);




plt.show()    








fig2 = plt.figure('Wave Equation Animation')
plt.clf()

ax2 = fig2.gca()


for i in np.arange(tt.size):     
    plt.cla()
        
    ax2.plot(x_pl,f(x_pl))

    ax2.plot(x_pl,u(x_pl,tt[i]))
    

    
    
    ax2.set_xlabel('x')
    # ax.set_xlim(pl*l1, pl*r1)
    ax2.set_ylabel('T')
    # ax.set_ylim(pl*l2, pl*r2)
    
    plt.axis('equal');
            
    plt.savefig('FrameStore/WaveEquation1D/Wave_'+str(i).zfill(6)+'.png',format='png');
    
    print('Done Frame ' + str(i) + '/' + str(tt.size-1))
