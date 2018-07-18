

import numpy as np
import argparse
import matplotlib.pyplot as plt
import os
import sys

def square_neighbors(L):
    N = L*L
    site_dic = {}
    x_y_dic = {}
    for j in range(N):
        row = j//L
        column = j-row*L
        site_dic[(row, column)] =j
        x_y_dic[j] = (row,column)
    nbr = []
    for j in range(N):
        row, column = x_y_dic[j]
        right_nbr = site_dic[row, (column+1)%L]
        up_nbr = site_dic[(row+1)%L, column]
        left_nbr = site_dic[row, (column-1+L)%L]
        down_nbr = site_dic[(row-1+L)%L,column]
        nbr.append((right_nbr, up_nbr, left_nbr, down_nbr))
    nbr = tuple(nbr)
    return nbr, site_dic, x_y_dic

def triangular_neighbors(L):
    N = L*L
    site_dic = {}
    x_y_dic = {}
    for j in range(N):
        row = j//L
        column = j-row*L
        site_dic[(row, column)] =j
        x_y_dic[j] = (row,column)
    nbr = []
    for j in range(N):
        row, column = x_y_dic[j]
        right_nbr = site_dic[row, (column+1)%L]
        up_nbr = site_dic[(row+1)%L, column]
        left_nbr = site_dic[row, (column-1)%L]
        down_nbr = site_dic[(row-1+L)%L,column]
        up_left_nbr = site_dic[(row+1+L)%L,(column-1)%L]
        dn_right_nbr = site_dic[(row-1+L)%L,(column+1)%L]
        nbr.append((right_nbr, up_nbr, left_nbr, down_nbr,up_left_nbr,dn_right_nbr))
    nbr = tuple(nbr)
    return nbr, site_dic, x_y_dic


def from_S_to_latt(S,L,site,xy):
    lat = np.zeros((L,L))
    for i in range(len(S)):
        ij = xy[i]
        lat[ij[0],ij[1]] = S[i]
    return lat

def save_latt(T, S,L,site,xy, keyword="Square"):
    lat = from_S_to_latt(S,L,site,xy)
    g=open(keyword+"/"+keyword+"lat"+str(T)+".dat","w")
    for i in range(L):
        for j in range(L):
            g.write("{:} ".format(lat[i,j]))
        g.write("\n")
    g.close()

def check_dir_create(name):
    if(not os.path.isdir("./"+name)):
        os.mkdir(name)


# ** Parse. Get the input for the simulation .

parser = argparse.ArgumentParser(description='Parameters for Monte-Carlo simulation:')
parser.add_argument("-L",dest ="L", help="Number of spin sites on one direction.", type=int, required=True)
parser.add_argument("-g",dest ="Lattice", help="Lattice geometry. Either Square or Triangular", type=str, default="Square")
parser.add_argument("-Ntrial",dest ="Ntrial", help="number of Trial for a Monte-Carlo simulation. Default value =100.",type = int, default = 100)
parser.add_argument("-Tmin",dest ="Tmin", help="Minimal temperature. Default value = 0.1.",type = float,default = 0.1)
parser.add_argument("-Tmax",dest ="Tmax", help="Maximal temperature. Default value = 10.",type=float, default = 10)
parser.add_argument("-nT",dest ="nT", help="number of step in temperature. Default value = 10.",type = int, default = 10)


args = parser.parse_args()

if(args.Lattice != "Square" and args.Lattice!="Triangular"):
    print(args.Lattice)
    print("Not sure about the geometry. Please set the -g parameter to 'Square' or 'Triangular'.")
    sys.exit()


# * Initialization of variables
L = args.L
N= L*L
T_num = np.linspace(args.Tmin,args.Tmax,args.nT)
N_trials = args.Ntrial
if(args.Lattice=="Triangular"):
    nbr, site_dic ,x_y_dic = triangular_neighbors(L)
    keyW = "triangular"
else:
    nbr, site_dic ,x_y_dic = square_neighbors(L)
    keyW = "square"

check_dir_create(keyW)

Magnetization = []

N_cluster_size = []
Mean_cluster_size = []
E_n = []
E_sq_n = []
E=[]
Cv = []
M = np.zeros(args.nT)
Mvar = np.zeros(args.nT)
mag = open(keyW+"/"+keyW+"Magnetization.dat","w")
cl = open(keyW+"/"+keyW+"Cluster_size.dat","w")


for i in range(args.nT):
    t = T_num[i]
    print("temperature = ",t)
    beta = 1./t
    p = 1-np.exp(-2*beta)
    S = [np.random.choice([-1,1]) for _ in range(N)]
    N_cluster_size = np.zeros(N_trials)
    Magnetization = np.zeros(N_trials)
    for itera in range(N_trials):
        k=np.random.randint(0,N)
        Pocket = [k]
        Cluster = [k]
        N_cluster = 1
        while Pocket !=[]:
            s = np.random.choice(Pocket)
            for l in nbr[s]:
                if (S[l] == S[s]):
                    if (l not in Cluster):
                        if (np.random.random_sample()<p):
                            N_cluster +=1
                            Pocket.append(l)
                            Cluster.append(l)
            Pocket.remove(s)
        N_cluster_size[itera] = len(Cluster)
        for s in Cluster:
            S[s] = - S[s]
        Magnetization[itera] =abs(np.sum([S[s] for s in range(N)]))
    save_latt(t,S,L,site_dic,x_y_dic, keyW)
    Mean_cluster_size.append(np.mean(N_cluster_size))
    M[i] = np.mean(Magnetization)/N
    Mvar[i] = np.var(Magnetization/N)/N

    mag.write("{:.10}\t{:.10}\t{:.10}\n".format(t, M[-1], Mvar[-1]))
    cl.write("{:.10}\t{:.10}\n".format(t, Mean_cluster_size[-1]))


mag.close()
cl.close()
# ** Plotting Data

plt.figure()
plt.errorbar(T_num, M, yerr=np.sqrt(Mvar),fmt='x--')
plt.title("Magnetization as a function of temperature for the "+str(L)+"x"+str(L)+" "+keyW+" lattice.")
plt.show()
plt.savefig(keyW+"Magnetization.pdf")
