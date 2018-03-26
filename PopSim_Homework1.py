from __future__ import division
import numpy as np
from math import factorial
import networkx as nx
import csv
import matplotlib.pyplot as plt
import time
from matplotlib.pyplot import pause
from math import *
from math import log
import pandas as pd
import powerlaw


degfile=open('degree.csv','w')
matfile=open('matrix.csv','w')
def read_nodes():
    nodes=[]
    node_file= open('Sorted_services.txt','r')
    for node in node_file:
        temp=node.split(';')
        nodes.append(temp[1])
    return (nodes)

# calculate C(n,k)
def comb(n,k):
    return factorial(n) / factorial(k) / factorial(n - k)

        
# set initial fully connected network as 4 nodes
m0 = 4
init_network_degree = 2*comb(m0,2)

pause_time = 0.1
show_degr = 1
draw_powerlaw = 1
# flag for weather plot the network
plot_network=0

# node position map for plotting
node_pos = {}

# draw the initial network G with random positions
def plot_initial_network(G):
    for i in G.nodes():
       npos = nx.spring_layout(G)

    #Initial plot with only start nodes
    fig = plt.figure('PoP-Sim Network')
    fig.text(0, 0.97, 'starting nodes: green ',style='italic',fontsize=13)
    fig.text(0, 0.94, 'New node: blue ',style='italic',fontsize=13)
    fig.text(0, 0.91, 'Previously added nodes: red', style='italic',fontsize=13)
    nx.draw_networkx(G,npos, node_color = 'green')
    plt.draw()

def plot_new_edges(G, new_edges, i,n):
    plt.clf()
    # Create a color map for nodes to differenctiate between: stating nodes (green), new node (blue) and already added nodes (red)
    color_map = []
    for j in G.nodes():
        if int(j) < m0:
            color_map.append('green')
        elif j == n[m0 + i] :
            color_map.append('blue')
        else: color_map.append('red')
    # Define new node's position and draw the graph
    node_pos= nx.spring_layout(G)
    nx.draw_networkx(G, node_pos, node_color=color_map)
    nx.draw_networkx_edges(G, node_pos,new_edges, width = 2.0 , edge_color = 'b' )
    fig = plt.figure('PoP-Sim Network')
    fig.text(0, 0.97, 'starting nodes: green                 Iteration: '+ str(i+1),style='italic',fontsize=14)
    fig.text(0, 0.94, 'New node: blue ['+str(m0 + i) + ']',style='italic',fontsize=14)
    fig.text(0, 0.91, 'Previously added nodes: red', style='italic',fontsize=14)

    plt.draw()
    pause(pause_time)

def Pop_Sim(m,gamma):
    sim_matrix=gen_matrix()
    n = read_nodes() # read API nodes
    # initialize graph
    G = nx.Graph()
    # Connect all the initial m0 nodes to the graph
    # need to use the nodelist sorted by birth date

    for i in range(m0):
        G.add_node(n[i])
        for j in G:
            if (n[i] != j):
                G.add_edge(n[i], j)

     # could be removed if plotting is not needed
    if plot_network:
        plot_initial_network(G)

    # Adding new nodes
    N=len(sim_matrix)
    for i in range(N-m0):
        # Compute duration of calculations for consistent timing
        loop_start_time = time.time()

       # select neighbors the new node will connect to according to the  hyperbolic distance
        neighbors = choose_neighbour(G,gamma,m,sim_matrix,int(n[i+m0]))
        # A Check to make sure the correct number of neighbors are chosen
        if (len(neighbors) != m):
            print ("Error, number of neighbors is not as expected")
            return
        # Add the new node to the graph
        G.add_node(n[m0 + i])
        # Save new edges in a list for drawing purposed
        new_edges = []
        for nb in neighbors:
            G.add_edge(n[m0 + i], nb)
            new_edges.append((n[m0 + i], nb))
        degfile.write(str(new_edges)+'\n')


        if plot_network:
            plot_new_edges(G, new_edges, i,n)

    plt.close()

    loop_duration = time.time() - loop_start_time
    # Pause for the needed time, taking the calculation time into account
    if pause_time - loop_duration > 0:
        pause(pause_time - loop_duration)

    if draw_powerlaw:
        gen_degree(G)
    #nx.draw(G)

    if show_degr:
        print ('Press any key to continue')
        raw_input()
        powerlaww(G)
    else:
        print ('Press any key to exit')
        raw_input()
def gen_matrix():
    #wk = RWR(ApisList)
    #sim=wk.main()   #run RWR
    rwr_matrix=np.genfromtxt('real_rwr.csv', delimiter=';')
    np.fill_diagonal(rwr_matrix, 0)#zero diagnal of rwr matrix
    r=rwr_matrix
    rwr_matrix=normalize_matrix(rwr_matrix)# normalize rwr matrix
    uniform_matrix=(1/len(rwr_matrix))+np.zeros((len(rwr_matrix),len(rwr_matrix)))
    np.fill_diagonal(uniform_matrix, 0)
    sim_matrix = rwr_matrix + uniform_matrix#  add two matrix
    sim_matrix =( sim_matrix / 2 )
    return(sim_matrix)

def normalize_matrix(matrix):
    row_sums = matrix.sum(axis=1)
    new_matrix = matrix/row_sums[:, np.newaxis]
    return(new_matrix)
         
         
def choose_neighbour(G,gamma,m,sim,t):
    dist=hyperbolic_dist(G,gamma,sim,t)
    print(dist)
    connect_to = np.random.choice(G.nodes(),m, replace=False, p=dist)
    return(connect_to)


def hyperbolic_dist(G,gamma,sim,t):
    dist = []
    beta = 1 / (gamma - 1)
    for node in G.nodes():
        rs=(beta* log(int(node)) + (1 - beta) *(log(t)))
        theta_st=(sim[int(node)-1][t-1] *(pi/2))
        dist.append(abs(1-(rs + log(t)+log(theta_st/2))))
    sum_dist = sum(dist)
    dist = [i / sum_dist for i in dist]
    return(dist)

#save degree in a text file for further analysis 
def gen_degree(G):
    plt.close() 
    degfile = open('degreelist.csv','w')
    deg_list=[]
    for n in G.nodes():
        deg_list.append(G.degree(n))
    degrees=np.array(deg_list)
    print(degrees)
    for d in deg_list:
        degfile.write(str(d)+'\n')
    return(degrees)


#powerlaw fitting and degree exponent.
def powerlaww(G):
    degree_list=gen_degree(G)
    fit=powerlaw.Fit(degree_list,linear_bins=True)
    print("The degree of exponent is:", fit.power_law.alpha)
    print("gamma's statndard error is:",fit.power_law.sigma)
    print("Likelihood of PL and Exponetial:",fit.distribution_compare('power_law','exponential'))
    fig1=fit.plot_pdf(color='r', linewidth=0, marker='*')
    fig2=fit.plot_ccdf(color='b', linewidth=2)
    fit.power_law.plot_pdf(color='b',linestyle='--')
    fit.plot_ccdf(color='r', linewidth=1)
 #   fit.power_law.plot_ccdf(color='r',linestyle='--',ax=fig2)
 #   powerlaw.plot_pdf(rands,linear_bins=True,color='r')
    plt.pause(100)




if __name__ == "__main__":
    #read_nodes()
   Pop_Sim(3,3)
   #powerlaww()

    
