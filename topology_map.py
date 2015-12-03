import networkx as nx
import exponent
import equations
import numpy as np
import generate_networks
import math
import loglikelihood
import warnings
from matplotlib import pyplot as plt

warnings.filterwarnings("ignore")

def main():
    np.random.seed(50)
    r1,t1 = study("drug_net.txt")
    print t1[:5]
    r2,t2 = study("drug_net.txt",theta_init = t1)
    print t2[:5]
    graphe = nx.read_weighted_edgelist("drug_net.txt")
    graphe = nx.convert_node_labels_to_integers(graphe)
    degree = graphe.degree()
    tdiff = math.pi - np.abs(math.pi - np.abs(t1 - t2))
    degres = [1,2,3,4,5,6]

    diffavg = [np.average([tdiff[node] for node in degree if degree[node]==deg]) for deg in degres]
    plt.plot(range(len(tdiff)),tdiff,degree.values(),'o')
    plt.show()





    #Etude 2 :

def study(name,theta_init = None) :
    # graphe = nx.read_edgelist("Internet_hyperbolic_map/AS_ARK200906_topology.txt")
    graphe = nx.read_weighted_edgelist(name)
    gamma = exponent.get_exponent(graphe)
    # gamma = 2.1
    # Verify gamma :
    # exponent.plot_gamma(graphe)

    Nobs = graphe.number_of_nodes()
    k_bar_obs = np.mean(graphe.degree().values())
    k_max_obs = max(graphe.degree().values())
    C_obs = nx.average_clustering(graphe, count_zeros=False)
    print Nobs, k_bar_obs, k_max_obs, C_obs

    alpha, kappa0, kappaC, P0, k_bar = equations.solve(gamma, k_bar_obs, k_max_obs)
    # Verify equations
    #equations.verify(alpha, kappa0, kappaC, P0, k_bar, gamma, k_bar_obs, k_max_obs)
    N = Nobs / (1 - P0)
    print gamma, alpha, kappa0, kappaC, P0, k_bar, N
    '''
    alpha = 0.58
    kappa0 = 0.9
    kappaC = 4790
    P0 = 0.33
    N = 35685
    k_bar = 9.86
    '''
    # Estimating beta
    '''
    beta,path_followed = generate_networks.approximated_beta(N,k_bar,kappa0,gamma,C_obs,
                                             start=1.01,  step = 0.01,nb_networks=100)
    print beta
    '''
    # Verify beta:
    # print path_followed

    beta = 1.02
    # Creating map:
    graphe = nx.convert_node_labels_to_integers(graphe)
    constants = kappa0,gamma,beta,Nobs,k_bar
    if theta_init is None : theta_init = 2*math.pi*np.random.randint(0,Nobs,Nobs)/Nobs
    theta_opt = loglikelihood.optimisation_likelihood(graphe,constants,theta_init=theta_init.copy())
    #theta_opt = loglikelihood.optimisation_complete_par_etapes_likelihood(graphe,[5,4,3,2,0],theta_init.copy(),constants)
    r_opt = loglikelihood.rayon(graphe,constants)
    #theta_opt = loglikelihood.optimisation_par_etapes_likelihood(graphe,[6,5,4,3,1],3,theta_init.copy(),constants)

    return r_opt,theta_opt

    #Verify loglikelihood
    #loglikelihood.verifyloglikelihood(theta_init,graphe,constants)
    #loglikelihood.verifyloglikelihood(theta_opt,graphe,constants)


main()
