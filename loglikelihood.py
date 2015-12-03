import networkx as nx
import math
import numpy as np
import operator


def optimisation_likelihood(graphe, constantes,theta_init =None,subgroup_of_nodes=None):
    """On maximise la likelihood, en la maximisant pour chaque i separement"""
    kappa0, gamma, beta, Nobs, k_bar = constantes

    kappa = {node: max(kappa0, degree - gamma / beta) for node, degree in graphe.degree().items()}
    mu = beta * math.sin(math.pi / beta) / (2 * math.pi * k_bar)
    a = (nx.to_numpy_matrix(graphe)).A

    if theta_init is None :
        theta_current = 2 * math.pi * np.random.randint(0, Nobs, Nobs) / Nobs
    else :
        theta_current = theta_init


    def localloglikelihood(theta_angle, node):
        """maximisation de la likelihood sur node"""
        th2 = theta_current
        dth = math.pi - np.abs(math.pi - np.abs(theta_angle - th2))
        khi = Nobs * dth / 2 * math.pi * mu * kappa[node] * kappa.values()
        p = 1 / (1 + khi ** beta)

        lp = np.log(p)
        lp[a[node] == 0] = 0
        l_p = np.log(1 - p)
        l_p[a[node] == 1] = 0
        alp = lp + l_p
        return np.sum(alp)

    theta_test = np.linspace(0, 2 * math.pi, num=Nobs)
    def maximize_locally(node):
        """"""
        lll_vector = map(lambda theta_angle: localloglikelihood(theta_angle, node), theta_test)
        theta_optimal = theta_test[np.argmax(lll_vector)]
        return theta_optimal


    if subgroup_of_nodes is None : subgroup_of_nodes = graphe.nodes()
    for node, degree in sorted(graphe.degree_iter(subgroup_of_nodes), key=operator.itemgetter(1), reverse=True):
        theta_new = maximize_locally(node)
        theta_current[node] = theta_new


    return theta_current


def optimisation_complete_par_etapes_likelihood(graphe, bornes, theta0, constantes):

    kappa0, gamma, beta, Nobs, k_bar = constantes

    def optimisation_etape_likelihood(nodes_subgraph,theta_global):
            mapping_old_new = {}
            mapping_new_old = {}
            for i in range(len(nodes_subgraph)):
                old = nodes_subgraph[i]
                mapping_old_new[old] = i
                mapping_new_old[i] = old
            graphe_intermediaire = graphe.subgraph(nodes_subgraph)
            graphe_intermediaire = nx.relabel_nodes(graphe_intermediaire, mapping_old_new)
            Nlocal = len(nodes_subgraph)
            local_variables = kappa0, gamma, beta, Nlocal, k_bar

            theta_init= np.asarray([theta_global[mapping_new_old[node]]for node in graphe_intermediaire])
            theta_local_opt = optimisation_likelihood(graphe_intermediaire,local_variables,theta_init=theta_init)

            theta_global_opt = theta_global

            for node in graphe_intermediaire :
                theta_global_opt[mapping_new_old[node]] = theta_local_opt[node]
            return theta_global_opt

    theta_current = theta0
    for borne in bornes :
        nodes_subgraph = [node for node, degree in sorted(graphe.degree_iter(), key=operator.itemgetter(1), reverse=True) if
                          degree >= borne]

        theta_current = optimisation_etape_likelihood(nodes_subgraph,theta_current)



    return theta_current

def optimisation_par_etapes_likelihood(graphe, bornes,borne_lim, theta0, constantes):

    kappa0, gamma, beta, Nobs, k_bar = constantes

    def optimisation_etape_likelihood(nodes_subgraph, nodes_evaluated,theta_global):
        #rename nodes of new network
            mapping_old_new = {}
            mapping_new_old = {}
            for i in range(len(nodes_subgraph)):
                old = nodes_subgraph[i]
                mapping_old_new[old] = i
                mapping_new_old[i] = old
        #new network
            graphe_intermediaire = graphe.subgraph(nodes_subgraph)
            graphe_intermediaire = nx.relabel_nodes(graphe_intermediaire, mapping_old_new)
        #new "constants"
            nodes_evaluated_local = map(lambda node :mapping_old_new[node],nodes_evaluated)
            Nlocal = len(nodes_subgraph)
            local_variables = kappa0, gamma, beta, Nlocal, k_bar
        #optimize on subgraph
            theta_init= np.asarray([theta_global[mapping_new_old[node]]for node in graphe_intermediaire])
            theta_local_opt = optimisation_likelihood(graphe_intermediaire,local_variables,
                                                      theta_init=theta_init,
                                                      subgroup_of_nodes=nodes_evaluated_local)
        #update theta values
            theta_global_opt = theta_global

            for node in graphe_intermediaire :
                theta_global_opt[mapping_new_old[node]] = theta_local_opt[node]
            return theta_global_opt

    theta_current = theta0
    borne_1 = float("inf")
    for borne in bornes :
        nodes_subgraph = [node for node, degree in sorted(graphe.degree_iter(), key=operator.itemgetter(1), reverse=True) if
                          degree >= borne]
        if borne > borne_lim :
            theta_current = optimisation_etape_likelihood(nodes_subgraph,nodes_subgraph,theta_current)
        else :
            nodes_evaluated = [node for node, degree in sorted(graphe.degree_iter(), key=operator.itemgetter(1), reverse=True) if
                           borne_1 > degree >= borne]
            theta_current = optimisation_etape_likelihood(nodes_subgraph,nodes_evaluated,theta_current)
        borne_1 = borne



    return theta_current



def verifyloglikelihood(theta_vector, graphe, constantes):
    """On calcule somme (aij * logpij + (1-aij)* log(1-pij))"""
    # aij
    a = (nx.to_numpy_matrix(graphe)).A

    # pij
    kappa0, gamma, beta, Nobs, k_bar = constantes
    mu = beta * math.sin(math.pi / beta) / (2 * math.pi * k_bar)
    kappa = {node: max(kappa0, degree - gamma / beta) for node, degree in graphe.degree().items()}
    kp = np.outer(np.ones(Nobs), kappa.values())
    kp2 = np.outer(kappa.values(), np.ones(Nobs))
    mukpkp2 = 2 * math.pi * mu * kp * kp2
    th = np.outer(np.ones(Nobs), theta_vector)
    th2 = np.outer(theta_vector, np.ones(Nobs))
    dth = math.pi - np.abs(math.pi - np.abs(th - th2))
    khi = Nobs * dth / mukpkp2
    p = 1 / (1 + khi ** beta)

    # log(pij)
    lp = np.triu(np.log(p), k=1)
    l_p = np.triu(np.log(1 - p), k=1)

    # aij*log(pij)
    lp[a == 0] = 0
    l_p[a == 1] = 0
    print np.sum(lp + l_p)
    return np.sum(lp + l_p)

def rayon(graphe, constantes):
    kappa0, gamma, beta, Nobs, k_bar = constantes
    kappa = np.asarray([ max(kappa0, degree - gamma / beta) for node, degree in graphe.degree().items()])
    mu = beta * math.sin(math.pi / beta) / (2 * math.pi * k_bar)
    rayon = 2*np.log(Nobs/(math.pi*mu*kappa*kappa0))
    return rayon