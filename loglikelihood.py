# coding=utf-8
import networkx as nx
import math
import numpy as np
import operator


def optimisation_likelihood(graphe, constantes, theta_init=None, subgroup_of_nodes=None):
    """On maximise la likelihood, en la maximisant pour chaque i separement
    theta_init est la configuration initiale de l'optimisation
    subgroup_of_nodes est le sous groupe pour lesquels on maximise la likelihood, les autres restent inchanges"""
    kappa0, gamma, beta, Nobs, k_bar = constantes

    kappa = {node: max(kappa0, degree - gamma / beta) for node, degree in graphe.degree().items()}
    mu = beta * math.sin(math.pi / beta) / (2 * math.pi * k_bar)
    a = (nx.to_numpy_matrix(graphe)).A

    if theta_init is None:
        theta_current = 2 * math.pi * np.random.randint(0, Nobs, Nobs) / Nobs
    else:
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

    if subgroup_of_nodes is None: subgroup_of_nodes = graphe.nodes()

    best_theta = theta_init
    best_likelihood = verifyloglikelihood(theta_init, graphe, constantes, silent=False)
    nb_retours = 0
    while (nb_retours < 2):
        for node, degree in sorted(graphe.degree_iter(subgroup_of_nodes), key=operator.itemgetter(1), reverse=True):
            theta_new = maximize_locally(node)
            if degree == 2:
                neighbor, neighbor_degree = max(graphe.degree_iter(graphe.neighbors(node)), key=operator.itemgetter(1))
                theta_new2 = theta_current[neighbor]
                print 2, theta_new, theta_new2
            if degree == 1:
                theta_new1 = theta_current[graphe.neighbors(node)[0]]
                print 1, theta_new, theta_new1

            theta_current[node] = theta_new
        current_likelihood = verifyloglikelihood(theta_current, graphe, constantes, silent=False)
        if current_likelihood < best_likelihood:
            nb_retours += 1
        else:
            best_likelihood = current_likelihood
            best_theta = theta_current
            nb_retours = 0
    return best_theta

'''
def optimisation_complete_par_etapes_likelihood(graphe, bornes, theta0, constantes):
    """On maximise la likelihood en commencant par des sous graphes contenant les noeuds principaux,
    on optimise sur le sous graphe complet à chaque étape"""
    kappa0, gamma, beta, Nobs, k_bar = constantes

    def optimisation_etape_likelihood(nodes_subgraph, theta_global):
        """On maximise la likelihood sur le sous graphe"""
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

        theta_init = np.asarray([theta_global[mapping_new_old[node]] for node in graphe_intermediaire])
        theta_local_opt = optimisation_likelihood(graphe_intermediaire, local_variables, theta_init=theta_init)

        theta_global_opt = theta_global

        for node in graphe_intermediaire:
            theta_global_opt[mapping_new_old[node]] = theta_local_opt[node]
        return theta_global_opt

    theta_current = theta0
    for borne in bornes:
        nodes_subgraph = [node for node, degree in
                          sorted(graphe.degree_iter(), key=operator.itemgetter(1), reverse=True) if
                          degree >= borne]

        theta_current = optimisation_etape_likelihood(nodes_subgraph, theta_current)

    return theta_current
'''

def subgraph_mapping(graphe, nodes_subgraph):
    mapping_old_new = {}
    mapping_new_old = {}
    for i in range(len(nodes_subgraph)):
        old = nodes_subgraph[i]
        mapping_old_new[old] = i
        mapping_new_old[i] = old
        # new network
    graphe_intermediaire = graphe.subgraph(nodes_subgraph)
    graphe_intermediaire = nx.relabel_nodes(graphe_intermediaire, mapping_old_new)
    return graphe_intermediaire, mapping_old_new, mapping_new_old


def optimisation_par_etapes_likelihood(graphe, bornes, borne_lim, theta0, constantes, method="degree"):
    """On maximise la likelihood en commencant par des sous graphes contenant les noeuds principaux,
    on optimise uniquement sur els nouveaux noeuds à chaque étape"""
    kappa0, gamma, beta, Nobs, k_bar = constantes

    def optimisation_etape_likelihood(nodes_subgraph, nodes_evaluated, theta_global):
        # rename nodes of new network
        graphe_intermediaire, mapping_old_new, mapping_new_old = subgraph_mapping(graphe, nodes_subgraph)
        # new "constants"
        nodes_evaluated_local = map(lambda node: mapping_old_new[node], nodes_evaluated)
        Nlocal = len(nodes_subgraph)
        local_variables = kappa0, gamma, beta, Nlocal, k_bar
        # optimize on subgraph
        theta_init = np.asarray([theta_global[mapping_new_old[node]] for node in graphe_intermediaire])
        theta_local_opt = optimisation_likelihood(graphe_intermediaire, local_variables,
                                                  theta_init=theta_init,
                                                  subgroup_of_nodes=nodes_evaluated_local)
        # update theta values
        theta_global_opt = theta_global

        for node in graphe_intermediaire:
            theta_global_opt[mapping_new_old[node]] = theta_local_opt[node]
        return theta_global_opt
    '''
    if method == "degree":
        theta_current = theta0
        borne_sup = float("inf")
        degree_dict = graphe.degree_iter()
        for borne_inf in bornes:
            nodes_subgraph = [node for node, degree in degree_dict if degree >= borne_inf]
            if borne_inf > borne_lim:
                """on optimise sur tous les noeuds du sous graphe"""
                theta_current = optimisation_etape_likelihood(nodes_subgraph, nodes_subgraph, theta_current)
            else:
                """on optimise uniquement sur les noeuds qui viennent juste d'etre ajoutés"""
                nodes_evaluated = [node for node, degree in degree_dict if borne_sup > degree >= borne_inf]
                theta_current = optimisation_etape_likelihood(nodes_subgraph, nodes_evaluated, theta_current)
            borne_sup = borne_inf
    '''
    if method == "core_number":
        theta_current = theta0

        core_dict = nx.core_number(graphe)
        bornes = sorted(set(core_dict.values()), reverse=True)
        borne_sup = float("inf")
        theta_precedent_centre = [theta_current[node] for node, core_num in core_dict.iteritems() if core_num == 27]
        theta_courant_centre = [theta_current[node] for node, core_num in core_dict.iteritems() if core_num == 27]
        for borne_inf in bornes:
            nodes_subgraph = [node for node, core_number in core_dict.iteritems() if core_number >= borne_inf]
            if borne_inf > borne_lim:
                """on optimise sur tous les noeuds du sous graphe"""
                theta_current = optimisation_etape_likelihood(nodes_subgraph, nodes_subgraph, theta_current)
                theta_courant_centre = [theta_current[node] for node, core_num in core_dict.iteritems() if
                                        core_num == 27]
                print borne_inf, np.mean(np.abs(
                    math.pi - np.absolute(np.asarray(theta_courant_centre) - np.asarray(theta_precedent_centre)))), \
                    np.max(np.abs(
                        math.pi - np.absolute(np.asarray(theta_courant_centre) - np.asarray(theta_precedent_centre))))
                theta_precedent_centre = theta_courant_centre
            else:
                """on optimise uniquement sur les noeuds qui viennent juste d'etre ajoutés"""
                nodes_evaluated = [node for node, core_number in core_dict.iteritems() if
                                   borne_sup > core_number >= borne_inf]
                theta_current = optimisation_etape_likelihood(nodes_subgraph, nodes_evaluated, theta_current)
            borne_sup = borne_inf

    return theta_current
    """On maximise la likelihood en commencant par des sous graphes contenant les noeuds principaux,
    on optimise sur le sous graphe complet à chaque étape"""


def verifyloglikelihood(theta_vector, graphe, constantes, silent=False):
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
    if not silent: print np.sum(lp + l_p), np.sum(lp), np.sum(l_p)
    return np.sum(lp + l_p)


def rayon(graphe, constantes):
    kappa0, gamma, beta, Nobs, k_bar = constantes
    kappa = np.asarray([max(kappa0, degree - gamma / beta) for node, degree in graphe.degree().items()])
    mu = beta * math.sin(math.pi / beta) / (2 * math.pi * k_bar)
    rayon = 2 * np.log(Nobs / (math.pi * mu * kappa * kappa0))
    return rayon
