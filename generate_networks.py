import math
import powerlaw
import numpy as np
import networkx as nx


def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step


def clust(N, k_bar, kappa0, gamma, start=1.1, stop=10, step=0.1, nb_networks=100):
    def mapping_function(beta):
        mean, std = clust_beta_given(beta, N, k_bar, kappa0, gamma, nb_networks)
        return [beta, mean, std]

    return map(mapping_function, drange(start, stop, step))


def clust_beta_given(beta, N, k_bar, kappa0, gamma, nb_networks):
    def mapping_function(_):
        return get_clustering_from_random_network(beta, N, k_bar, kappa0, gamma)

    result = map(mapping_function, range(nb_networks))
    return np.mean(result), np.std(result)


def approximated_beta(N, k_bar, kappa0, gamma, c_real, start=1.01, step=0.1, nb_networks=100):
    path_followed = []

    def rec_approximated_beta(beta_min, beta, beta_max):
        if beta_max is not None and (beta_max - beta_min) < step:
            return beta
        c_mean, c_std = clust_beta_given(beta, N, k_bar, kappa0, gamma, nb_networks)
        path_followed.append([beta, c_mean, c_std])
        print beta, c_mean, c_std
        if c_mean > c_real:
            return rec_approximated_beta(beta_min, (beta + beta_min) / 2, beta)
        if beta_max is not None:
            return rec_approximated_beta(beta, (beta + beta_max) / 2, beta_max)
        return rec_approximated_beta(beta, 2 * beta, None)

    beta_optimal = rec_approximated_beta(start, 1 / (1 - c_real), None)
    return beta_optimal, path_followed


def get_clustering_from_random_network(beta, N, k_bar, kappa0, gamma):
    """
    :rtype: float
    """
    graphe = generate_random_network(beta, N, k_bar, kappa0, gamma)
    if all(v == 0 for v in nx.clustering(graphe).values()):
        return 0
    return nx.average_clustering(graphe, count_zeros=False)


def generate_random_network(beta, N, k_bar, kappa0, gamma):
    mu = beta * math.sin(math.pi / beta) / (2 * math.pi * k_bar)
    graphe_random = nx.Graph()
    kappa = powerlaw.Power_Law(xmin=kappa0, parameters=[gamma]).generate_random(N)
    theta = np.random.uniform(0, 2 * math.pi, N)
    # print 1+(len(kappa)/sum(np.log(kappa)))
    th = np.outer(np.ones(N), theta)
    # print "th"
    th2 = np.outer(theta, np.ones(N))
    # print "th2"
    kp = np.outer(np.ones(N), kappa)
    # print "kp"
    kp2 = np.outer(kappa, np.ones(N))
    # print "kp2"
    dth = math.pi - np.abs(math.pi - np.abs(th - th2))
    # print "dth"
    khi = N * dth / (2 * math.pi * mu * kp * kp2)
    # print "khi"
    r = np.random.uniform(0, 1, [N, N])
    # print "r"
    r2 = np.triu(np.log(1 / r - 1) / beta, k=1)
    p2 = np.triu(np.log(khi), k=1)
    # print np.where(r2>p2) ie ln(1/r-1) / beta > ln(khi) <=> r < 1/(1+khi^beta)
    graphe_random.add_edges_from(np.transpose(np.where(r2 > p2)))
    return graphe_random
    # print graphe_random.edges()
    # print graphe_random.number_of_nodes()
    # print graphe_random.number_of_edges()
    # print exponent.get_exponent(graphe_random)
    # print graphe_random.degree().values()
    # print np.mean(graphe_random.degree().values())
    # print max(graphe_random.degree().values())
