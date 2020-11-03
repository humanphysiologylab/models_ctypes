import random
import copy
import numpy as np

from src.helpers import cauchy_inverse_cdf, uniform_vector


def tournament_selection(p, k=2):
    return max(random.sample(p, k), key=lambda x: x['fitness'])


def one_point_crossover(g1, g2):
    assert (len(g1) == len(g2))
    assert (len(g1) >= 2)

    n = len(g1)
    i = random.randint(0, n - 1)

    g_new = list(g1[:i]) + list(g2[i:])
    return g_new


def two_point_crossover(g1, g2):
    assert (len(g1) == len(g2))
    assert (len(g1) >= 3)

    n = len(g1)
    i1 = random.randint(0, n - 2)
    i2 = random.randint(i1, n - 1)

    g_new = list(g1[:i1]) + list(g2[i1: i2]) + list(g1[i2:])
    return g_new


def uniform_crossover(g1, g2):
    assert (len(g1) == len(g2))
    r = np.random.randint(0, 1 + 1, size=len(g1)).astype(bool)
    g_new = np.select([r, ~r], [g1, g2])
    return g_new


def sbx_crossover(parent1, parent2, bounds, cross_rate=0.9):
    # adopted realcross from NSGA-II: Non-dominated Sorting Genetic Algorithm - II
    # Authors: Dr. Kalyanmoy Deb, Sameer Agrawal, Amrit Pratap, T Meyarivan
    # Paper Title: A Fast and Elitist multi-objective Genetic Algorithm: NSGA-II
    # Journal: IEEE Transactions on Evolutionary Computation (IEEE-TEC)
    # Year: 2002
    # Volume: 6
    # Number: 2
    # Pages: 182-197

    assert (len(parent1) == len(parent2))
    child1, child2 = np.empty_like(parent1), np.empty_like(parent1)

    eta_c = 10  # The order of the polynomial for the SBX crossover

    if np.random.random() <= cross_rate:
        for j in range(len(parent1)):

            if np.random.random() <= 0.5:

                y1 = parent1[j]
                y2 = parent2[j]

                if np.abs(y1 - y2) > 1e-14:

                    if y1 > y2:
                        y1, y2 = y2, y1

                    yl, yu = bounds[j]
                    rand = np.random.random()
                    beta = 1.0 + (2.0 * (y1 - yl) / (y2 - y1))
                    alpha = 2.0 - np.power(beta, -(eta_c + 1.0))
                    if rand <= (1.0 / alpha):
                        betaq = np.power((rand * alpha), (1.0 / (eta_c + 1.0)))
                    else:
                        betaq = np.power((1.0 / (2.0 - rand * alpha)), (1.0 / (eta_c + 1.0)))
                    c1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1))
                    beta = 1.0 + (2.0 * (yu - y2) / (y2 - y1))
                    alpha = 2.0 - np.power(beta, -(eta_c + 1.0))
                    if rand <= (1.0 / alpha):
                        betaq = np.power((rand * alpha), (1.0 / (eta_c + 1.0)))
                    else:
                        betaq = np.power((1.0 / (2.0 - rand * alpha)), (1.0 / (eta_c + 1.0)))
                    c2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1))
                    if c1 < yl:
                        c1 = yl
                    if c2 < yl:
                        c2 = yl
                    if c1 > yu:
                        c1 = yu
                    if c2 > yu:
                        c2 = yu
                    if np.random.random() <= 0.5:
                        child1[j] = c2
                        child2[j] = c1
                    else:
                        child1[j] = c1
                        child2[j] = c2
                else:
                    child1[j] = y1
                    child2[j] = y2

            else:
                child1[j] = parent1[j]
                child2[j] = parent2[j]
    else:
        child1 = copy.deepcopy(parent1)
        child2 = copy.deepcopy(parent2)
    return child1, child2


def cauchy_mutation(genes, gamma=1, bounds=None):  # do not change gamma=1

    if bounds is None:
        bounds = [[None, None]] * len(genes)
    assert (len(genes) == len(bounds))

    genes_new = []
    shift = cauchy_inverse_cdf(gamma)
    shift_vec = shift * uniform_vector(len(genes))  # vector mutation

    for gene, (lb, ub), shift in zip(genes, bounds, shift_vec):
        if lb is not None:  # bounce, TODO: rewrite for general case
            assert (lb <= gene <= ub), f"Violated: {lb} <= {gene} <= {ub}"
            b = ub - gene
            ptp = ub - lb
            while shift < 0:
                shift += 2 * ptp
            shift %= (2 * ptp)
            shift = np.abs(np.abs(shift - b) - ptp) - (ptp - b)

        gene_new = gene + shift
        genes_new.append(gene_new)

    return genes_new


def selection(population):
    return tournament_selection(population)


def crossover(genes1, genes2, bounds=None):
    return sbx_crossover(genes1, genes2, bounds=bounds)


def mutation(genes, bounds):
    return cauchy_mutation(genes, bounds=bounds)


def do_step(population, new_size, elite_size, bounds):

    if new_size is None:
        new_size = len(population)
    new_population = []

    if elite_size > 0:
        new_population += sorted(population, key=lambda organism: organism['fitness'], reverse=True)[:elite_size]

    while len(new_population) < new_size:
        parent1, parent2 = population[0], population[0]
        while parent1 is parent2:
            parent1 = selection(population)
            parent2 = selection(population)

        offspring_genes = sbx_crossover(parent1['genes'], parent2['genes'], bounds=bounds)
        #  offspring_genes = [uniform_crossover(parent1['genes'], parent2['genes'])]
        for genes in offspring_genes:  # TODO: rewrite for different crossover types
            child = dict(genes=genes)
            child['genes'] = mutation(child['genes'], bounds=bounds)
            child['state'] = copy.deepcopy(parent1['state'])
            new_population.append(copy.deepcopy(child))

    new_population = new_population[:new_size]  # TODO: sbx_crossover safe processing

    return new_population
