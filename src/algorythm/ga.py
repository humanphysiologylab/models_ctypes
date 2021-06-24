import random
import copy
import numpy as np
from numba import njit

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


@njit
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
        child1 = parent1.copy()
        child2 = parent2.copy()
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


def cauchy_mutation_population(population, bounds, gamma, mutation_rate, inplace=False):
    n_genes = len(population[0]['genes'])
    assert (n_genes == len(bounds))

    genes = np.concatenate([organism['genes'].flatten() for organism in population])

    p = np.random.random(len(population))
    shifts = gamma * np.tan(np.pi * (p - 0.5))

    mut_mask = np.random.random(len(population)) <= mutation_rate
    shifts = shifts * mut_mask

    shifts = np.tile(shifts, (n_genes, 1)).T.flatten()

    u = np.random.randn(n_genes * len(population)).reshape((n_genes, len(population)))
    u = u / np.linalg.norm(u, axis=1)[:, None]
    u = u.flatten()

    shifts = shifts * u
    lb, ub = np.tile(bounds, (len(population), 1)).T

    assert (np.all(lb <= genes) and np.all(genes <= ub)), "genes are outside bounds"

    ptp = ub - lb
    b = ub - genes

    shifts = np.remainder(shifts, 2 * ptp)
    shifts = np.abs(np.abs(shifts - b) - ptp) - (ptp - b)

    genes = genes + shifts
    genes = np.reshape(genes, (len(population), n_genes))

    if inplace:
        mutants = population
    else:
        mutants = copy.deepcopy(population)

    for i in range(len(mutants)):
        mutants[i]['genes'] = genes[i]

    return mutants


def selection(population):
    return tournament_selection(population)


def crossover(genes1, genes2, bounds=None):
    return sbx_crossover(genes1, genes2, bounds=bounds)


def mutation(genes, bounds):
    return cauchy_mutation(genes, bounds=bounds)


def do_step(population, new_size, elite_size, bounds, **kw):

    crossover_rate = kw.get('crossover_rate', 1)
    mutation_rate = kw.get('mutation_rate', 1)
    gamma = kw.get('gamma', 1)

    if new_size is None:
        new_size = len(population)
    new_population = []

    while len(new_population) < new_size - elite_size:
        parent1, parent2 = population[0], population[0]
        while parent1 is parent2:
            parent1 = selection(population)
            parent2 = selection(population)

        if np.random.random() <= crossover_rate:
            offspring_genes = sbx_crossover(parent1['genes'], parent2['genes'], bounds=bounds)
            #  offspring_genes = [uniform_crossover(parent1['genes'], parent2['genes'])]
            for genes in offspring_genes:  # TODO: rewrite for different crossover types
                child = dict(genes=genes)
                child['state'] = parent1['state']
                new_population.append(copy.deepcopy(child))
        else:  # no crossover
            child1 = copy.deepcopy(parent1)
            child2 = copy.deepcopy(parent2)
            new_population += [child1, child2]

    new_population = new_population[:new_size - elite_size]  # TODO: sbx_crossover safe processing

    new_population = cauchy_mutation_population(new_population, bounds=bounds, gamma=gamma,
                                                mutation_rate=mutation_rate, inplace=True)

    new_population += sorted(population, key=lambda organism: organism['fitness'], reverse=True)[:elite_size]

    return new_population

@njit
def transform_genes_bounds(genes, bounds, gammas, mask_multipliers):
    assert (len(genes) == len(bounds) == len(gammas) == len(mask_multipliers))

    genes_transformed = np.zeros_like(genes)
    bounds_transformed = np.zeros_like(bounds)

    scaler_dimensional = 1 / np.sqrt(len(genes))
    for i in range(len(genes)):
        lb, ub = bounds[i]
        gene = genes[i]
        if mask_multipliers[i]:  # log10 scale
            bounds_transformed[i, 1] = np.log10(ub / lb) * 1 / (gammas[i] / scaler_dimensional)
            genes_transformed[i] = np.log10(gene)
            lb_temp = np.log10(lb)
            ub_temp = np.log10(ub)
        else:  # linear scale
            genes_transformed[i] = gene
            bounds_transformed[i, 1] = (ub - lb) * 1 / (gammas[i] / scaler_dimensional)
            lb_temp = lb
            ub_temp = ub
        genes_transformed[i] = (genes_transformed[i] - lb_temp) / (ub_temp - lb_temp) * bounds_transformed[i, 1]

    return genes_transformed, bounds_transformed


@njit
def transform_genes_bounds_back(genes_transformed, bounds_transformed, bounds_back, mask_multipliers):
    assert (len(genes_transformed) == len(bounds_transformed) == len(mask_multipliers))

    genes_back = np.zeros_like(genes_transformed)

    for i in range(len(genes_transformed)):  # log10 scale
        lb_back, ub_back = bounds_back[i]
        lb_tran, ub_tran = bounds_transformed[i]
        gene = genes_transformed[i]
        if mask_multipliers[i]:
            genes_back[i] = np.log10(lb_back) + (gene - lb_tran) / (ub_tran - lb_tran) * (
                        np.log10(ub_back) - np.log10(lb_back))
            genes_back[i] = np.power(10, genes_back[i])
        else:  # linear scale
            genes_back[i] = lb_back + (gene - lb_tran) / (ub_tran - lb_tran) * (ub_back - lb_back)

    return genes_back