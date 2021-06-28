import pickle
import os

import numpy as np
import pandas as pd

from tqdm.auto import tqdm


def collect_results(case, dirname_results, load_dump=False, voigt=False):

    idx = pd.IndexSlice

    if voigt:
        group, cell, suffix = case.split('/')
        group, cell = int(group[-1]), int(cell[-1])
    else:
        group, cell, suffix = None, None, case

    config_path = os.path.join(dirname_results, case)
    with open(os.path.join(config_path, "config_backup.pickle"), 'rb') as f:
        config = pickle.load(f)

    genes = config['runtime']['genes_dict']
    n_genes = sum(map(len, genes.values()))

    n_organisms = config['runtime']['n_organisms']

    dump_filename = os.path.join(config_path, 'dump.bin')
    dump_last_filename = os.path.join(config_path, 'dump_last.npy')
    dump_elite_filename = os.path.join(config_path, 'dump_elite.npy')
    genes_best_filename = os.path.join(config_path, 'genes_best.csv')

    if os.path.isfile(genes_best_filename):
        genes_best = pd.read_csv(genes_best_filename, index_col=[0, 1]).iloc[:, -1]
    else:
        genes_best = None

    if os.path.isfile(dump_last_filename):
        dump_last = np.load(dump_last_filename)
        dump_last = dump_last.reshape(-1, n_genes + 1)
    else:
        dump_last = None

    if os.path.isfile(dump_elite_filename):
        dump_elite = np.fromfile(dump_elite_filename)
        dump_elite = dump_elite.reshape(-1, n_genes + 1)
    else:
        dump_elite = None

    if load_dump and os.path.isfile(dump_filename):
        dump = np.fromfile(dump_filename, dtype=np.half)
        dump = dump.reshape(-1, n_genes + 1)
        n_epochs = len(dump) // n_organisms
    else:
        dump = None
        n_epochs = 1

    columns_tuples = [('common', 'fitness')] + [(item[0], key) for item in genes.items() for key in item[1]]
    m_index_columns = pd.MultiIndex.from_tuples(columns_tuples)

    index_tuples = [(epoch, org) for epoch in range(n_epochs) for org in range(n_organisms)]
    m_index_index = pd.MultiIndex.from_tuples(index_tuples)

    if dump is not None:
        dump = pd.DataFrame(np.roll(dump, axis=1, shift=1),
                            index=m_index_index,
                            columns=m_index_columns)
    if dump_last is None:
        if dump is not None:
            dump_last = dump.loc[idx[n_epochs - 1, :]]
            np.save(dump_last_filename, np.roll(dump_last.values, axis=1, shift=-1))
    else:
        dump_last = pd.DataFrame(np.roll(dump_last, axis=1, shift=1),
                                 columns=m_index_columns)


    phenotype_model_last = {}
    phenotype_model_recreated = {}
    state = {}

    for exp_cond_name in config['experimental_conditions']:

        if exp_cond_name == 'common':
            continue

        filename = os.path.join(config_path, "phenotype", f"phenotype_{exp_cond_name}.csv")
        if os.path.isfile(filename):
            try:
                phenotype_model_last[exp_cond_name] = pd.read_csv(filename)
            except pd.errors.EmptyDataError as e:
                print(f'{filename} is empty')
                continue

        filename_state = filename.replace('phenotype', 'state')
        if os.path.isfile(filename_state):
            state[exp_cond_name] = pd.read_csv(filename_state, index_col=0).iloc[:, -1]


    output_dict = dict(trio = (group, cell, suffix),
                       genes = genes,
                       dump_last = dump_last,
                       dump_elite = dump_elite,
                       dump = dump,
                       phenotype_model_last = phenotype_model_last,
                       config = config,
                       genes_best = genes_best,
                       state=state)

    return output_dict
