import pickle
import os

import numpy as np
import pandas as pd

from tqdm.auto import tqdm


def collect_results(case, dirname_results, load_dump=False):

    idx = pd.IndexSlice

    group, cell, suffix = case.split('/')
    group, cell = int(group[-1]), int(cell[-1])

    config_path = os.path.join(dirname_results, case)
    with open(os.path.join(config_path, "config_backup.pickle"), 'rb') as f:
        config = pickle.load(f)

    genes = config['runtime']['genes_dict']
    n_genes = sum(map(len, genes.values()))

    n_organisms = config['runtime']['n_organisms']

    dump_filename = os.path.join(config_path, 'dump.bin')
    dump_last_filename = os.path.join(config_path, 'dump_last.npy')

    if os.path.isfile(dump_last_filename):
        dump_last = np.load(dump_last_filename)
        dump_last = dump_last.reshape(-1, n_genes + 1)
    else:
        dump_last = None

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
        dump_last = dump.loc[idx[n_epochs - 1, :]]
        np.save(dump_last_filename, np.roll(dump_last.values, axis=1, shift=-1))
    else:
        dump_last = pd.DataFrame(np.roll(dump_last, axis=1, shift=1),
                                 columns=m_index_columns)


    phenotype_model_last = {}

    for exp_cond_name in tqdm(config['experimental_conditions'], desc='phenotype'):

        if exp_cond_name == 'common':
            continue

        filename = os.path.join(config_path, "phenotype", f"phenotype_{exp_cond_name}.csv")
        if os.path.isfile(filename):
            try:
                phenotype_model_last[exp_cond_name] = pd.read_csv(filename)
            except pd.errors.EmptyDataError as e:
                print(f'{filename} is empty')
                continue

    output_dict = dict(trio = (group, cell, suffix),
                       genes = genes,
                       dump_last = dump_last,
                       dump = dump,
                       phenotype_model_last = phenotype_model_last,
                       config = config)

    return output_dict


def create_C_S(organism, config, exp_cond_name):

    legend = config['runtime']['legend']
    genes_dict = config['runtime']['genes_dict']
    constants_dict = config['runtime']['constants_dict']

    genes_current = organism['genes'][['common', exp_cond_name]]
    constants_dict_current = {**constants_dict['common'],
                              **constants_dict[exp_cond_name],
                             }

    C = legend['constants'].copy()
    S = organism['state'][exp_cond_name].copy()

    for i in range(len(genes_current)):
        g_name = genes_current.index.get_level_values(1).to_list()[i]

        if g_name in legend['constants'].index:
            for ecn in ['common', exp_cond_name]:
                if g_name in genes_dict[ecn]:
                    if genes_dict[ecn][g_name]['is_multiplier']:
                        C[g_name] *= genes_current[ecn, g_name]
                    else:
                        C[g_name] = genes_current[ecn, g_name]

        if g_name in legend['states'].index:
            for ecn in ['common', exp_cond_name]:
                if g_name in genes_dict[ecn]:
                    if genes_dict[ecn][g_name]['is_multiplier']:
                        S[g_name] *= genes_current[ecn, g_name]
                    else:
                        S[g_name] = genes_current[ecn, g_name]

    for c_name, c in constants_dict_current.items():
        print(c_name, c)
        if c_name in legend['constants'].index:
            C[c_name] = c
        if c_name in legend['states'].index:
            S[c_name] = c

    return C, S
