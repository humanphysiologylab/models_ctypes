import sys
import pickle

import numpy as np
import pandas as pd

from scipy.integrate import solve_ivp

from src.helpers import RMSE, calculate_mean_abs_noise, calculate_RMSE_weightened, \
                        calculate_RMSE_balanced, calculate_composite_RMSE_V_CaT, autoscaling, \
                        value_from_bounds

from deprecated import deprecated


def create_genes_dict_from_config(config):
    genes_dict = {ec_name:
                      {p_name: dict(bounds=p['bounds'],
                                    gamma=p['gamma'],
                                    is_multiplier=p.get("is_multiplier", False))
                       for p_name, p in ec['params'].items() if isinstance(p, dict)}
                  for ec_name, ec in config['experimental_conditions'].items()}

    return genes_dict


def create_constants_dict_from_config(config):
    constants_dict = {
        ec_name: {p_name: value for p_name, value in ec['params'].items() if isinstance(value, (int, float))}
        for ec_name, ec in config['experimental_conditions'].items()}

    return constants_dict


def generate_organism(genes_dict, genes_m_index, state):

    genes = [value_from_bounds(gene_params['bounds'], log_scale=gene_params['is_multiplier'])
             if gene_name not in state.index else state[exp_cond_name][gene_name]
             for exp_cond_name, exp_cond in genes_dict.items() for gene_name, gene_params in exp_cond.items()]

    genes = pd.Series(data=genes, index=genes_m_index)

    organism = dict(genes=genes,
                    state=state.copy())
    return organism


def init_population(config):

    population = [generate_organism(config['runtime']['genes_dict'],
                                    config['runtime']['m_index'],
                                    config['runtime']['states_initial'])
                  for _ in range(config['runtime']['n_organisms'])]
    return population


def init_population_from_backup(backup, config):

    genes_dict = config['runtime']['genes_dict']

    assert (len(backup) == config['runtime']['n_organisms'])
    assert (len(backup[0]['genes']) == sum(map(len, genes_dict.values())))
    assert (backup[0]['state'].shape == config['runtime']['states_initial'].shape)

    population = [dict(genes=pd.Series(data=organism['genes'], index=config['runtime']['m_index']),
                       state=pd.DataFrame(data=organism['state'], columns=config['runtime']['states_initial'].columns,
                                          index=config['runtime']['states_initial'].index))
                  for organism in backup]

    return population


def run_model_ctypes(S, C, stim_protocol, config):

    stim_period = C[config['stim_period_legend_name']]
    t_sampling = config['t_sampling']

    if 'n_beats' in config and 't_run' not in config:
        n_beats = config['n_beats']
    elif 'n_beats' not in config and 't_run' in config:
        t_run = config['t_run']
        n_beats = np.ceil(t_run / stim_period).astype(int)
        if n_beats % 2 == 0:
            n_beats += 1
    else:
        print('Invalid config: check n_beats & t_run entries.',
              file=sys.stderr, fflush=True)
        exit()

    tol = config.get('tol', 1e-6)

    n_samples_per_stim = int(np.round(stim_period / t_sampling))
    output = np.empty((n_samples_per_stim * n_beats + 1, len(S)))

    if stim_protocol is not None:
        stim_protocol_Ist = stim_protocol[config['column_stim_protocol']].values.copy()
        if config.get('sparsed_protocol', False):
            stim_protocol_t = stim_protocol['t'].values.copy()
        else:
            stim_protocol_t = None
    else:
        stim_protocol_t = None
        stim_protocol_Ist = None

    if config['run_chain']:
        chain_length = 30
        v_threshold = 1e-1
        t_safe = 10
        status = run_model_ctypes.model.run_chain(S.values.copy(), C.values.copy(),
                                                  chain_length, v_threshold, t_safe,
                                                  n_beats, t_sampling, tol, output)
    else:

        status = run_model_ctypes.model.run(S.values.copy(), C.values.copy(),
                                            n_beats, t_sampling, tol, output,
                                            None, None,
                                            stim_protocol_Ist,  stim_protocol_t
                                            )



    n_beats_save = config.get('n_beats_save', 1)

    output = output[-n_samples_per_stim * n_beats_save - 1:].T

    return status, output


@deprecated(reason="use ctypes models")
def run_model_scipy(S, C, config):

        stim_period = C[config['stim_period_legend_name']]
        t_sampling = config['t_sampling']

        if 'n_beats' in config and 't_run' not in config:
            n_beats = config['n_beats']
        elif 'n_beats' not in config and 't_run' in config:
            t_run = config['t_run']
            n_beats = np.ceil(t_run / stim_period).astype(int)
            if n_beats % 2 == 0:
                n_beats += 1
        else:
            print('Invalid config: check n_beats & t_run entries.',
                  file=sys.stderr, fflush=True)
            exit()

        n_samples_per_stim = int(stim_period / t_sampling)

        t_space = np.linspace(0, stim_period * n_beats, n_samples_per_stim * n_beats + 1, endpoint=True)
        t_span = 0, t_space[-1]

        def fun(t, y, params):  # TODO: rewrite
            ydot = np.zeros_like(y)
            run_model_scipy.model.fun(t, y, ydot, params)
            return ydot

        sol = solve_ivp(fun, y0=S,
                        t_span=t_span, t_eval=t_space,
                        args=(C.values.copy(),),
                        method='LSODA',# rtol=1e-9,
                        max_step=1. * t_sampling,
                )

        output = sol.y[:, -n_samples_per_stim - 1:]
        status = sol.status

        return status, output


def update_S_C_from_genes_current(S, C, genes_current, exp_cond_name, config):

    legend = config['runtime']['legend']
    genes_dict = config['runtime']['genes_dict']
    constants_dict = config['runtime']['constants_dict']
    constants_dict_current = {**constants_dict['common'],
                              **constants_dict[exp_cond_name]}

    for i, g_name in enumerate(genes_current.index.get_level_values('g_name')):

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
        if c_name in legend['constants'].index:
            C[c_name] = c
        if c_name in legend['states'].index:
            S[c_name] = c


def update_genes_from_state(organism, genes_current, config, exp_cond_name):

    legend = config['runtime']['legend']
    genes_dict = config['runtime']['genes_dict']

    for i, g_name in enumerate(genes_current.index.get_level_values('g_name')):

        if g_name in legend['states'].index:
            for ecn in ['common', exp_cond_name]:
                if g_name in genes_dict[ecn]:
                    if genes_dict[ecn][g_name]['is_multiplier']:
                        organism['genes'][ecn, g_name] = organism['state'][exp_cond_name][g_name] \
                                                         / legend['states'][g_name]
                    else:
                        organism['genes'][ecn, g_name] = organism['state'][exp_cond_name][g_name]


def update_phenotype_state(organism, config):

    organism['phenotype'] = dict()

    legend = config['runtime']['legend']

    for exp_cond_name in config['experimental_conditions']:

        if exp_cond_name == 'common':
            continue

        genes_current = organism['genes'][['common', exp_cond_name]]

        S, C = organism['state'][exp_cond_name].copy(), legend['constants'].copy()
        update_S_C_from_genes_current(S, C, genes_current, exp_cond_name, config)

        if config.get('use_scipy', False):
            status, res = run_model_scipy(S, C, config)
            if (status != 0) or np.any(np.isnan(res)):
                return 1
        else:
            stim_protocol = config['experimental_conditions'][exp_cond_name]['stim_protocol']
            status, res = run_model_ctypes(S, C, stim_protocol, config)
            if (status != 2) or np.any(np.isnan(res)):
                return 1

        organism['phenotype'][exp_cond_name] = pd.DataFrame(res.T, columns=legend['states'].index)

        organism['state'][exp_cond_name] = organism['phenotype'][exp_cond_name].iloc[-1]
        update_genes_from_state(organism, genes_current, config, exp_cond_name)

    return 0


def calculate_n_samples_per_stim(exp_cond_name, config):
    stim_period_legend_name = config['stim_period_legend_name']
    stim_period = config['experimental_conditions'][exp_cond_name]['params'][stim_period_legend_name]
    t_sampling = config['t_sampling']
    n_samples_per_stim = int(np.round(stim_period / t_sampling))
    return n_samples_per_stim


def update_fitness(organism, config):

    loss = 0

    columns_control = config.get("columns_control", ["V"])
    columns_model = config.get("columns_model", ["V"])

    if config['loss'] == 'V_CaT_shared':

        phenotype_model_list = []
        phenotype_control_list = []

        for exp_cond_name, exp_cond in config['experimental_conditions'].items():

            if exp_cond_name == 'common':
                continue

            n_samples_per_stim = calculate_n_samples_per_stim(exp_cond_name, config)

            phenotype_control = exp_cond['phenotype'][columns_control].copy()[-n_samples_per_stim - 1:]
            phenotype_model   = organism['phenotype'][exp_cond_name][columns_model].copy()[-n_samples_per_stim - 1:]

            phenotype_model   = phenotype_model[:len(phenotype_control)]

            phenotype_model_list.append(phenotype_model.values)
            phenotype_control_list.append(phenotype_control.values)

        cat_model_concat = np.concatenate([x[:, 1] for x in phenotype_model_list])
        cat_control_concat = np.concatenate([x[:, 1] for x in phenotype_control_list])

        cat_control_concat_scaled, _, (alpha, beta) = autoscaling(signal_to_scale=cat_control_concat,
                                                                  signal_reference=cat_model_concat)

        if alpha <= 0:
            organism['fitness'] = np.NINF
            return

        cumlen = 0
        for i, x in enumerate(phenotype_control_list):

            phenotype_control_list[i][:, 1] = cat_control_concat_scaled[cumlen: cumlen + len(x)]
            cumlen += len(x)

            phenotype_control = phenotype_control_list[i]
            phenotype_model   = phenotype_model_list[i]

            phenotype_control = (phenotype_control - phenotype_model.min(axis=0)) / phenotype_model.ptp(axis=0)
            phenotype_model   = (phenotype_model   - phenotype_model.min(axis=0)) / phenotype_model.ptp(axis=0)

            weights = 1 / calculate_mean_abs_noise(phenotype_control)
            weights /= sum(weights)

            # rmse = calculate_RMSE_weightened(phenotype_control, phenotype_model, weights)
            error = (phenotype_control - phenotype_model)**2
            error[:10] *= 10
            error = np.sqrt(np.mean(error, axis=0))
            error *= weights
            rmse = np.sum(error)

            loss += rmse

        organism['fitness'] = -loss
        return


    for exp_cond_name, exp_cond in config['experimental_conditions'].items():

        if exp_cond_name == 'common':
            continue

        n_samples_per_stim = calculate_n_samples_per_stim(exp_cond_name, config)

        phenotype_control = exp_cond['phenotype'][columns_control][-n_samples_per_stim - 1:]
        phenotype_model   = organism['phenotype'][exp_cond_name][columns_model][-n_samples_per_stim - 1:]

        phenotype_model   = phenotype_model[:len(phenotype_control)]

        if config.get('align_depolarization', False):

            column_v_model = 'V' if 'V' in phenotype_model else 'v'
            column_v_control = 'V' if 'V' in phenotype_control else 'v'

            v_model = phenotype_model[column_v_model].to_numpy()
            v_control = phenotype_control[column_v_control].to_numpy()

            v_level = np.mean(v_model) # 0 # mV
            shift_model = np.where(v_model > v_level)[0][0]
            shift_control = np.where(v_control > v_level)[0]

            if not len(shift_control):
                loss = np.inf
                organism['fitness'] = -loss
                return

            shift_control = shift_control[0]

            shift = shift_model - shift_control

            phenotype_control = np.roll(phenotype_control, shift, axis=1)

        if config['loss'] == 'RMSE':
            loss += RMSE(phenotype_control, phenotype_model)

        elif config['loss'] == 'RMSE_balanced':
            loss += calculate_RMSE_balanced(phenotype_control, phenotype_model)

        elif config['loss'] == 'RMSE_weightened':

            phenotype_model = phenotype_model.to_numpy()
            phenotype_control = phenotype_control.to_numpy()

            phenotype_model   = (phenotype_model   - phenotype_control.min(axis=0)) / phenotype_control.ptp(axis=0)
            phenotype_control = (phenotype_control - phenotype_control.min(axis=0)) / phenotype_control.ptp(axis=0)

            weights = 1 / calculate_mean_abs_noise(phenotype_control)
            weights /= sum(weights)
            loss += calculate_RMSE_weightened(phenotype_control, phenotype_model, weights)

        elif config['loss'] == 'composite_RMSE_V_CaT':
            loss += calculate_composite_RMSE_V_CaT(phenotype_control.to_numpy(),
                                                   phenotype_model.to_numpy())

        elif config['loss'] == 'composite_RMSE_V_CaT_noisy':

            phenotype_model = phenotype_model.to_numpy()
            phenotype_control = phenotype_control.to_numpy()

            phenotype_model   = (phenotype_model   - phenotype_control.min(axis=0)) / phenotype_control.ptp(axis=0)
            phenotype_control = (phenotype_control - phenotype_control.min(axis=0)) / phenotype_control.ptp(axis=0)

            weights = 1 / calculate_mean_abs_noise(phenotype_control)
            weights /= sum(weights)

            rmse_v = RMSE(phenotype_control[:, 0], phenotype_model[:, 0])

            ca_exp_scaled, rmse_ca, coeffs = autoscaling(signal_to_scale=phenotype_control[:, 1],
                                                         signal_reference=phenotype_model[:, 1])

            loss += rmse_v * weights[0] + rmse_ca * weights[1]

        else:
            print(f'Unknown loss {config["loss"]}',
                  file=sys.stderr, fflush=True)
            exit()

    organism['fitness'] = -loss


def generate_bounds_gammas_mask_multipliers(genes_dict):

    bounds, gammas, mask_multipliers = map(np.array, [[x[item] for y in genes_dict.values() for x in y.values()]
                                                      for item in ['bounds', 'gamma', 'is_multiplier']])
    return bounds, gammas, mask_multipliers


def save_epoch(population, kw_output):

    dump_filename = kw_output['dump_filename']
    dump_last_filename = kw_output['dump_last_filename']
    dump_elite_filename = kw_output['dump_elite_filename']
    backup_filename = kw_output['backup_filename']
    log_filename = kw_output['log_filename']

    population = sorted(population, key=lambda organism: organism['fitness'], reverse=True)

    genes = np.array([organism['genes'] for organism in population])
    fitness = np.array([organism['fitness'] for organism in population])

    dump_current = np.hstack([genes, fitness[:, None]])

    with open(dump_filename, 'ba+') as file_dump:
        dump_current.astype(np.half).tofile(file_dump)

    np.save(dump_last_filename, dump_current)

    with open(backup_filename, "wb") as file_backup:
        pickle.dump(population, file_backup)

    with open(log_filename, "a") as file_log:
        s = np.array2string(dump_current[0, :],
                            formatter={'float_kind': lambda x: "%.3f" % x},
                            max_line_width=1000)[1:-1] + "\n"
        file_log.write(s)

    with open(dump_elite_filename, 'ba+') as f:
        dump_current[0, :].tofile(f)
