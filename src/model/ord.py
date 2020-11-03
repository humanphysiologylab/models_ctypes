import numpy as np
import pandas as pd
from numba import njit

size_algebraic = 200
size_states = 49
size_constants = 206

legend_states = pd.DataFrame(
    [{'name': 'v', 'component': 'membrane', 'units': 'millivolt'},
    {'name': 'CaMKt', 'component': 'CaMK', 'units': 'millimolar'},
    {'name': 'cass', 'component': 'intracellular_ions', 'units': 'millimolar'},
    {'name': 'nai', 'component': 'intracellular_ions', 'units': 'millimolar'},
    {'name': 'nass', 'component': 'intracellular_ions', 'units': 'millimolar'},
    {'name': 'ki', 'component': 'intracellular_ions', 'units': 'millimolar'},
    {'name': 'kss', 'component': 'intracellular_ions', 'units': 'millimolar'},
    {'name': 'cansr', 'component': 'intracellular_ions', 'units': 'millimolar'},
    {'name': 'cajsr', 'component': 'intracellular_ions', 'units': 'millimolar'},
    {'name': 'cai', 'component': 'intracellular_ions', 'units': 'millimolar'},
    {'name': 'm', 'component': 'INa', 'units': 'dimensionless'},
    {'name': 'hf', 'component': 'INa', 'units': 'dimensionless'},
    {'name': 'hs', 'component': 'INa', 'units': 'dimensionless'},
    {'name': 'j', 'component': 'INa', 'units': 'dimensionless'},
    {'name': 'hsp', 'component': 'INa', 'units': 'dimensionless'},
    {'name': 'jp', 'component': 'INa', 'units': 'dimensionless'},
    {'name': 'mL', 'component': 'INaL', 'units': 'dimensionless'},
    {'name': 'hL', 'component': 'INaL', 'units': 'dimensionless'},
    {'name': 'hLp', 'component': 'INaL', 'units': 'dimensionless'},
    {'name': 'a', 'component': 'Ito', 'units': 'dimensionless'},
    {'name': 'iF', 'component': 'Ito', 'units': 'dimensionless'},
    {'name': 'iS', 'component': 'Ito', 'units': 'dimensionless'},
    {'name': 'ap', 'component': 'Ito', 'units': 'dimensionless'},
    {'name': 'iFp', 'component': 'Ito', 'units': 'dimensionless'},
    {'name': 'iSp', 'component': 'Ito', 'units': 'dimensionless'},
    {'name': 'd', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'ff', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'fs', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'fcaf', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'fcas', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'jca', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'ffp', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'fcafp', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'nca', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'IC1', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'IC2', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'C1', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'C2', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'O', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'IO', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'IObound', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'Obound', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'Cbound', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'D', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'xs1', 'component': 'IKs', 'units': 'dimensionless'},
    {'name': 'xs2', 'component': 'IKs', 'units': 'dimensionless'},
    {'name': 'xk1', 'component': 'IK1', 'units': 'dimensionless'},
    {'name': 'Jrelnp', 'component': 'ryr', 'units': 'dimensionless'},
    {'name': 'Jrelp', 'component': 'ryr', 'units': 'dimensionless'}])
    
    
legend_constants = pd.DataFrame(
    [{'name': 'celltype', 'component': 'environment', 'units': 'dimensionless'},
    {'name': 'nao', 'component': 'extracellular', 'units': 'millimolar'},
    {'name': 'cao', 'component': 'extracellular', 'units': 'millimolar'},
    {'name': 'ko', 'component': 'extracellular', 'units': 'millimolar'},
    {'name': 'R', 'component': 'physical_constants', 'units': 'joule_per_kilomole_kelvin'},
    {'name': 'T', 'component': 'physical_constants', 'units': 'kelvin'},
    {'name': 'F', 'component': 'physical_constants', 'units': 'coulomb_per_mole'},
    {'name': 'zna', 'component': 'physical_constants', 'units': 'dimensionless'},
    {'name': 'zca', 'component': 'physical_constants', 'units': 'dimensionless'},
    {'name': 'zk', 'component': 'physical_constants', 'units': 'dimensionless'},
    {'name': 'L', 'component': 'cell_geometry', 'units': 'centimeter'},
    {'name': 'rad', 'component': 'cell_geometry', 'units': 'centimeter'},
    {'name': 'i_Stim_Start', 'component': 'membrane', 'units': 'millisecond'},
    {'name': 'i_Stim_End', 'component': 'membrane', 'units': 'millisecond'},
    {'name': 'i_Stim_Amplitude', 'component': 'membrane', 'units': 'microA_per_microF'},
    {'name': 'i_Stim_Period', 'component': 'membrane', 'units': 'millisecond'},
    {'name': 'i_Stim_PulseDuration', 'component': 'membrane', 'units': 'millisecond'},
    {'name': 'KmCaMK', 'component': 'CaMK', 'units': 'millimolar'},
    {'name': 'aCaMK', 'component': 'CaMK', 'units': 'per_millimolar_per_millisecond'},
    {'name': 'bCaMK', 'component': 'CaMK', 'units': 'per_millisecond'},
    {'name': 'CaMKo', 'component': 'CaMK', 'units': 'dimensionless'},
    {'name': 'KmCaM', 'component': 'CaMK', 'units': 'millimolar'},
    {'name': 'cmdnmax_b', 'component': 'intracellular_ions', 'units': 'millimolar'},
    {'name': 'kmcmdn', 'component': 'intracellular_ions', 'units': 'millimolar'},
    {'name': 'trpnmax', 'component': 'intracellular_ions', 'units': 'millimolar'},
    {'name': 'kmtrpn', 'component': 'intracellular_ions', 'units': 'millimolar'},
    {'name': 'BSRmax', 'component': 'intracellular_ions', 'units': 'millimolar'},
    {'name': 'KmBSR', 'component': 'intracellular_ions', 'units': 'millimolar'},
    {'name': 'BSLmax', 'component': 'intracellular_ions', 'units': 'millimolar'},
    {'name': 'KmBSL', 'component': 'intracellular_ions', 'units': 'millimolar'},
    {'name': 'csqnmax', 'component': 'intracellular_ions', 'units': 'millimolar'},
    {'name': 'kmcsqn', 'component': 'intracellular_ions', 'units': 'millimolar'},
    {'name': 'cm', 'component': 'intracellular_ions', 'units': 'microF_per_centimeter_squared'},
    {'name': 'PKNa', 'component': 'reversal_potentials', 'units': 'dimensionless'},
    {'name': 'mssV1', 'component': 'INa', 'units': 'millivolt'},
    {'name': 'mssV2', 'component': 'INa', 'units': 'millivolt'},
    {'name': 'mtV1', 'component': 'INa', 'units': 'millivolt'},
    {'name': 'mtV2', 'component': 'INa', 'units': 'millivolt'},
    {'name': 'mtD1', 'component': 'INa', 'units': 'dimensionless'},
    {'name': 'mtD2', 'component': 'INa', 'units': 'dimensionless'},
    {'name': 'mtV3', 'component': 'INa', 'units': 'millivolt'},
    {'name': 'mtV4', 'component': 'INa', 'units': 'millivolt'},
    {'name': 'hssV1', 'component': 'INa', 'units': 'millivolt'},
    {'name': 'hssV2', 'component': 'INa', 'units': 'millivolt'},
    {'name': 'Ahf', 'component': 'INa', 'units': 'dimensionless'},
    {'name': 'GNa', 'component': 'INa', 'units': 'milliS_per_microF'},
    {'name': 'shift_INa_inact', 'component': 'INa', 'units': 'millivolt'},
    {'name': 'thL', 'component': 'INaL', 'units': 'millisecond'},
    {'name': 'GNaL_b', 'component': 'INaL', 'units': 'milliS_per_microF'},
    {'name': 'Gto_b', 'component': 'Ito', 'units': 'milliS_per_microF'},
    {'name': 'Kmn', 'component': 'ICaL', 'units': 'millimolar'},
    {'name': 'k2n', 'component': 'ICaL', 'units': 'per_millisecond'},
    {'name': 'PCa_b', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'GKr_b', 'component': 'IKr', 'units': 'milliS_per_microF'},
    {'name': 'A1', 'component': 'IKr', 'units': 'per_millisecond'},
    {'name': 'B1', 'component': 'IKr', 'units': 'per_millivolt'},
    {'name': 'q1', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'A2', 'component': 'IKr', 'units': 'per_millisecond'},
    {'name': 'B2', 'component': 'IKr', 'units': 'per_millivolt'},
    {'name': 'q2', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'A3', 'component': 'IKr', 'units': 'per_millisecond'},
    {'name': 'B3', 'component': 'IKr', 'units': 'per_millivolt'},
    {'name': 'q3', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'A4', 'component': 'IKr', 'units': 'per_millisecond'},
    {'name': 'B4', 'component': 'IKr', 'units': 'per_millivolt'},
    {'name': 'q4', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'A11', 'component': 'IKr', 'units': 'per_millisecond'},
    {'name': 'B11', 'component': 'IKr', 'units': 'per_millivolt'},
    {'name': 'q11', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'A21', 'component': 'IKr', 'units': 'per_millisecond'},
    {'name': 'B21', 'component': 'IKr', 'units': 'per_millivolt'},
    {'name': 'q21', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'A31', 'component': 'IKr', 'units': 'per_millisecond'},
    {'name': 'B31', 'component': 'IKr', 'units': 'per_millivolt'},
    {'name': 'q31', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'A41', 'component': 'IKr', 'units': 'per_millisecond'},
    {'name': 'B41', 'component': 'IKr', 'units': 'per_millivolt'},
    {'name': 'q41', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'A51', 'component': 'IKr', 'units': 'per_millisecond'},
    {'name': 'B51', 'component': 'IKr', 'units': 'per_millivolt'},
    {'name': 'q51', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'A52', 'component': 'IKr', 'units': 'per_millisecond'},
    {'name': 'B52', 'component': 'IKr', 'units': 'per_millivolt'},
    {'name': 'q52', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'A53', 'component': 'IKr', 'units': 'per_millisecond'},
    {'name': 'B53', 'component': 'IKr', 'units': 'per_millivolt'},
    {'name': 'q53', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'A61', 'component': 'IKr', 'units': 'per_millisecond'},
    {'name': 'B61', 'component': 'IKr', 'units': 'per_millivolt'},
    {'name': 'q61', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'A62', 'component': 'IKr', 'units': 'per_millisecond'},
    {'name': 'B62', 'component': 'IKr', 'units': 'per_millivolt'},
    {'name': 'q62', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'A63', 'component': 'IKr', 'units': 'per_millisecond'},
    {'name': 'B63', 'component': 'IKr', 'units': 'per_millivolt'},
    {'name': 'q63', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'Kmax', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'Ku', 'component': 'IKr', 'units': 'per_millisecond'},
    {'name': 'n', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'halfmax', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'Kt', 'component': 'IKr', 'units': 'per_millisecond'},
    {'name': 'Vhalf', 'component': 'IKr', 'units': 'millivolt'},
    {'name': 'Temp', 'component': 'IKr', 'units': 'dimensionless'},
    {'name': 'GKs_b', 'component': 'IKs', 'units': 'milliS_per_microF'},
    {'name': 'txs1_max', 'component': 'IKs', 'units': 'millisecond'},
    {'name': 'GK1_b', 'component': 'IK1', 'units': 'milliS_per_microF'},
    {'name': 'kna1', 'component': 'INaCa_i', 'units': 'per_millisecond'},
    {'name': 'kna2', 'component': 'INaCa_i', 'units': 'per_millisecond'},
    {'name': 'kna3', 'component': 'INaCa_i', 'units': 'per_millisecond'},
    {'name': 'kasymm', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'wna', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'wca', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'wnaca', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'kcaon', 'component': 'INaCa_i', 'units': 'per_millisecond'},
    {'name': 'kcaoff', 'component': 'INaCa_i', 'units': 'per_millisecond'},
    {'name': 'qna', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'qca', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'KmCaAct', 'component': 'INaCa_i', 'units': 'millimolar'},
    {'name': 'Gncx_b', 'component': 'INaCa_i', 'units': 'milliS_per_microF'},
    {'name': 'k1p', 'component': 'INaK', 'units': 'per_millisecond'},
    {'name': 'k1m', 'component': 'INaK', 'units': 'per_millisecond'},
    {'name': 'k2p', 'component': 'INaK', 'units': 'per_millisecond'},
    {'name': 'k2m', 'component': 'INaK', 'units': 'per_millisecond'},
    {'name': 'k3p', 'component': 'INaK', 'units': 'per_millisecond'},
    {'name': 'k3m', 'component': 'INaK', 'units': 'per_millisecond'},
    {'name': 'k4p', 'component': 'INaK', 'units': 'per_millisecond'},
    {'name': 'k4m', 'component': 'INaK', 'units': 'per_millisecond'},
    {'name': 'Knai0', 'component': 'INaK', 'units': 'millimolar'},
    {'name': 'Knao0', 'component': 'INaK', 'units': 'millimolar'},
    {'name': 'delta', 'component': 'INaK', 'units': 'millivolt'},
    {'name': 'Kki', 'component': 'INaK', 'units': 'per_millisecond'},
    {'name': 'Kko', 'component': 'INaK', 'units': 'per_millisecond'},
    {'name': 'MgADP', 'component': 'INaK', 'units': 'millimolar'},
    {'name': 'MgATP', 'component': 'INaK', 'units': 'millimolar'},
    {'name': 'Kmgatp', 'component': 'INaK', 'units': 'millimolar'},
    {'name': 'H', 'component': 'INaK', 'units': 'millimolar'},
    {'name': 'eP', 'component': 'INaK', 'units': 'dimensionless'},
    {'name': 'Khp', 'component': 'INaK', 'units': 'millimolar'},
    {'name': 'Knap', 'component': 'INaK', 'units': 'millimolar'},
    {'name': 'Kxkur', 'component': 'INaK', 'units': 'millimolar'},
    {'name': 'Pnak_b', 'component': 'INaK', 'units': 'milliS_per_microF'},
    {'name': 'GKb_b', 'component': 'IKb', 'units': 'milliS_per_microF'},
    {'name': 'PNab', 'component': 'INab', 'units': 'milliS_per_microF'},
    {'name': 'PCab', 'component': 'ICab', 'units': 'milliS_per_microF'},
    {'name': 'GpCa', 'component': 'IpCa', 'units': 'milliS_per_microF'},
    {'name': 'KmCap', 'component': 'IpCa', 'units': 'millimolar'},
    {'name': 'bt', 'component': 'ryr', 'units': 'millisecond'},
    {'name': 'Jrel_scaling_factor', 'component': 'ryr', 'units': 'dimensionless'},
    {'name': 'Jup_b', 'component': 'SERCA', 'units': 'dimensionless'},
    {'name': 'frt', 'component': 'membrane', 'units': 'per_millivolt'},
    {'name': 'cmdnmax', 'component': 'intracellular_ions', 'units': 'millimolar'},
    {'name': 'Ahs', 'component': 'INa', 'units': 'dimensionless'},
    {'name': 'thLp', 'component': 'INaL', 'units': 'millisecond'},
    {'name': 'GNaL', 'component': 'INaL', 'units': 'milliS_per_microF'},
    {'name': 'Gto', 'component': 'Ito', 'units': 'milliS_per_microF'},
    {'name': 'Aff', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'PCa', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'tjca', 'component': 'ICaL', 'units': 'millisecond'},
    {'name': 'v0', 'component': 'ICaL', 'units': 'millivolt'},
    {'name': 'GKr', 'component': 'IKr', 'units': 'milliS_per_microF'},
    {'name': 'GKs', 'component': 'IKs', 'units': 'milliS_per_microF'},
    {'name': 'GK1', 'component': 'IK1', 'units': 'milliS_per_microF'},
    {'name': 'vcell', 'component': 'cell_geometry', 'units': 'microliter'},
    {'name': 'GKb', 'component': 'IKb', 'units': 'milliS_per_microF'},
    {'name': 'v0', 'component': 'INab', 'units': 'millivolt'},
    {'name': 'v0', 'component': 'ICab', 'units': 'millivolt'},
    {'name': 'a_rel', 'component': 'ryr', 'units': 'millisecond'},
    {'name': 'btp', 'component': 'ryr', 'units': 'millisecond'},
    {'name': 'upScale', 'component': 'SERCA', 'units': 'dimensionless'},
    {'name': 'ffrt', 'component': 'membrane', 'units': 'coulomb_per_mole_millivolt'},
    {'name': 'Afs', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'PCap', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'PCaNa', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'PCaK', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'B_1', 'component': 'ICaL', 'units': 'per_millivolt'},
    {'name': 'B_2', 'component': 'ICaL', 'units': 'per_millivolt'},
    {'name': 'B_3', 'component': 'ICaL', 'units': 'per_millivolt'},
    {'name': 'Ageo', 'component': 'cell_geometry', 'units': 'centimeter_squared'},
    {'name': 'B', 'component': 'INab', 'units': 'per_millivolt'},
    {'name': 'B', 'component': 'ICab', 'units': 'per_millivolt'},
    {'name': 'a_relp', 'component': 'ryr', 'units': 'millisecond'},
    {'name': 'PCaNap', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'PCaKp', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'Acap', 'component': 'cell_geometry', 'units': 'centimeter_squared'},
    {'name': 'vmyo', 'component': 'cell_geometry', 'units': 'microliter'},
    {'name': 'vnsr', 'component': 'cell_geometry', 'units': 'microliter'},
    {'name': 'vjsr', 'component': 'cell_geometry', 'units': 'microliter'},
    {'name': 'vss', 'component': 'cell_geometry', 'units': 'microliter'},
    {'name': 'h10_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h11_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h12_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k1_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k2_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k5_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'Gncx', 'component': 'INaCa_i', 'units': 'milliS_per_microF'},
    {'name': 'h10_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h11_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h12_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k1_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k2_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k5_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'b1', 'component': 'INaK', 'units': 'dimensionless'},
    {'name': 'a2', 'component': 'INaK', 'units': 'dimensionless'},
    {'name': 'a4', 'component': 'INaK', 'units': 'dimensionless'},
    {'name': 'Pnak', 'component': 'INaK', 'units': 'milliS_per_microF'},
    {'name': None, 'component': None, 'units': None}])
    
legend_algebraic = pd.DataFrame([
    {'name': 'Istim', 'component': 'membrane', 'units': 'microA_per_microF'},
    {'name': 'mss', 'component': 'INa', 'units': 'dimensionless'},
    {'name': 'hss', 'component': 'INa', 'units': 'dimensionless'},
    {'name': 'hLss', 'component': 'INaL', 'units': 'dimensionless'},
    {'name': 'hLssp', 'component': 'INaL', 'units': 'dimensionless'},
    {'name': 'ass', 'component': 'Ito', 'units': 'dimensionless'},
    {'name': 'iss', 'component': 'Ito', 'units': 'dimensionless'},
    {'name': 'dss', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'fss', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'km2n', 'component': 'ICaL', 'units': 'per_millisecond'},
    {'name': 'xs1ss', 'component': 'IKs', 'units': 'dimensionless'},
    {'name': 'xk1ss', 'component': 'IK1', 'units': 'dimensionless'},
    {'name': 'vfrt', 'component': 'membrane', 'units': 'dimensionless'},
    {'name': 'tm', 'component': 'INa', 'units': 'millisecond'},
    {'name': 'thf', 'component': 'INa', 'units': 'millisecond'},
    {'name': 'ths', 'component': 'INa', 'units': 'millisecond'},
    {'name': 'jss', 'component': 'INa', 'units': 'dimensionless'},
    {'name': 'ta', 'component': 'Ito', 'units': 'millisecond'},
    {'name': 'delta_epi', 'component': 'Ito', 'units': 'dimensionless'},
    {'name': 'fcass', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'anca', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'td', 'component': 'ICaL', 'units': 'millisecond'},
    {'name': 'tff', 'component': 'ICaL', 'units': 'millisecond'},
    {'name': 'tfs', 'component': 'ICaL', 'units': 'millisecond'},
    {'name': 'xs2ss', 'component': 'IKs', 'units': 'dimensionless'},
    {'name': 'txs1', 'component': 'IKs', 'units': 'millisecond'},
    {'name': 'txk1', 'component': 'IK1', 'units': 'millisecond'},
    {'name': 'tj', 'component': 'INa', 'units': 'millisecond'},
    {'name': 'hssp', 'component': 'INa', 'units': 'dimensionless'},
    {'name': 'mLss', 'component': 'INaL', 'units': 'dimensionless'},
    {'name': 'tiF_b', 'component': 'Ito', 'units': 'millisecond'},
    {'name': 'assp', 'component': 'Ito', 'units': 'dimensionless'},
    {'name': 'tfcaf', 'component': 'ICaL', 'units': 'millisecond'},
    {'name': 'tfcas', 'component': 'ICaL', 'units': 'millisecond'},
    {'name': 'tffp', 'component': 'ICaL', 'units': 'millisecond'},
    {'name': 'txs2', 'component': 'IKs', 'units': 'millisecond'},
    {'name': 'CaMKb', 'component': 'CaMK', 'units': 'millimolar'},
    {'name': 'thsp', 'component': 'INa', 'units': 'millisecond'},
    {'name': 'tjp', 'component': 'INa', 'units': 'millisecond'},
    {'name': 'tmL', 'component': 'INaL', 'units': 'millisecond'},
    {'name': 'tiS_b', 'component': 'Ito', 'units': 'millisecond'},
    {'name': 'tfcafp', 'component': 'ICaL', 'units': 'millisecond'},
    {'name': 'CaMKa', 'component': 'CaMK', 'units': 'millimolar'},
    {'name': 'tiF', 'component': 'Ito', 'units': 'millisecond'},
    {'name': 'Bcai', 'component': 'intracellular_ions', 'units': 'dimensionless'},
    {'name': 'tiS', 'component': 'Ito', 'units': 'millisecond'},
    {'name': 'Bcass', 'component': 'intracellular_ions', 'units': 'dimensionless'},
    {'name': 'dti_develop', 'component': 'Ito', 'units': 'dimensionless'},
    {'name': 'Bcajsr', 'component': 'intracellular_ions', 'units': 'dimensionless'},
    {'name': 'dti_recover', 'component': 'Ito', 'units': 'dimensionless'},
    {'name': 'ENa', 'component': 'reversal_potentials', 'units': 'millivolt'},
    {'name': 'tiFp', 'component': 'Ito', 'units': 'millisecond'},
    {'name': 'tiSp', 'component': 'Ito', 'units': 'millisecond'},
    {'name': 'EK', 'component': 'reversal_potentials', 'units': 'millivolt'},
    {'name': 'EKs', 'component': 'reversal_potentials', 'units': 'millivolt'},
    {'name': 'h', 'component': 'INa', 'units': 'dimensionless'},
    {'name': 'hp', 'component': 'INa', 'units': 'dimensionless'},
    {'name': 'fINap', 'component': 'INa', 'units': 'dimensionless'},
    {'name': 'INa', 'component': 'INa', 'units': 'microA_per_microF'},
    {'name': 'fINaLp', 'component': 'INaL', 'units': 'dimensionless'},
    {'name': 'INaL', 'component': 'INaL', 'units': 'microA_per_microF'},
    {'name': 'AiF', 'component': 'Ito', 'units': 'dimensionless'},
    {'name': 'AiS', 'component': 'Ito', 'units': 'dimensionless'},
    {'name': 'i', 'component': 'Ito', 'units': 'dimensionless'},
    {'name': 'ip', 'component': 'Ito', 'units': 'dimensionless'},
    {'name': 'fItop', 'component': 'Ito', 'units': 'dimensionless'},
    {'name': 'Ito', 'component': 'Ito', 'units': 'microA_per_microF'},
    {'name': 'f', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'Afcaf', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'Afcas', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'fca', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'fp', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'fcap', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'A_1', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'U_1', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'PhiCaL', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'A_2', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'U_2', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'PhiCaNa', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'A_3', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'U_3', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'PhiCaK', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'fICaLp', 'component': 'ICaL', 'units': 'dimensionless'},
    {'name': 'ICaL', 'component': 'ICaL', 'units': 'microA_per_microF'},
    {'name': 'ICaNa', 'component': 'ICaL', 'units': 'microA_per_microF'},
    {'name': 'Jrel_inf_temp', 'component': 'ryr', 'units': 'dimensionless'},
    {'name': 'Jrel_temp', 'component': 'ryr', 'units': 'dimensionless'},
    {'name': 'ICaK', 'component': 'ICaL', 'units': 'microA_per_microF'},
    {'name': 'Jrel_inf', 'component': 'ryr', 'units': 'dimensionless'},
    {'name': 'Jrel_infp', 'component': 'ryr', 'units': 'dimensionless'},
    {'name': 'IKr', 'component': 'IKr', 'units': 'microA_per_microF'},
    {'name': 'tau_rel_temp', 'component': 'ryr', 'units': 'millisecond'},
    {'name': 'tau_relp_temp', 'component': 'ryr', 'units': 'millisecond'},
    {'name': 'KsCa', 'component': 'IKs', 'units': 'dimensionless'},
    {'name': 'tau_rel', 'component': 'ryr', 'units': 'millisecond'},
    {'name': 'tau_relp', 'component': 'ryr', 'units': 'millisecond'},
    {'name': 'IKs', 'component': 'IKs', 'units': 'microA_per_microF'},
    {'name': 'rk1', 'component': 'IK1', 'units': 'millisecond'},
    {'name': 'IK1', 'component': 'IK1', 'units': 'microA_per_microF'},
    {'name': 'hca', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'hna', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h1_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h2_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h3_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h4_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h5_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h6_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h7_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h8_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h9_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k3p_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k3pp_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k3_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k4p_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k4pp_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k4_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k6_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k7_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k8_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'x1_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'x2_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'x3_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'x4_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'E1_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'E2_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'E3_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'E4_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'allo_i', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'JncxNa_i', 'component': 'INaCa_i', 'units': 'millimolar_per_millisecond'},
    {'name': 'JncxCa_i', 'component': 'INaCa_i', 'units': 'millimolar_per_millisecond'},
    {'name': 'INaCa_i', 'component': 'INaCa_i', 'units': 'microA_per_microF'},
    {'name': 'h1_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h2_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h3_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h4_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h5_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h6_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h7_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h8_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'h9_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k3p_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k3pp_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k3_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k4p_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k4pp_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k4_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k6_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k7_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'k8_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'x1_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'x2_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'x3_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'x4_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'E1_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'E2_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'E3_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'E4_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'allo_ss', 'component': 'INaCa_i', 'units': 'dimensionless'},
    {'name': 'JncxNa_ss', 'component': 'INaCa_i', 'units': 'millimolar_per_millisecond'},
    {'name': 'JncxCa_ss', 'component': 'INaCa_i', 'units': 'millimolar_per_millisecond'},
    {'name': 'INaCa_ss', 'component': 'INaCa_i', 'units': 'microA_per_microF'},
    {'name': 'Knai', 'component': 'INaK', 'units': 'millimolar'},
    {'name': 'Knao', 'component': 'INaK', 'units': 'millimolar'},
    {'name': 'P', 'component': 'INaK', 'units': 'dimensionless'},
    {'name': 'a1', 'component': 'INaK', 'units': 'dimensionless'},
    {'name': 'b2', 'component': 'INaK', 'units': 'dimensionless'},
    {'name': 'a3', 'component': 'INaK', 'units': 'dimensionless'},
    {'name': 'b3', 'component': 'INaK', 'units': 'dimensionless'},
    {'name': 'b4', 'component': 'INaK', 'units': 'dimensionless'},
    {'name': 'x1', 'component': 'INaK', 'units': 'dimensionless'},
    {'name': 'x2', 'component': 'INaK', 'units': 'dimensionless'},
    {'name': 'x3', 'component': 'INaK', 'units': 'dimensionless'},
    {'name': 'x4', 'component': 'INaK', 'units': 'dimensionless'},
    {'name': 'E1', 'component': 'INaK', 'units': 'dimensionless'},
    {'name': 'E2', 'component': 'INaK', 'units': 'dimensionless'},
    {'name': 'E3', 'component': 'INaK', 'units': 'dimensionless'},
    {'name': 'E4', 'component': 'INaK', 'units': 'dimensionless'},
    {'name': 'JnakNa', 'component': 'INaK', 'units': 'millimolar_per_millisecond'},
    {'name': 'JnakK', 'component': 'INaK', 'units': 'millimolar_per_millisecond'},
    {'name': 'INaK', 'component': 'INaK', 'units': 'microA_per_microF'},
    {'name': 'xkb', 'component': 'IKb', 'units': 'dimensionless'},
    {'name': 'IKb', 'component': 'IKb', 'units': 'microA_per_microF'},
    {'name': 'A', 'component': 'INab', 'units': 'microA_per_microF'},
    {'name': 'JdiffK', 'component': 'diff', 'units': 'millimolar_per_millisecond'},
    {'name': 'U', 'component': 'INab', 'units': 'dimensionless'},
    {'name': 'INab', 'component': 'INab', 'units': 'microA_per_microF'},
    {'name': 'A', 'component': 'ICab', 'units': 'microA_per_microF'},
    {'name': 'JdiffNa', 'component': 'diff', 'units': 'millimolar_per_millisecond'},
    {'name': 'U', 'component': 'ICab', 'units': 'dimensionless'},
    {'name': 'ICab', 'component': 'ICab', 'units': 'microA_per_microF'},
    {'name': 'IpCa', 'component': 'IpCa', 'units': 'microA_per_microF'},
    {'name': 'Jdiff', 'component': 'diff', 'units': 'millimolar_per_millisecond'},
    {'name': 'fJrelp', 'component': 'ryr', 'units': 'dimensionless'},
    {'name': 'Jrel', 'component': 'ryr', 'units': 'millimolar_per_millisecond'},
    {'name': 'Jupnp', 'component': 'SERCA', 'units': 'millimolar_per_millisecond'},
    {'name': 'Jupp', 'component': 'SERCA', 'units': 'millimolar_per_millisecond'},
    {'name': 'fJupp', 'component': 'SERCA', 'units': 'dimensionless'},
    {'name': 'Jleak', 'component': 'SERCA', 'units': 'millimolar_per_millisecond'},
    {'name': 'Jup', 'component': 'SERCA', 'units': 'millimolar_per_millisecond'},
    {'name': 'Jtr', 'component': 'trans_flux', 'units': 'millimolar_per_millisecond'}])


assert(size_states == legend_states.shape[0])
assert(size_constants == legend_constants.shape[0])
assert(size_algebraic == legend_algebraic.shape[0])

legend = dict(states = legend_states,
              constants = legend_constants,
              algebraic = legend_algebraic)
              

def init_states_constants():
    
    S = np.zeros(size_states)
    C = np.zeros(size_constants)
    
    C[0] = 0
    C[1] = 140
    C[2] = 1.8
    C[3] = 5.4
    C[4] = 8314
    C[5] = 310
    C[6] = 96485
    C[7] = 1
    C[8] = 2
    C[9] = 1
    C[10] = 0.01
    C[11] = 0.0011
    S[0] = -88.00190465
    C[12] = 10
    C[13] = 100000000000000000
    C[14] = -80
    C[15] = 1000
    C[16] = 0.5
    C[17] = 0.15
    C[18] = 0.05
    C[19] = 0.00068
    C[20] = 0.05
    C[21] = 0.0015
    S[1] = 0.0125840447
    S[2] = 8.49e-05
    C[22] = 0.05
    C[23] = 0.00238
    C[24] = 0.07
    C[25] = 0.0005
    C[26] = 0.047
    C[27] = 0.00087
    C[28] = 1.124
    C[29] = 0.0087
    C[30] = 10
    C[31] = 0.8
    S[3] = 7.268004498
    S[4] = 7.268089977
    S[5] = 144.6555918
    S[6] = 144.6555651
    S[7] = 1.619574538
    S[8] = 1.571234014
    S[9] = 8.6e-05
    C[32] = 1
    C[33] = 0.01833
    C[34] = 39.57
    C[35] = 9.871
    C[36] = 11.64
    C[37] = 34.77
    C[38] = 6.765
    C[39] = 8.552
    C[40] = 77.42
    C[41] = 5.955
    S[10] = 0.007344121102
    C[42] = 82.9
    C[43] = 6.086
    C[44] = 0.99
    S[11] = 0.6981071913
    S[12] = 0.6980895801
    C[45] = 75
    C[46] = 0
    S[13] = 0.6979908432
    S[14] = 0.4549485525
    S[15] = 0.6979245865
    S[16] = 0.0001882617273
    C[47] = 200
    S[17] = 0.5008548855
    S[18] = 0.2693065357
    C[48] = 0.019957499999999975
    C[49] = 0.02
    S[19] = 0.001001097687
    S[20] = 0.9995541745
    S[21] = 0.5865061736
    S[22] = 0.0005100862934
    S[23] = 0.9995541823
    S[24] = 0.6393399482
    C[50] = 0.002
    C[51] = 1000
    C[52] = 0.0001007
    S[25] = 2.34e-9
    S[26] = 0.9999999909
    S[27] = 0.9102412777
    S[28] = 0.9999999909
    S[29] = 0.9998046777
    S[30] = 0.9999738312
    S[31] = 0.9999999909
    S[32] = 0.9999999909
    S[33] = 0.002749414044
    C[53] = 0.04658545454545456
    S[34] = 0.999637
    S[35] = 6.83208e-05
    S[36] = 1.80145e-08
    S[37] = 8.26619e-05
    S[38] = 0.00015551
    S[39] = 5.67623e-05
    S[40] = 0
    S[41] = 0
    S[42] = 0
    S[43] = 0
    C[54] = 0.0264
    C[55] = 4.631E-05
    C[56] = 4.843
    C[57] = 4.986E-06
    C[58] = -0.004226
    C[59] = 4.23
    C[60] = 0.001214
    C[61] = 0.008516
    C[62] = 4.962
    C[63] = 1.854E-05
    C[64] = -0.04641
    C[65] = 3.769
    C[66] = 0.0007868
    C[67] = 1.535E-08
    C[68] = 4.942
    C[69] = 5.455E-06
    C[70] = -0.1688
    C[71] = 4.156
    C[72] = 0.005509
    C[73] = 7.771E-09
    C[74] = 4.22
    C[75] = 0.001416
    C[76] = -0.02877
    C[77] = 1.459
    C[78] = 0.4492
    C[79] = 0.008595
    C[80] = 5
    C[81] = 0.3181
    C[82] = 3.613E-08
    C[83] = 4.663
    C[84] = 0.149
    C[85] = 0.004668
    C[86] = 2.412
    C[87] = 0.01241
    C[88] = 0.1725
    C[89] = 5.568
    C[90] = 0.3226
    C[91] = -0.0006575
    C[92] = 5
    C[93] = 0.008978
    C[94] = -0.02215
    C[95] = 5.682
    C[96] = 0
    C[97] = 0
    C[98] = 1
    C[99] = 1
    C[100] = 0
    C[101] = 1
    C[102] = 37
    C[103] = 0.006358000000000001
    C[104] = 817.3
    S[44] = 0.2707758025
    S[45] = 0.0001928503426
    C[105] = 0.3239783999999998
    S[46] = 0.9967597594
    C[106] = 15
    C[107] = 5
    C[108] = 88.12
    C[109] = 12.5
    C[110] = 6e4
    C[111] = 6e4
    C[112] = 5e3
    C[113] = 1.5e6
    C[114] = 5e3
    C[115] = 0.5224
    C[116] = 0.167
    C[117] = 150e-6
    C[118] = 0.0008
    C[119] = 949.5
    C[120] = 182.4
    C[121] = 687.2
    C[122] = 39.4
    C[123] = 1899
    C[124] = 79300
    C[125] = 639
    C[126] = 40
    C[127] = 9.073
    C[128] = 27.78
    C[129] = -0.155
    C[130] = 0.5
    C[131] = 0.3582
    C[132] = 0.05
    C[133] = 9.8
    C[134] = 1.698e-7
    C[135] = 1e-7
    C[136] = 4.2
    C[137] = 1.698e-7
    C[138] = 224
    C[139] = 292
    C[140] = 30
    C[141] = 0.003
    C[142] = 3.75e-10
    C[143] = 2.5e-8
    C[144] = 0.0005
    C[145] = 0.0005
    C[146] = 4.75
    S[47] = 2.5e-7
    S[48] = 3.12e-7
    C[147] = 1.0
    C[148] = 1.0
    C[149] = C[6]/(C[4]*C[5])
    C[150] = C[22]*1.30000 if C[0] == 1 else C[22]
    C[151] = 1.00000-C[44]
    C[152] = 3.00000*C[47]
    C[153] = C[48]*0.6 if C[0] == 1.0 else  C[48]
    C[154] = C[49]*4.00000 if C[0] == 1 else C[49]*4.00000 if C[0] == 2 else C[49]
    C[155] = 0.600000
    C[156] = C[52]*1.20000 if C[0] == 1 else C[52]*2.50000 if C[0] == 2 else C[52]
    C[157] = 75.0000
    C[158] = 0
    C[159] = C[53]*1.30000 if C[0] == 1 else C[53]*0.800000 if C[0] == 2 else C[53]
    C[160] = C[103]*1.40000 if C[0] == 1 else C[103]
    C[161] = C[105]*1.20000 if C[0] == 1 else C[105]*1.30000 if C[0] == 2 else C[105]
    C[162] = 1000.00*3.14000*C[11]*C[11]*C[10]
    C[163] = C[141]*0.600000 if C[0] == 1 else C[141]
    C[164] = 0.00000
    C[165] = 0.00000
    C[166] = 0.500000*C[146]
    C[167] = 1.25000*C[146]
    C[168] = 1.30000 if C[0] == 1 else 1.00000
    C[205] = 0.00000
    C[169] = C[6]*C[149]
    C[170] = 1.00000-C[155]
    C[171] = 1.10000*C[156]
    C[172] = 0.00125000*C[156]
    C[173] = 0.000357400*C[156]
    C[174] = 2.00000*C[149]
    C[175] = C[149]
    C[176] = C[149]
    C[177] = 2.00000*3.14000*C[11]*C[11]+2.00000*3.14000*C[11]*C[10]
    C[178] = C[149]
    C[179] = 2.00000*C[149]
    C[180] = 0.500000*C[167]
    C[181] = 0.00125000*C[171]
    C[182] = 0.000357400*C[171]
    C[183] = 2.00000*C[177]
    C[184] = 0.680000*C[162]
    C[185] = 0.0552000*C[162]
    C[186] = 0.00480000*C[162]
    C[187] = 0.0200000*C[162]
    C[188] = C[109]+1.00000+(C[1]/C[106])*(1.00000+C[1]/C[107])
    C[189] = (C[1]*C[1])/(C[188]*C[106]*C[107])
    C[190] = 1.00000/C[188]
    C[191] = C[190]*C[2]*C[113]
    C[192] = C[114]
    C[193] = C[114]
    C[194] = C[118]*1.10000 if C[0] == 1 else C[118]*1.40000 if C[0] == 2 else C[118]
    C[195] = C[109]+1.00000+(C[1]/C[106])*(1.00000+C[1]/C[107])
    C[196] = (C[1]*C[1])/(C[195]*C[106]*C[107])
    C[197] = 1.00000/C[195]
    C[198] = C[197]*C[2]*C[113]
    C[199] = C[114]
    C[200] = C[114]
    C[201] = C[120]*C[132]
    C[202] = C[121]
    C[203] = ((C[125]*C[133])/C[134])/(1.00000+C[133]/C[134])
    C[204] = C[140]*0.900000 if C[0] == 1 else C[140]*0.700000 if C[0] == 2 else C[140]
    return S, C


@njit
def compute_rates_algebraic(t, S, C, R=None, A=None):

    if R is None:
        R = np.zeros_like(S)
    if A is None:
        A = np.zeros(size_algebraic)

    R[43] = C[205]
    R[34] = (-(C[66]*np.exp(C[67]*S[0])*S[34]*np.exp(((C[102]-20.0000)*np.log(C[68]))/10.0000)-C[69]*np.exp(C[70]*S[0])*S[35]*np.exp(((C[102]-20.0000)*np.log(C[71]))/10.0000))+C[78]*np.exp(C[79]*S[0])*S[36]*np.exp(((C[102]-20.0000)*np.log(C[80]))/10.0000))-C[87]*np.exp(C[88]*S[0])*S[34]*np.exp(((C[102]-20.0000)*np.log(C[89]))/10.0000)
    R[35] = (((C[66]*np.exp(C[67]*S[0])*S[34]*np.exp(((C[102]-20.0000)*np.log(C[68]))/10.0000)-C[69]*np.exp(C[70]*S[0])*S[35]*np.exp(((C[102]-20.0000)*np.log(C[71]))/10.0000))-(C[60]*np.exp(C[61]*S[0])*S[35]*np.exp(((C[102]-20.0000)*np.log(C[62]))/10.0000)-C[63]*np.exp(C[64]*S[0])*S[39]*np.exp(((C[102]-20.0000)*np.log(C[65]))/10.0000)))+C[81]*np.exp(C[82]*S[0])*S[37]*np.exp(((C[102]-20.0000)*np.log(C[83]))/10.0000))-C[90]*np.exp(C[91]*S[0])*S[35]*np.exp(((C[102]-20.0000)*np.log(C[92]))/10.0000)
    R[36] = -(C[54]*np.exp(C[55]*S[0])*S[36]*np.exp(((C[102]-20.0000)*np.log(C[56]))/10.0000)-C[57]*np.exp(C[58]*S[0])*S[37]*np.exp(((C[102]-20.0000)*np.log(C[59]))/10.0000))-(C[78]*np.exp(C[79]*S[0])*S[36]*np.exp(((C[102]-20.0000)*np.log(C[80]))/10.0000)-C[87]*np.exp(C[88]*S[0])*S[34]*np.exp(((C[102]-20.0000)*np.log(C[89]))/10.0000))
    R[37] = ((C[54]*np.exp(C[55]*S[0])*S[36]*np.exp(((C[102]-20.0000)*np.log(C[56]))/10.0000)-C[57]*np.exp(C[58]*S[0])*S[37]*np.exp(((C[102]-20.0000)*np.log(C[59]))/10.0000))-(C[72]*np.exp(C[73]*S[0])*S[37]*np.exp(((C[102]-20.0000)*np.log(C[74]))/10.0000)-C[75]*np.exp(C[76]*S[0])*S[38]*np.exp(((C[102]-20.0000)*np.log(C[77]))/10.0000)))-(C[81]*np.exp(C[82]*S[0])*S[37]*np.exp(((C[102]-20.0000)*np.log(C[83]))/10.0000)-C[90]*np.exp(C[91]*S[0])*S[35]*np.exp(((C[102]-20.0000)*np.log(C[92]))/10.0000))
    R[38] = ((C[72]*np.exp(C[73]*S[0])*S[37]*np.exp(((C[102]-20.0000)*np.log(C[74]))/10.0000)-C[75]*np.exp(C[76]*S[0])*S[38]*np.exp(((C[102]-20.0000)*np.log(C[77]))/10.0000))-(C[84]*np.exp(C[85]*S[0])*S[38]*np.exp(((C[102]-20.0000)*np.log(C[86]))/10.0000)-C[93]*np.exp(C[94]*S[0])*S[39]*np.exp(((C[102]-20.0000)*np.log(C[95]))/10.0000)))-(((C[96]*C[97]*np.exp(C[98]*np.log(S[43])))/(np.exp(C[98]*np.log(S[43]))+C[99]))*S[38]-C[97]*S[41])
    R[39] = (((C[60]*np.exp(C[61]*S[0])*S[35]*np.exp(((C[102]-20.0000)*np.log(C[62]))/10.0000)-C[63]*np.exp(C[64]*S[0])*S[39]*np.exp(((C[102]-20.0000)*np.log(C[65]))/10.0000))+C[84]*np.exp(C[85]*S[0])*S[38]*np.exp(((C[102]-20.0000)*np.log(C[86]))/10.0000))-C[93]*np.exp(C[94]*S[0])*S[39]*np.exp(((C[102]-20.0000)*np.log(C[95]))/10.0000))-(((C[96]*C[97]*np.exp(C[98]*np.log(S[43])))/(np.exp(C[98]*np.log(S[43]))+C[99]))*S[39]-((C[97]*C[84]*np.exp(C[85]*S[0])*np.exp(((C[102]-20.0000)*np.log(C[86]))/10.0000))/(C[93]*np.exp(C[94]*S[0])*np.exp(((C[102]-20.0000)*np.log(C[95]))/10.0000)))*S[40])
    R[40] = ((((C[96]*C[97]*np.exp(C[98]*np.log(S[43])))/(np.exp(C[98]*np.log(S[43]))+C[99]))*S[39]-((C[97]*C[84]*np.exp(C[85]*S[0])*np.exp(((C[102]-20.0000)*np.log(C[86]))/10.0000))/(C[93]*np.exp(C[94]*S[0])*np.exp(((C[102]-20.0000)*np.log(C[95]))/10.0000)))*S[40])+(C[100]/(1.00000+np.exp(-(S[0]-C[101])/6.78900)))*S[42])-C[100]*S[40]
    R[41] = ((((C[96]*C[97]*np.exp(C[98]*np.log(S[43])))/(np.exp(C[98]*np.log(S[43]))+C[99]))*S[38]-C[97]*S[41])+(C[100]/(1.00000+np.exp(-(S[0]-C[101])/6.78900)))*S[42])-C[100]*S[41]
    R[42] = -((C[100]/(1.00000+np.exp(-(S[0]-C[101])/6.78900)))*S[42]-C[100]*S[41])-((C[100]/(1.00000+np.exp(-(S[0]-C[101])/6.78900)))*S[42]-C[100]*S[40])
    A[3] = 1.00000/(1.00000+np.exp((S[0]+87.6100)/7.48800))
    R[17] = (A[3]-S[17])/C[47]
    A[4] = 1.00000/(1.00000+np.exp((S[0]+93.8100)/7.48800))
    R[18] = (A[4]-S[18])/C[152]
    A[1] = 1.00000/(1.00000+np.exp(-(S[0]+C[34])/C[35]))
    A[13] = 1.00000/(C[38]*np.exp((S[0]+C[36])/C[37])+C[39]*np.exp(-(S[0]+C[40])/C[41]))
    R[10] = (A[1]-S[10])/A[13]
    A[2] = 1.00000/(1.00000+np.exp(((S[0]+C[42])-C[46])/C[43]))
    A[14] = 1.00000/(1.43200e-05*np.exp(-((S[0]+1.19600)-C[46])/6.28500)+6.14900*np.exp(((S[0]+0.509600)-C[46])/20.2700))
    R[11] = (A[2]-S[11])/A[14]
    A[15] = 1.00000/(0.00979400*np.exp(-((S[0]+17.9500)-C[46])/28.0500)+0.334300*np.exp(((S[0]+5.73000)-C[46])/56.6600))
    R[12] = (A[2]-S[12])/A[15]
    A[5] = 1.00000/(1.00000+np.exp(-(S[0]-14.3400)/14.8200))
    A[17] = 1.05150/(1.00000/(1.20890*(1.00000+np.exp(-(S[0]-18.4099)/29.3814)))+3.50000/(1.00000+np.exp((S[0]+100.000)/29.3814)))
    R[19] = (A[5]-S[19])/A[17]
    A[7] = 1.00000/(1.00000+np.exp(-(S[0]+3.94000)/4.23000))
    A[21] = 0.600000+1.00000/(np.exp(-0.0500000*(S[0]+6.00000))+np.exp(0.0900000*(S[0]+14.0000)))
    R[25] = (A[7]-S[25])/A[21]
    A[8] = 1.00000/(1.00000+np.exp((S[0]+19.5800)/3.69600))
    A[22] = 7.00000+1.00000/(0.00450000*np.exp(-(S[0]+20.0000)/10.0000)+0.00450000*np.exp((S[0]+20.0000)/10.0000))
    R[26] = (A[8]-S[26])/A[22]
    A[23] = 1000.00+1.00000/(3.50000e-05*np.exp(-(S[0]+5.00000)/4.00000)+3.50000e-05*np.exp((S[0]+5.00000)/6.00000))
    R[27] = (A[8]-S[27])/A[23]
    A[19] = A[8]
    R[30] = (A[19]-S[30])/C[157]
    A[9] = S[30]*1.00000
    A[20] = 1.00000/(C[51]/A[9]+np.power(1.00000+C[50]/S[2], 4.00000))
    R[33] = A[20]*C[51]-S[33]*A[9]
    A[10] = 1.00000/(1.00000+np.exp(-(S[0]+11.6000)/8.93200))
    A[25] = C[104]+1.00000/(0.000232600*np.exp((S[0]+48.2800)/17.8000)+0.00129200*np.exp(-(S[0]+210.000)/230.000))
    R[44] = (A[10]-S[44])/A[25]
    A[11] = 1.00000/(1.00000+np.exp(-(S[0]+2.55380*C[3]+144.590)/(1.56920*C[3]+3.81150)))
    A[26] = 122.200/(np.exp(-(S[0]+127.200)/20.3600)+np.exp((S[0]+236.800)/69.3300))
    R[46] = (A[11]-S[46])/A[26]
    A[36] = (C[20]*(1.00000-S[1]))/(1.00000+C[21]/S[2])
    R[1] = C[18]*A[36]*(A[36]+S[1])-C[19]*S[1]
    A[16] = A[2]
    A[27] = 2.03800+1.00000/(0.0213600*np.exp(-((S[0]+100.600)-C[46])/8.28100)+0.305200*np.exp(((S[0]+0.994100)-C[46])/38.4500))
    R[13] = (A[16]-S[13])/A[27]
    A[31] = 1.00000/(1.00000+np.exp(-(S[0]-24.3400)/14.8200))
    R[22] = (A[31]-S[22])/A[17]
    A[32] = 7.00000+1.00000/(0.0400000*np.exp(-(S[0]-4.00000)/7.00000)+0.0400000*np.exp((S[0]-4.00000)/7.00000))
    R[28] = (A[19]-S[28])/A[32]
    A[33] = 100.000+1.00000/(0.000120000*np.exp(-S[0]/3.00000)+0.000120000*np.exp(S[0]/7.00000))
    R[29] = (A[19]-S[29])/A[33]
    A[34] = 2.50000*A[22]
    R[31] = (A[8]-S[31])/A[34]
    A[24] = A[10]
    A[35] = 1.00000/(0.0100000*np.exp((S[0]-50.0000)/20.0000)+0.0193000*np.exp(-(S[0]+66.5400)/31.0000))
    R[45] = (A[24]-S[45])/A[35]
    A[28] = 1.00000/(1.00000+np.exp(((S[0]+89.1000)-C[46])/6.08600))
    A[37] = 3.00000*A[15]
    R[14] = (A[28]-S[14])/A[37]
    A[38] = 1.46000*A[27]
    R[15] = (A[16]-S[15])/A[38]
    A[29] = 1.00000/(1.00000+np.exp(-(S[0]+42.8500)/5.26400))
    A[39] = A[13]
    R[16] = (A[29]-S[16])/A[39]
    A[41] = 2.50000*A[32]
    R[32] = (A[19]-S[32])/A[41]
    A[6] = 1.00000/(1.00000+np.exp((S[0]+43.9400)/5.71100))
    A[18] = 1.00000-0.950000/(1.00000+np.exp((S[0]+70.0000)/5.00000)) if C[0] == 1 else 1
    A[30] = 4.56200+1.00000/(0.393300*np.exp(-(S[0]+100.000)/100.000)+0.0800400*np.exp((S[0]+50.0000)/16.5900))
    A[43] = A[30]*A[18]
    R[20] = (A[6]-S[20])/A[43]
    A[40] = 23.6200+1.00000/(0.00141600*np.exp(-(S[0]+96.5200)/59.0500)+1.78000e-08*np.exp((S[0]+114.100)/8.07900))
    A[45] = A[40]*A[18]
    R[21] = (A[6]-S[21])/A[45]
    A[47] = 1.35400+0.000100000/(np.exp((S[0]-167.400)/15.8900)+np.exp(-(S[0]-12.2300)/0.215400))
    A[49] = 1.00000-0.500000/(1.00000+np.exp((S[0]+70.0000)/20.0000))
    A[51] = A[47]*A[49]*A[43]
    R[23] = (A[6]-S[23])/A[51]
    A[52] = A[47]*A[49]*A[45]
    R[24] = (A[6]-S[24])/A[52]
    A[67] = C[155]*S[26]+C[170]*S[27]
    A[68] = 0.300000+0.600000/(1.00000+np.exp((S[0]-10.0000)/10.0000))
    A[69] = 1.00000-A[68]
    A[70] = A[68]*S[28]+A[69]*S[29]
    A[71] = C[155]*S[31]+C[170]*S[27]
    A[72] = A[68]*S[32]+A[69]*S[29]
    A[12] = S[0]*C[149]
    A[73] = (4.00000*C[169]*(S[2]*np.exp(2.00000*A[12])-0.341000*C[2]))/C[174]
    A[74] = C[174]*(S[0]-C[158])
    A[75] = A[73]*(1.00000-0.500000*A[74]) if ((-1e-7 <= A[74]) and (A[74] <= 1e-7)) else (A[73]*A[74])/(np.exp(A[74])-1.00000)
    A[42] = A[36]+S[1]
    A[82] = 1.00000/(1.00000+C[17]/A[42])
    A[83] = (1.00000-A[82])*C[156]*A[75]*S[25]*(A[67]*(1.00000-S[33])+S[30]*A[70]*S[33])+A[82]*C[171]*A[75]*S[25]*(A[71]*(1.00000-S[33])+S[30]*A[72]*S[33])
    A[85] = (C[166]*-A[83])/(1.00000+1.00000*(np.power(1.50000/S[8], 8.00000)))
    A[88] = A[85]*1.70000 if C[0] == 2.00000 else A[85]
    A[91] = C[146]/(1.00000+0.0123000/S[8])
    A[94] = 0.001 if A[91] < 0.001 else A[91]
    R[47] = (A[88]-S[47])/A[94]
    A[86] = (C[180]*-A[83])/(1.00000+np.power(1.50000/S[8], 8.00000))
    A[89] = A[86]*1.70000 if C[0] == 2 else A[86]
    A[92] = C[167]/(1.00000+0.0123000/S[8])
    A[95] = 0.001 if A[92] < 0.001 else A[92]
    R[48] = (A[89]-S[48])/A[95]
    A[53] = ((C[4]*C[5])/C[6])*np.log(C[3]/S[5])
    A[61] = 1.00000/(1.00000+np.exp((S[0]-213.600)/151.200))
    A[62] = 1.00000-A[61]
    A[63] = A[61]*S[20]+A[62]*S[21]
    A[64] = A[61]*S[23]+A[62]*S[24]
    A[65] = 1.00000/(1.00000+C[17]/A[42])
    A[66] = C[154]*(S[0]-A[53])*((1.00000-A[65])*S[19]*A[63]+A[65]*S[22]*A[64])
    A[90] = C[159]*(np.power(C[3]/5.40000, 1.0/2))*S[38]*(S[0]-A[53])
    A[54] = ((C[4]*C[5])/C[6])*np.log((C[3]+C[33]*C[1])/(S[5]+C[33]*S[3]))
    A[93] = 1.00000+0.600000/(1.00000+np.power(3.80000e-05/S[9], 1.40000))
    A[96] = C[160]*A[93]*S[44]*S[45]*(S[0]-A[54])
    A[97] = 1.00000/(1.00000+np.exp(((S[0]+105.800)-2.60000*C[3])/9.49300))
    A[98] = C[161]*(np.power(C[3], 1.0/2))*A[97]*S[46]*(S[0]-A[53])
    A[162] = C[128]*np.exp(((1.00000-C[129])*S[0]*C[6])/(3.00000*C[4]*C[5]))
    A[166] = (C[123]*(np.power(C[3]/C[131], 2.00000)))/((np.power(1.00000+C[1]/A[162], 3.00000)+np.power(1.00000+C[3]/C[131], 2.00000))-1.00000)
    A[163] = C[136]/(1.00000+C[135]/C[137]+S[3]/C[138]+S[5]/C[139])
    A[167] = (C[124]*A[163]*C[135])/(1.00000+C[133]/C[134])
    A[161] = C[127]*np.exp((C[129]*S[0]*C[6])/(3.00000*C[4]*C[5]))
    A[164] = (C[119]*(np.power(S[3]/A[161], 3.00000)))/((np.power(1.00000+S[3]/A[161], 3.00000)+np.power(1.00000+S[5]/C[130], 2.00000))-1.00000)
    A[165] = (C[122]*(np.power(C[1]/A[162], 3.00000)))/((np.power(1.00000+C[1]/A[162], 3.00000)+np.power(1.00000+C[3]/C[131], 2.00000))-1.00000)
    A[168] = (C[126]*(np.power(S[5]/C[130], 2.00000)))/((np.power(1.00000+S[3]/A[161], 3.00000)+np.power(1.00000+S[5]/C[130], 2.00000))-1.00000)
    A[169] = C[203]*A[164]*C[202]+A[165]*A[168]*A[167]+C[202]*A[168]*A[167]+A[167]*A[164]*C[202]
    A[170] = A[165]*C[201]*A[168]+A[164]*C[202]*A[166]+A[166]*C[201]*A[168]+C[202]*A[166]*A[168]
    A[171] = C[202]*A[166]*C[203]+A[167]*A[165]*C[201]+A[165]*C[201]*C[203]+A[166]*C[203]*C[201]
    A[172] = A[168]*A[167]*A[165]+A[166]*C[203]*A[164]+A[165]*C[203]*A[164]+A[167]*A[165]*A[164]
    A[173] = A[169]/(A[169]+A[170]+A[171]+A[172])
    A[174] = A[170]/(A[169]+A[170]+A[171]+A[172])
    A[177] = 3.00000*(A[173]*A[166]-A[174]*A[167])
    A[175] = A[171]/(A[169]+A[170]+A[171]+A[172])
    A[176] = A[172]/(A[169]+A[170]+A[171]+A[172])
    A[178] = 2.00000*(A[176]*C[201]-A[175]*A[164])
    A[179] = C[204]*(C[7]*A[177]+C[9]*A[178])
    A[180] = 1.00000/(1.00000+np.exp(-(S[0]-14.4800)/18.3400))
    A[181] = C[163]*A[180]*(S[0]-A[53])
    A[0] = C[14] if ((t >= C[12]) and (t <= C[13]) and (((t-C[12])-np.floor((t-C[12])/C[15])*C[15]) <= C[16])) else 0
    A[183] = (S[6]-S[5])/2.00000
    R[5] = (-((A[66]+A[90]+A[96]+A[98]+A[181]+A[0])-2.00000*A[179])*C[32]*C[183])/(C[6]*C[184])+(A[183]*C[187])/C[184]
    A[79] = (0.750000*C[169]*(S[6]*np.exp(A[12])-C[3]))/C[176]
    A[80] = C[176]*(S[0]-C[158])
    A[81] = A[79]*(1.00000-0.500000*A[80]) if ((-1e-7 <= A[80]) and (A[80] <= 1e-7)) else (A[79]*A[80])/(np.exp(A[80])-1.00000)
    A[87] = (1.00000-A[82])*C[173]*A[81]*S[25]*(A[67]*(1.00000-S[33])+S[30]*A[70]*S[33])+A[82]*C[182]*A[81]*S[25]*(A[71]*(1.00000-S[33])+S[30]*A[72]*S[33])
    R[6] = (-A[87]*C[32]*C[183])/(C[6]*C[187])-A[183]
    A[50] = ((C[4]*C[5])/C[6])*np.log(C[1]/S[3])
    A[55] = C[44]*S[11]+C[151]*S[12]
    A[56] = C[44]*S[11]+C[151]*S[14]
    A[57] = 1.00000/(1.00000+C[17]/A[42])
    A[58] = C[45]*(S[0]-A[50])*(np.power(S[10], 3.00000))*((1.00000-A[57])*A[55]*S[13]+A[57]*A[56]*S[15])
    A[59] = 1.00000/(1.00000+C[17]/A[42])
    A[60] = C[153]*(S[0]-A[50])*S[16]*((1.00000-A[59])*S[17]+A[59]*S[18])
    A[127] = 1.00000/(1.00000+np.power(C[117]/S[9], 2.00000))
    A[100] = np.exp((C[115]*S[0]*C[6])/(C[4]*C[5]))
    A[107] = 1.00000+(C[1]/C[108])*(1.00000+1.00000/A[100])
    A[108] = C[1]/(C[108]*A[100]*A[107])
    A[111] = A[108]*C[112]
    A[101] = 1.00000+(S[3]/C[108])*(1.00000+A[100])
    A[102] = (S[3]*A[100])/(C[108]*A[101])
    A[114] = A[102]*C[112]
    A[104] = 1.00000+(S[3]/C[106])*(1.00000+S[3]/C[107])
    A[105] = (S[3]*S[3])/(A[104]*C[106]*C[107])
    A[117] = A[105]*A[102]*C[110]
    A[118] = A[108]*C[189]*C[110]
    A[109] = 1.00000/A[107]
    A[110] = A[109]*C[111]
    A[112] = A[110]+A[111]
    A[99] = np.exp((C[116]*S[0]*C[6])/(C[4]*C[5]))
    A[103] = 1.00000/A[101]
    A[113] = (A[103]*C[111])/A[99]
    A[115] = A[113]+A[114]
    A[106] = 1.00000/A[104]
    A[116] = A[106]*S[9]*C[113]
    A[119] = C[192]*A[115]*(A[117]+A[116])+C[193]*A[117]*(C[192]+A[112])
    A[120] = C[191]*A[117]*(A[115]+C[193])+A[115]*A[116]*(C[191]+A[118])
    A[121] = C[191]*A[112]*(A[117]+A[116])+A[118]*A[116]*(C[192]+A[112])
    A[122] = C[192]*A[118]*(A[115]+C[193])+A[112]*C[193]*(C[191]+A[118])
    A[123] = A[119]/(A[119]+A[120]+A[121]+A[122])
    A[124] = A[120]/(A[119]+A[120]+A[121]+A[122])
    A[125] = A[121]/(A[119]+A[120]+A[121]+A[122])
    A[126] = A[122]/(A[119]+A[120]+A[121]+A[122])
    A[128] = (3.00000*(A[126]*A[117]-A[123]*A[118])+A[125]*A[114])-A[124]*A[111]
    A[129] = A[124]*C[192]-A[123]*C[191]
    A[130] = 0.800000*C[194]*A[127]*(C[7]*A[128]+C[8]*A[129])
    A[182] = (C[142]*C[169]*(S[3]*np.exp(A[12])-C[1]))/C[178]
    A[184] = C[178]*(S[0]-C[164])
    A[185] = A[182]*(1.00000-0.500000*A[184]) if ((-1e-7 <= A[184]) and (A[184] <= 1e-7)) else (A[182]*A[184])/(np.exp(A[184])-1.00000)
    A[187] = (S[4]-S[3])/2.00000
    R[3] = (-(A[58]+A[60]+3.00000*A[130]+3.00000*A[179]+A[185])*C[183]*C[32])/(C[6]*C[184])+(A[187]*C[187])/C[184]
    A[76] = (0.750000*C[169]*(S[4]*np.exp(A[12])-C[1]))/C[175]
    A[77] = C[175]*(S[0]-C[158])
    A[78] = A[76]*(1.00000-0.500000*A[77]) if ((-1.00000e-07 <= A[77]) and (A[77] <= 1.00000e-07)) else (A[76]*A[77])/(np.exp(A[77])-1.00000)
    A[84] = (1.00000-A[82])*C[172]*A[78]*S[25]*(A[67]*(1.00000-S[33])+S[30]*A[70]*S[33])+A[82]*C[181]*A[78]*S[25]*(A[71]*(1.00000-S[33])+S[30]*A[72]*S[33])
    A[157] = 1.00000/(1.00000+np.power(C[117]/S[2], 2.00000))
    A[137] = 1.00000+(C[1]/C[108])*(1.00000+1.00000/A[100])
    A[138] = C[1]/(C[108]*A[100]*A[137])
    A[141] = A[138]*C[112]
    A[131] = 1.00000+(S[4]/C[108])*(1.00000+A[100])
    A[132] = (S[4]*A[100])/(C[108]*A[131])
    A[144] = A[132]*C[112]
    A[134] = 1.00000+(S[4]/C[106])*(1.00000+S[4]/C[107])
    A[135] = (S[4]*S[4])/(A[134]*C[106]*C[107])
    A[147] = A[135]*A[132]*C[110]
    A[148] = A[138]*C[196]*C[110]
    A[139] = 1.00000/A[137]
    A[140] = A[139]*C[111]
    A[142] = A[140]+A[141]
    A[133] = 1.00000/A[131]
    A[143] = (A[133]*C[111])/A[99]
    A[145] = A[143]+A[144]
    A[136] = 1.00000/A[134]
    A[146] = A[136]*S[2]*C[113]
    A[149] = C[199]*A[145]*(A[147]+A[146])+C[200]*A[147]*(C[199]+A[142])
    A[150] = C[198]*A[147]*(A[145]+C[200])+A[145]*A[146]*(C[198]+A[148])
    A[151] = C[198]*A[142]*(A[147]+A[146])+A[148]*A[146]*(C[199]+A[142])
    A[152] = C[199]*A[148]*(A[145]+C[200])+A[142]*C[200]*(C[198]+A[148])
    A[153] = A[149]/(A[149]+A[150]+A[151]+A[152])
    A[154] = A[150]/(A[149]+A[150]+A[151]+A[152])
    A[155] = A[151]/(A[149]+A[150]+A[151]+A[152])
    A[156] = A[152]/(A[149]+A[150]+A[151]+A[152])
    A[158] = (3.00000*(A[156]*A[147]-A[153]*A[148])+A[155]*A[144])-A[154]*A[141]
    A[159] = A[154]*C[199]-A[153]*C[198]
    A[160] = 0.200000*C[194]*A[157]*(C[7]*A[158]+C[8]*A[159])
    R[4] = (-(A[84]+3.00000*A[160])*C[32]*C[183])/(C[6]*C[187])-A[187]
    A[190] = (C[144]*S[9])/(C[145]+S[9])
    A[186] = (C[143]*4.00000*C[169]*(S[9]*np.exp(2.00000*A[12])-0.341000*C[2]))/C[179]
    A[188] = C[179]*(S[0]-C[165])
    A[189] = A[186]*(1.00000-0.500000*A[188]) if ((-1.00000e-07 <= A[188]) and (A[188] <= 1.00000e-07)) else (A[186]*A[188])/(np.exp(A[188])-1.00000)
    R[0] = -(A[58]+A[60]+A[66]+A[83]+A[84]+A[87]+A[90]+A[96]+A[98]+A[130]+A[160]+A[179]+A[185]+A[181]+A[190]+A[189]+A[0])
    A[191] = (S[2]-S[9])/0.200000
    A[192] = 1.00000/(1.00000+C[17]/A[42])
    A[193] = C[147]*((1.00000-A[192])*S[47]+A[192]*S[48])
    A[46] = 1.00000/(1.00000+(C[26]*C[27])/(np.power(C[27]+S[2], 2.00000))+(C[28]*C[29])/(np.power(C[29]+S[2], 2.00000)))
    R[2] = A[46]*(((-(A[83]-2.00000*A[160])*C[32]*C[183])/(2.00000*C[6]*C[187])+(A[193]*C[186])/C[187])-A[191])
    A[194] = (C[168]*0.00437500*S[9])/(S[9]+0.000920000)
    A[195] = (C[168]*2.75000*0.00437500*S[9])/((S[9]+0.000920000)-0.000170000)
    A[196] = 1.00000/(1.00000+C[17]/A[42])
    A[197] = (0.00393750*S[7])/15.0000
    A[198] = C[148]*(((1.00000-A[196])*A[194]+A[196]*A[195])-A[197])
    A[44] = 1.00000/(1.00000+(C[150]*C[23])/(np.power(C[23]+S[9], 2.00000))+(C[24]*C[25])/(np.power(C[25]+S[9], 2.00000)))
    R[9] = A[44]*(((-((A[190]+A[189])-2.00000*A[130])*C[32]*C[183])/(2.00000*C[6]*C[184])-(A[198]*C[185])/C[184])+(A[191]*C[187])/C[184])
    A[199] = (S[7]-S[8])/100.000
    R[7] = A[198]-(A[199]*C[186])/C[185]
    A[48] = 1.00000/(1.00000+(C[30]*C[31])/(np.power(C[31]+S[8], 2.00000)))
    R[8] = A[48]*(A[199]-A[193])
    return R, A
    
# let's JIT
# compute_rates_algebraic(0, *init_states_constants())


def compute_rates(t, S, C, R=None, A=None):
    return compute_rates_algebraic(t, S, C, R, A)[0]
