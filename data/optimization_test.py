from bayes_opt import BayesianOptimization
from bayes_opt import UtilityFunction
from bayes_opt import SequentialDomainReductionTransformer
import numpy as np
import os
import time
import pathlib
import logging

path = (str)(pathlib.Path().resolve()) + '/share_folder/'


P_sys_target = 80.0
P_dis_target = 50.0
SV_target = 30.0


def before_launch():
    global P_sys_target
    P_sys_target = input("Define P_sys_target : \n")
    P_sys_target = float(P_sys_target)
    global P_dis_target
    P_dis_target = input("Define P_dis_target : \n")
    P_dis_target = float(P_dis_target)
    global SV_target
    SV_target = input("Define SV_target : \n")
    SV_target = float(SV_target)
    global start_coeff
    start_coeff = 1
    # start_coeff = input("Define start_coeff : \n")
    # start_coeff = int(start_coeff)  

def pre_register_func(target):
    return - max([abs(P_sys_target - target[0]) / P_sys_target, abs(P_dis_target - target[1]) / P_dis_target, abs(SV_target - target[2]) / SV_target])

def add_next_point(point : dict):
    num_params = len(point)
    with open(path + "next_point" + ".csv", "w") as file:
        for i in range(num_params):
            key = 'p' + str(i)
            file.write(str(point.get(key)) + '\n')

# for restoring from backup
def add_target_and_point(target : list, point : dict, ind : int):
    num_params = len(point)
    with open(path + "params_" + str(ind) + ".csv", "w") as file:
        for i in range(num_params):
            key = 'p' + str(i)
            file.write(str(point.get(key)) + '\n')

    num_targets = len(target)
    with open(path + "targets_" + str(ind) + ".csv", "w") as file:
        for i in range(num_targets):
            file.write(str(target[i]) + '\n')

before_launch()

log_filename = 'bayes_' + f'{(int)(time.time())}' + '' + '.log'

logging.basicConfig(level = logging.INFO, filename = f'log/{log_filename}', filemode = 'w', format = '%(asctime)s - %(message)s', datefmt = '%d/%m/%y %H:%M:%S')
logging.info(f'Bayesian optimization will be performed for case {P_sys_target}, {P_dis_target}, {SV_target}')
logging.info(f'starting index = {start_coeff}')

acq_func = UtilityFunction(kind = "ei")

params_to_register = np.zeros(7)
targets_to_register = np.zeros(3)
np.set_printoptions(suppress = True)

best_eps = -42.0
best_target = np.zeros(3)

for i in range(start_coeff):
    with open(path + "params_" + str(i) + ".csv", "r") as file:
        param_lines = file.readlines()
        for j in range(len(param_lines)):
            params_to_register[j] = np.float64(param_lines[j])

    with open(path + "targets_" + str(i) + ".csv", "r") as file:
        target_lines = file.readlines()
        for j in range(len(target_lines)):
            targets_to_register[j] = np.float64(target_lines[j])
    eps_to_push = pre_register_func(targets_to_register)

    if i == 0:
        parameter_bounds = {
            'p0' : (params_to_register[0] / 2.50, params_to_register[0] * 2.50),
            'p1' : (params_to_register[1] / 2.50, params_to_register[1] * 2.50),
            'p2' : (8.0, 25.0),
            'p3' : (params_to_register[3] / 1.50, params_to_register[3] * 1.50),
            'p4' : (params_to_register[4] / 1.50, params_to_register[4] * 1.50),
            'p5' : (params_to_register[5] / 2.00, params_to_register[5] * 2.00),
            'p6' : (params_to_register[6] / 2.00, params_to_register[6] * 2.00),
        }
        logging.info(f'Parameter bounds deParameter bounds definedfined {parameter_bounds}')

        window = [
            0.001, # 0.001
            0.001, # 0.001
            5.0,
            500.0, # 5%
            2.5, # 5%
            0.00005,
            500.0
        ]
        bounds_transformer = SequentialDomainReductionTransformer(
            eta = 0.95,
            gamma_osc = 0.8,
            minimum_window = window
        )
        optimizer = BayesianOptimization(
            f = None,
            pbounds = parameter_bounds,
            random_state = 42,
            allow_duplicate_points = True,
            bounds_transformer = bounds_transformer
        )
        optimizer.set_gp_params(alpha = 1e-3, n_restarts_optimizer = 25)

    if eps_to_push > -0.04:
        print("fucking finita")
        print(params_to_register, targets_to_register)
    if eps_to_push > best_eps:
        best_target = targets_to_register.copy()
        best_eps = eps_to_push
        logging.info(f'New best error : {-best_eps}')
        logging.info(f'New best target : {best_target}')
    optimizer.register(params = params_to_register, target = eps_to_push)
    logging.info(f'Parameters registered : {params_to_register}')
    logging.info(f'Target registered : {-eps_to_push}')
    # optimizer.set_bounds(optimizer._bounds_transformer.transform(optimizer._space))
    # logging.info(f'Parameter bounds updated')
    # logging.info(optimizer._space.bounds)
    logging.info("----------------------------------------------------------------------------")

stamp = os.stat(path + "next_target.csv").st_mtime
i = 0

while (eps_to_push < - 0.04):

    next_point = optimizer.suggest(acq_func)
    add_next_point(next_point)
    logging.info(f'New point added {next_point}')

    while (os.stat(path + "next_target.csv").st_mtime < stamp + 1.0):
        time.sleep(10.0)
    stamp = os.stat(path + "next_target.csv").st_mtime

    with open(path + "next_target.csv", "r") as file:
        target_lines = file.readlines()
    
    if len(target_lines) != 1:
        for j in range(len(target_lines)):
            targets_to_register[j] = np.float64(target_lines[j])

        add_target_and_point(targets_to_register, next_point, i + start_coeff)
        logging.info(f'Target {targets_to_register}; point {next_point} written to csv file')
        i += 1

        eps_to_push = pre_register_func(targets_to_register)
        if eps_to_push > best_eps:
            best_target = targets_to_register.copy()
            best_eps = eps_to_push
            logging.info(f'New best error : {-best_eps}')
            logging.info(f'New best target : {best_target}')

        optimizer.register(params = next_point, target = eps_to_push)
        logging.info(f'Parameters registered : {next_point}')
        logging.info(f'Target registered : {-eps_to_push}')
        if i > 10:
            optimizer.set_bounds(optimizer._bounds_transformer.transform(optimizer._space))
            logging.info(f'Parameter bounds updated')
            logging.info(optimizer._space.bounds)
        logging.info("----------------------------------------------------------------------------")