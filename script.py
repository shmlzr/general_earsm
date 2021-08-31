#!/usr/bin/env python

"""
    A script to work with GCEARSM in OpenFOAM

"""

import numpy as np
from shutil import copy
from tempfile import mkstemp
from shutil import move
from os import remove, close
import pyfoam


def change_coef(case, theta):
    """
    Changes coefficient entries in case/constant/RASProperties
    :param case: OpenFOAM case name
    :param theta: numpy array of coefficients
    """
    copy(case + '/constant/RASProperties', case + '/constant/RASProperties_old')
    tmp = 10 ** 12
    tmp2 = 10 ** 12
    maxIter = -1
    cc = False
    j = 0
    file_path = case + '/constant/RASProperties'
    print(file_path)
    fh, abs_path = mkstemp()
    with open(abs_path, 'w') as new_file:
        with open(file_path) as file:
            for i, line in enumerate(file):
                if 'Theta' in line:
                    new_file.write('    Theta   (\n')
                    tmp = i + 1
                    tmp2 = tmp + np.size(theta)

                elif tmp <= i < tmp2:
                    new_file.write('                ' + str(theta[j]) + '\n')
                    j += 1

                else:
                    new_file.write(line)

    close(fh)
    remove(file_path)
    move(abs_path, file_path)


def coef_library():
    library = ['T1', 'T1*l1', 'T1*l2', 'T1*pow(l1,2)', 'T1*l1*l2', 'T1*pow(l2,2)', 'T1*pow(l1,3)', 'T1*pow(l1,2)*l2',
               'T1*l1*pow(l2,2)', 'T1*pow(l2,3)', 'T1*pow(l1,4)', 'T1*pow(l1,3)*l2', 'T1*pow(l1,2)*pow(l2,2)',
               'T1*l1*pow(l2,3)',
               'T1*pow(l2,4)', 'T1*pow(l1,5)', 'T1*pow(l1,4)*l2', 'T1*pow(l1,3)*pow(l2,2)', 'T1*pow(l1,2)*pow(l2,3)',
               'T1*l1*pow(l2,4)', 'T1*pow(l2,5)', 'T1*pow(l1,6)', 'T1*pow(l1,5)*l2', 'T1*pow(l1,4)*pow(l2,2)',
               'T1*pow(l1,3)*pow(l2,3)', 'T1*pow(l1,2)*pow(l2,4)', 'T1*l1*pow(l2,5)', 'T1*pow(l2,6)',
               'T2', 'T2*l1', 'T2*l2', 'T2*pow(l1,2)', 'T2*l1*l2', 'T2*pow(l2,2)', 'T2*pow(l1,3)', 'T2*pow(l1,2)*l2',
               'T2*l1*pow(l2,2)', 'T2*pow(l2,3)', 'T2*pow(l1,4)', 'T2*pow(l1,3)*l2', 'T2*pow(l1,2)*pow(l2,2)',
               'T2*l1*pow(l2,3)',
               'T2*pow(l2,4)', 'T2*pow(l1,5)', 'T2*pow(l1,4)*l2', 'T2*pow(l1,3)*pow(l2,2)', 'T2*pow(l1,2)*pow(l2,3)',
               'T2*l1*pow(l2,4)', 'T2*pow(l2,5)', 'T2*pow(l1,6)', 'T2*pow(l1,5)*l2', 'T2*pow(l1,4)*pow(l2,2)',
               'T2*pow(l1,3)*pow(l2,3)', 'T2*pow(l1,2)*pow(l2,4)', 'T2*l1*pow(l2,5)', 'T2*pow(l2,6)',
               'T3', 'T3*l1', 'T3*l2', 'T3*pow(l1,2)', 'T3*l1*l2', 'T3*pow(l2,2)', 'T3*pow(l1,3)', 'T3*pow(l1,2)*l2',
               'T3*l1*pow(l2,2)', 'T3*pow(l2,3)', 'T3*pow(l1,4)', 'T3*pow(l1,3)*l2', 'T3*pow(l1,2)*pow(l2,2)',
               'T3*l1*pow(l2,3)',
               'T3*pow(l2,4)', 'T3*pow(l1,5)', 'T3*pow(l1,4)*l2', 'T3*pow(l1,3)*pow(l2,2)', 'T3*pow(l1,2)*pow(l2,3)',
               'T3*l1*pow(l2,4)', 'T3*pow(l2,5)', 'T3*pow(l1,6)', 'T3*pow(l1,5)*l2', 'T3*pow(l1,4)*pow(l2,2)',
               'T3*pow(l1,3)*pow(l2,3)', 'T3*pow(l1,2)*pow(l2,4)', 'T3*l1*pow(l2,5)', 'T3*pow(l2,6)']
    return library


def print_model(inds, val):
    tmp = ''
    for i, ind in enumerate(inds):
        tmp += ' + ' + str(val[i]) + '*' + str(coef_library()[ind])
        # tmp.append(coef_library()[ind])
    print('Model:', tmp)


def create_theta(inds, val):
    """
        Creates a dense theta vector from sparse information. All elements of theta besides inds/val input are zero,
        i.e. inactive.
    :param inds: list of indices of active terms
    :param val: list of corresponding coefficient values
    """
    print_model(inds, val)
    theta = np.zeros(84)
    for i, ind in enumerate(inds):
        theta[ind] = val[i]
    return theta


def compute_cost(y_pred, y_true):
    mse = np.average((y_true - y_pred) ** 2.0)
    return mse

if __name__ == '__main__':
    # define theta
    # indices according to library and corresponding coefficient values
    theta_ind, theta_val = [0, 10, 22, 25], [0.1, 0.2, 0.3, 0.9]
    theta_dense = create_theta(theta_ind, theta_val)

    # create new_case
    pyfoam.create_case('case', 'new_case')

    # change coefficients in new_case
    change_coef('new_case', theta_dense)

    # run case with openfoam solver
    pyfoam.run_case('new_case', 'simpleFoam')

    # get velocity field
    U = pyfoam.getRANSVector('new_case', 18000, 'U')

    # define validation velocity field
    U_true = 0.*np.ones([3, 15600]) # this is a dummy

    # compute cost, e.g. mse
    mse = compute_cost(U, U_true)
    print('MSE =', mse)
