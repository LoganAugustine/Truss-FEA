# AERO 306 Truss FEA Project
# Open JDW pdf page 197
from Plotting import *

import numpy as np
from sympy import *


# Read .txt file
def read(file):
    f = open(file, "r")
    f_list = f.readlines()
    f.close()
    for i in range(len(f_list)):
        x = f_list[i].split(',')
        f_list[i] = x
    return f_list


# Node information
def node(new):
    for i in range(len(new)):
        if 'NumNodes' in new[i]:  # Identify position in list
            n = int(new[i][1])
    nodeinfo = []
    for j in range(n):
        row = []
        for k in range(len(new[j + 1])):
            row.append(int(new[j + 1][k]))
        nodeinfo.append(row)
    return n, nodeinfo


# Element information
def element(new):
    for i in range(len(new)):
        if 'NumElements' in new[i]:
            n = int(new[i][1])
    c = []
    el_con = []
    for i in range(len(new)):
        if 'Element' in new[i]:
            for j in range(n):
                for k in range(0, 3):
                    c.append(int(new[i + (j + 1)][k]))
                el_con.append(c)
                c = []
    lengths = []
    for i in range(len(new)):
        if 'ElemLengths' in new[i]:
            for j in range(n):
                lengths.append(float(new[i][j + 1]))
    thetas = []
    for i in range(len(new)):
        if 'ElemThetas' in new[i]:
            for j in range(n):
                thetas.append(float(new[i][j + 1]))
    return n, el_con, lengths, thetas


# Nodal DOFs
def nodalDOF(new, numElems):
    c = []
    dofs = []
    for i in range(len(new)):
        if 'Node DOFs\n' in new[i]:
            for j in range(0, numElems):
                for k in range(0, 4):
                    c.append(int(new[i + (j + 1)][k]))
                dofs.append(c)
                c = []
    return dofs


# Calculate K global
def KGlobal(E, A, numNodes, numElems, dofs, elements, thetas):
    # Create empty K_global
    K_global = np.zeros((numNodes * 2, numNodes * 2))
    dofMatrix = np.array(dofs)

    KGelem_list = []
    for i in range(numElems):
        c = cos(thetas[i])
        s = sin(thetas[i])
        KGelem = E * A / elements[i] * np.array([[c ** 2, c * s, -c ** 2, -c * s],
                                                 [c * s, s ** 2, -c * s, -s ** 2],
                                                 [-c ** 2, -c * s, c ** 2, c * s],
                                                 [-c * s, -s ** 2, c * s, s ** 2]])
        KGelem_list.append(KGelem)

        for j in range(0, 4):
            for k in range(0, 4):
                K_global[dofMatrix[i, j], dofMatrix[i, k]] += float('%.5g' % KGelem[j, k])

    return K_global, KGelem_list


# Apply kinematic constraints
def constraints(new):
    for i in range(len(new)):
        if 'PinDOFS' in new[i]:
            l = int(new[i][1])
        if 'DOFS with Pins' in new[i]:
            pins = []
            for j in range(0, l):
                pins.append(int(new[i][j + 1]))
    return pins


# Use explicit method
def explicitMethod(KG, pins, numNodes):
    for i in pins:
        for j in range(0, numNodes * 2):
            KG[i, j] = 0
        for k in range(0, numNodes * 2):
            KG[k, i] = 0
        KG[i, i] = 1
    return KG


# Call properties of structure
def properties(new):
    for i in range(len(new)):
        if 'Elastic' in new[i]:  # Identify position in list
            E = float(new[i][1])
        if 'Area' in new[i]:
            A = float(new[i][1])
    return E, A


# Write force matrix
def forceMatrix(new, numNodes):
    F = []
    for i in range(len(new)):
        if 'ForceMatrix' in new[i]:  # Identify position in list
            for j in range(numNodes * 2):
                F.append([float(new[i][j + 1])])
    return F


# Determine stresses per element (via JDW method)
def stress(thetas, sol, KGelem_list, dofs, A):
    els = []
    stress = []
    for i in range(len(dofs)):
        el = []
        c = cos(thetas[i])
        s = sin(thetas[i])

        for j in dofs[i]:
            el.append([sol[j]])
        els.append(el)

        F = Matrix(KGelem_list[i]) * Matrix(el)

        # T-rotational
        Trot = Matrix([[c, s, 0, 0],
                       [-s, c, 0, 0],
                       [0, 0, c, s],
                       [0, 0, -s, c]])
        # Get stress for each element
        fin = Trot * F / A

        fin = np.array(fin)

        for i in range(0, 4):
            stress_e = []
            stress_e.append(fin[2])
        stress.append(stress_e)

    return np.array(stress)


# Begin reading from file
user = input('Problem? (Enter 1, 2, or 3): ')  # Ask user for problem number
if user == '1':
    print('Please enjoy problem 1!')
    file = 'P1.txt'
elif user == '2':
    print('Please enjoy problem 2!')
    file = 'P2.txt'
elif user == '3':
    print('Please enjoy problem 3!')
    file = 'P3.txt'
else:
    print('Not a valid input. Please try again at your earliest convenience.')

# Open file and begin calling from functions
new = read(file)
numNodes, nodeinfo = node(new)
numElems, el_con, lengths, thetas = element(new)
E, A = properties(new)

dofs = nodalDOF(new, numElems)

KG, KGelem_list = KGlobal(E, A, numNodes, numElems, dofs, lengths, thetas)

pins = constraints(new)

# Apply the explicit method
KGe = Matrix(explicitMethod(KG, pins, numNodes))

forceMatrix = Matrix(forceMatrix(new, numNodes))

# Invert K_global
inv = KGe.inv()
# Multiply inverse by force matrix to get final solution
sol = inv * forceMatrix
# pprint(sol)

# Acquire stresses
stress = stress(thetas, sol, KGelem_list, dofs, A)


# Plotting
def plot(s, numNodes, e, n):
    sol = np.array(s)
    # Turn back into list
    q = []
    for i in range(0, numNodes * 2):
        q.append(sol[i][0])

    elems_info = e
    nodes_info = n

    vis = show_truss(q, numElems, elems_info, nodes_info)
    return vis


pl = plot(sol, numNodes, el_con, nodeinfo)

############################################

# Write to a separate file for convenience
if user == '1':
    f = open("P1_results.txt", "w")
    f.write('Results:\n')
    f.write('\n')
    f.write('2-D Displacements:\n')
    f.write(str(np.array(sol)))
    f.write('\n')
    f.write('\n')
    f.write('\n')
    f.write('Axial Stress:\n')
    f.write(str(stress))
    f.close()
    print()
    print('Remember to check the output / results file, as well as the plot. Goodbye!')
elif user == '2':
    f = open("P2_results.txt", "w")
    f.write('Results:\n')
    f.write('\n')
    f.write('2-D Displacements:\n')
    f.write(str(np.array(sol)))
    f.write('\n')
    f.write('\n')
    f.write('\n')
    f.write('Axial Stress:\n')
    f.write(str(stress))
    f.close()
    print()
    print('Remember to check the output / results file, as well as the plot. Goodbye!')
elif user == '3':
    f = open("P3_results.txt", "w")
    f.write('Results:\n')
    f.write('\n')
    f.write('2-D Displacements:\n')
    f.write(str(np.array(sol)))
    f.write('\n')
    f.write('\n')
    f.write('\n')
    f.write('Axial Stress:\n')
    f.write(str(stress))
    f.close()
    print()
    print('Remember to check the output / results file, as well as the plot. Goodbye!')
