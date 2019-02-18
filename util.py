import sys
import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
import pickle
from pats_md import *
# from pylab import *
from graphviz import Graph
import mdtraj as md
from itertools import combinations

def read_rmsd(name):
    f = open(name)
    rmsd = []
    while True:
        l = f.readline()
        if l[0] != '#' and l[0] != '@':
            break
    while l != '':
        rmsd.append(float(l.split()[1]))
        l = f.readline()
    f.close()
    return rmsd

# modify the rmsd.xvg in concatenated trajectory
def modify_rmsd(ip, op):
    f = open(ip)
    o = open(op, 'w')
    while True:
        l = f.readline()
        if l[0] != '#' and l[0] != '@':
            break
        else:
            o.write(l)

    cum_time = 0.0
    rmsd = l.split()[1]
    o.write(str(cum_time) + "    " + str(rmsd) + "\n")
    l = f.readline()
    cum_time += 1.0

    while l != '':
        time = float(l.split()[0])
        rmsd = l.split()[1]
        if time != 0.0:
            o.write(str(cum_time) + "    " + str(rmsd) + "\n")
            cum_time += 1.0
        l = f.readline()
    f.close()
    o.close()

def make_tree_pacs(log_file):
    log = pd.read_csv(log_file, header = None)
    log = np.array(log)
    G = Graph(format='svg')
    G.attr('node', shape='circle')
    G.graph_attr.update(size="320")
    G.node('0-0', '0-0')
    for i in range(5):
        state = '1-' + str(i)
        G.node(state, state)
        G.edge('0-0', state)

    for i in range(len(log)):
        log_i = log[i]
        for j, l in enumerate(log_i):
            state = str(i + 2) + '-' + str(j)
            pstate = str(i + 1) + '-' + str(int(l))
            G.node(state, state)
            G.edge(pstate, state)
    G.render('tree_pacs')


def draw_pats_tree(pkl,native, traj, out):
    fractions = frac_native_contacts(native, traj)
    with open(pkl, 'rb') as f:
        var_list = pickle.load(f)
        rootnode = var_list[0]
    num_state = dfs(rootnode)
    print('num_state  ' + str(num_state))
    
    G = Graph(format='png')
    G.attr("node", shape="circle", style="filled")
    G.graph_attr.update(size="30")
    make_graph(G, rootnode, fractions)
    G.render(out)

def dfs(nd):
    state = nd.state
    if len(nd.childNodes) == 0:
        return 1
    return sum([dfs(ch) for ch in nd.childNodes]) + 1    

def make_graph(G,nd, fractions):
    state = nd.state

    # r,g,b = 255 - int(255/N * state), 0, int(255/N * state)
    # r,g,b = 255 - int(255/N * state), 255, 255
    v = 255 / (max(fractions) - min(fractions)) * (fractions[state] - min(fractions))
    r,g,b, _  = [int(c * 255) for c in cm.cool(int(v))]
    r_hex,g_hex,b_hex = hex(r), hex(g), hex(b)
    color_hex = "#" + r_hex[2:].zfill(2) + g_hex[2:].zfill(2) + b_hex[2:].zfill(2)
    G.node(str(state), fillcolor=color_hex)
    parent_node = nd.parentNode
    if parent_node != None:
        parent_state = parent_node.state
        G.edge(str(parent_state), str(state))
    for child_node in nd.childNodes:
        make_graph(G,child_node, fractions)

def best_hummer_q(traj, native):
    """Compute the fraction of native contacts according the definition from
    Best, Hummer and Eaton [1]

    Parameters
    ----------
    traj : md.Trajectory
        The trajectory to do the computation for
    native : md.Trajectory
        The 'native state'. This can be an entire trajecory, or just a single frame.
        Only the first conformation is used

    Returns
    -------
    q : np.array, shape=(len(traj),)
        The fraction of native contacts in each frame of `traj`

    References
    ----------
    ..[1] Best, Hummer, and Eaton, "Native contacts determine protein folding
          mechanisms in atomistic simulations" PNAS (2013)
    """

    BETA_CONST = 50  # 1/nm
    LAMBDA_CONST = 1.8
    NATIVE_CUTOFF = 0.45  # nanometers

    # get the indices of all of the heavy atoms
    heavy = native.topology.select_atom_indices('heavy')
    # get the pairs of heavy atoms which are farther than 3
    # residues apart
    heavy_pairs = np.array(
        [(i,j) for (i,j) in combinations(heavy, 2)
            if abs(native.topology.atom(i).residue.index - \
                   native.topology.atom(j).residue.index) > 3])

    # compute the distances between these pairs in the native state
    heavy_pairs_distances = md.compute_distances(native[0], heavy_pairs)[0]
    # and get the pairs s.t. the distance is less than NATIVE_CUTOFF
    native_contacts = heavy_pairs[heavy_pairs_distances < NATIVE_CUTOFF]
    print("Number of native contacts", len(native_contacts))

    # now compute these distances for the whole trajectory
    r = md.compute_distances(traj, native_contacts)
    # and recompute them for just the native state
    r0 = md.compute_distances(native[0], native_contacts)

    q = np.mean(1.0 / (1 + np.exp(BETA_CONST * (r - LAMBDA_CONST * r0))), axis=1)
    return q 
   
def frac_native_contacts(n, t):
    native = md.load(n)
    # t = md.load('../merged_pats.gro')
    traj = md.load(t)
    q = best_hummer_q(traj, native)
    return q
