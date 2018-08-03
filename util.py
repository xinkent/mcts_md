import numpy as np
import pandas as pd
from graphviz import Graph

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


