import numpy as np
import os, glob
from math import *
from util import *
import random

MAX_child = 5

class Node:
    def __init__(self, move = None, parent = None, state = None, depth = None):
        self.parentNode = parent # "None" for the root node
        self.childNodes = []
        self.depth = 0
        self.rmsd = 0
        self.visits = 0
        self.state = state
        self.untriedMoves = MAX_child
        self.c = 1/50 

    def UCTSelectChild(self):
        s = sorted(self.childNodes, key = lambda c: c.rmsd/c.visits + self.c * sqrt(2*log(self.visits)/c.visits))[-1]
        return s

    def AddChild(self, s):
        n = Node(parent = self, state = s)
        self.untriedMoves -= 1
        self.childNodes.append(n)
        return n

    def Update(self, result):
        self.visits += 1
        self.rmsd += result

    def MDrun(self):
        state = self.state
        pstate = self.parentNode.state
        tmp = str(state) + '_tmp'
        if pstate == 0:
            os.system('gmx_mpi grompp -f md.mdp -c 0.gro -t 0.cpt -p topol.top -o %s.tpr -maxwarn 5' % tmp)
        else:
            os.system('gmx_mpi grompp -f md.mdp -t md_%d.trr -o %s.tpr -c md_%d.gro -maxwarn 5' %(pstate, tmp, pstate))
        # os.system('gmx_mpi mdrun -deffnm %s' % tmp) # pstate.trrからmdrun
        os.system('mpirun -np 4 gmx_mpi mdrun -deffnm %s -dlb auto' % tmp)
        os.system("echo 4 4 | gmx_mpi rms -s target_npt.gro -f %s.trr  -o rmsd_%d.xvg -tu ns" % (tmp, state)) # rmsdを測定
        rmsds = np.array(read_rmsd('rmsd_%d.xvg'%state))
        min_rmsd = np.min(rmsds)
        min_i = rmsds.argsort()[0]
        os.system('echo 0 | gmx_mpi trjconv -s %s.tpr -f %s.trr -o md_%s.trr -e %d' % (tmp, tmp, state, min_i)) # 最小値までのトラジェクトリーを切り出し
        os.system('echo 0 | gmx_mpi trjconv -s %s.tpr -f %s.trr -o md_%s.gro -e %d' % (tmp, tmp, state, min_i))
        for file in glob.glob("%s*" % tmp):
            os.remove(file)
        return min_rmsd

    def prog_widenning(self):
        return sqrt(self.visits) < len(self.childNodes)

    def __repr__(self):
        if self.parentNode == None:
                return "[S:" + str(self.state)  + " W/V:" + str(self.rmsd) + "/" + str(self.visits) + " U:" + str(self.untriedMoves) + "]"
        else:
            return "[S:" + str(self.state) + "  UCT:" + str(self.rmsd/self.visits + self.c * sqrt(2 * log(self.parentNode.visits)/self.visits)) + " W/V:" + str(self.rmsd) + "/" + str(self.visits) + " U:" + str(self.untriedMoves) + "]"

    def TreeToString(self, indent):
        s = self.IndentString(indent) + str(self)
        for c in self.childNodes:
             s += c.TreeToString(indent+1)
        return s

    def IndentString(self,indent):
        s = "\n"
        for i in range (1,indent+1):
            s += "| "
        return s

    def ChildrenToString(self):
        s = ""
        for c in self.childNodes:
             s += str(c) + "\n"
        return s


def UCT(rootstate, itermax, verbose = False):
    rootnode = Node(state = rootstate)
    n_state = 0
    max_rmsd = -10000
    max_node    = rootnode
    o = open('log_mcts.txt','w')
    o.close()
    for i in range(itermax):
        o = open('log_mcts.txt','a')
        node = rootnode
        state = rootstate
        # Select
        while node.untriedMoves == 0 and node.childNodes != []: # node is fully expanded and non-terminal
            node = node.UCTSelectChild()
            state = node.state

        # Expand
        if node.untriedMoves != 0: # if we can expand (i.e. state/node is non-terminal)
            state = n_state + 1
            node = node.AddChild(state) # add child and descend tree
            n_state += 1

        # Backpropagate
        result = -1 * node.MDrun()
        # result = get_random()
        if result > max_rmsd:
            max_rmsd = result
            max_node = node
        while node != None: # backpropagate from the expanded node and work back to the root node
            node.Update(result) # state is terminal. Update node with result from POV of node.playerJustMoved
            node = node.parentNode

        o.write(str(-1 * max_rmsd) + '\n')
        o.close()
        # if (verbose): ot.write(rootnode.TreeToString(0))
        # else: ot.write(rootnode.ChildrenToString())
    # return sorted(rootnode.childNodes, key = lambda c: c.visits)[-1].move # return the move that was most visited
    
    trjs = ""
    node = max_node
    while True:
        trjs = "md_" + str(node.state) + ".trr " + trjs
        node = node.parentNode
        if node.parentNode == None:
            break
    o_trj = open('trjs.txt', 'w')
    o_trj.write(trjs)
    o_trj.close()
    os.system("gmx_mpi trjcat -f " + trjs + " -o merged_mcts.trr -cat")
    os.system("echo 4 4 | gmx_mpi rms -s target_npt.gro -f merged.trr -tu ns -o rmsd_mcts_tmp.xvg")
    modify_rmsd('rmsd_mcts_tmp.xvg', 'rmsd_mcts.xvg')
    os.remove('rmsd_mcts_tmp.xvg')
    for file in (glob.glob("*#") + glob.glob("md_*") + glob.glob("rmsd_[0-9]*")):
        os.remove(file)
    # Output some information about the tree - can be omitted
    ot = open('tree.txt','w')
    if (verbose): ot.write(rootnode.TreeToString(0))
    else: ot.write(rootnode.ChildrenToString())
    ot.close()

def UCT_progressive_widenning(rootstate, itermax, verbose = False):
    rootnode = Node(state = rootstate)
    n_state = 0
    max_rmsd = -10000
    max_node    = rootnode
    o = open('log_prog.txt','w')
    ot = open('tree_prog.txt','w')
    o.close()
    for i in range(itermax):
        o = open('log_prog.txt','a')
        ot = open('tree_prog.txt', 'a')
        node = rootnode
        state = rootstate
        # Select
        while (node.untriedMoves == 0 or node.prog_widenning) and node.childNodes != []: # node is fully expanded and non-terminal
            node = node.UCTSelectChild()
            state = node.state

        # Expand
        if node.untriedMoves != 0: # if we can expand (i.e. state/node is non-terminal)
            state = n_state + 1
            node = node.AddChild(state) # add child and descend tree
            n_state += 1

        # Backpropagate
        # result = -1 * node.MDrun()
        result = get_random()
        if result > max_rmsd:
            max_rmsd = result
            max_node = node
        while node != None: # backpropagate from the expanded node and work back to the root node
            node.Update(result) # state is terminal. Update node with result from POV of node.playerJustMoved
            node = node.parentNode

        o.write(str(-1 * max_rmsd) + '\n')
        o.close()
        if (verbose): ot.write(rootnode.TreeToString(0))
        else: ot.write(rootnode.ChildrenToString())
        ot.close()
    # return sorted(rootnode.childNodes, key = lambda c: c.visits)[-1].move # return the move that was most visited

    return
    trjs = ""
    node = max_node
    while True:
        trjs = "md_" + str(node.state) + ".trr " + trjs
        node = node.parentNode
        if node.parentNode == None:
            break
    o_trj = open('trjs.txt', 'w')
    o_trj.write(trjs)
    o_trj.close()
    os.system("gmx_mpi trjcat -f " + trjs + " -o merged.trr -cat")
    os.system("echo 4 4 | gmx_mpi rms -s target_npt.gro -f merged.trr -tu ns -o rmsd_merged_tmp.xvg")
    modify_rmsd('rmsd_merged_tmp.xvg', 'merged_rmsd.xvg')
    # for file in (glob.glob("*#") + glob.glob("md_*") + glob.glob("rmsd_*")):
        # os.remove(file)
    # Output some information about the tree - can be omitted
    if (verbose): ot.write(rootnode.TreeToString(0))
    else: ot.write(rootnode.ChildrenToString())
    ot.close()

 
if __name__ == "__main__":
    UCT(0,10, verbose = True)
    
