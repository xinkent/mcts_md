import numpy as np
import os, glob
from math import *
from util import *
import random
import copy
from graphviz import Graph
import argparse
import pickle
MAX_child = 3
MAX_try = 3
MIN_RMSD = 0.1

parser = argparse.ArgumentParser()
parser.add_argument('--reactant',  '-r',                   default = '0')
parser.add_argument('--target',    '-t',                   default = 'target_processed')
parser.add_argument('--topol',     '-top',                 default = 'topol')
parser.add_argument('--steps',     '-s',     type = int,   default = 1000)
parser.add_argument('--c',         '-c',     type = float, default = 0.05)
parser.add_argument('--continue_', '-cn',    type = int,   default = 0)
parser.add_argument('--ntmpi',     '-ntmpi', type = int,   default = 1)
parser.add_argument('--ntomp',     '-ntomp', type = int,   default = 10) 
parser.add_argument('--delete',    '-del'  , type = int,   default = 0)
args = parser.parse_args()
reactant = args.reactant
target   = args.target
c_       = args.c
topol    = args.topol
ntmpi    = args.ntmpi
ntomp    = args.ntomp
delete   = args.delete
class Node:
    def __init__(self, move = None, parent = None, state = None, c = c_, depth = 0):
        self.parentNode = parent # "None" for the root node
        self.childNodes = []
        self.depth = depth
        self.child_depth = 0
        self.rmsd_sum = 0
        self.rmsd_max = -1000
        self.rmsd_depth_dict = {}
        self.visits = 0
        self.state = state
        self.untriedMoves = MAX_child
        self.c = c
        self.alpha = 1.5
        self.rmsd = 1000 # 仮
        self.try_num = 0

    def UCTSelectChild(self):
        child_rmsds = [ch.rmsd_max for ch in self.childNodes]
        rmsd_diff = max(child_rmsds) - min(child_rmsds)
        child_depths = [ch.child_depth for ch in self.childNodes]
        # depth_diff = max(child_depths) - min(child_depths)
        c = rmsd_diff * self.alpha
        s = sorted(self.childNodes, key = lambda ch: ch.rmsd_max + self.c * sqrt(2*log(self.visits)/ch.visits))[-1] # 通常盤
        # s = sorted(self.childNodes, key = lambda ch: ch.rmsd_max + c * sqrt(2*log(self.visits)/ch.visits))[-1] # RMSDの差に比例してCを定める方式
        # s = sorted(self.childNodes, key = lambda ch: ch.rmsd_max + self.c * (depth_diff/2) * sqrt(2*log(self.visits)/ch.visits))[-1]
        return s

    def UCTSelectChildByDepth(self):
        node_idx = np.arange(len(self.childNodes))
        while True:
            print('node_idx  ' + str(node_idx))
            child_depths = np.array([self.childNodes[idx].child_depth for idx in node_idx])
            min_depth = min(child_depths)
            min_depth_idx = np.where(child_depths == min(child_depths))[0]
            ucts = [self.childNodes[idx].rmsd_depth_dict[min_depth] + self.c * sqrt(2*log(self.visits)/self.childNodes[idx].visits) for idx in node_idx]
            max_uct_idx = np.argmax(ucts)
            if max_uct_idx in  min_depth_idx:
                return self.childNodes[node_idx[max_uct_idx]]
            else:
                node_idx = np.delete(node_idx, min_depth_idx)


    def CalcUCT(self):
        pnd = self.parentNode
        if pnd == None:
            return -1
        child_rmsds = [ch.rmsd_max for ch in pnd.childNodes]
        child_depths = [ch.child_depth for ch in pnd.childNodes]
        min_depth = min(child_depths)

        # rmsd_diff = max(child_rmsds) - min(child_rmsds)
        # c = rmsd_diff * self.alpha
        # uct = self.rmsd_max + self.c * sqrt(2*log(pnd.visits) / self.visits)
        uct = self.rmsd_depth_dict[min_depth] + self.c * sqrt(2*log(pnd.visits) / self.visits)
        return uct

    def MakeChild(self, s, d):
        n = Node(parent = self, state = s, c = self.c, depth = d)
        return n

    def AddChild(self, n):
        # n = Node(parent = self, state = s)
        self.untriedMoves -= 1
        self.childNodes.append(n)
        return n

    def DeleteChild(self, n):
        delete_i = -1
        for ch_i, ch in enumerate(self.childNodes):
            if ch.state == n.state:
                delete_i = ch_i
                break
        self.childNodes.pop(delete_i)
        # self.untriedMoves += 1

    def SearchMaxRmsd(self):
        child_rmsds = [ch.rmsd_max for ch in self.childNodes]
        if child_rmsds != []:
            self.max_rmsd = max(child_rmsds)
        else:
            self.max_rmsd = -1000

    def Update(self, result, d):
        self.visits += 1
        self.rmsd_sum += result
        if result > self.rmsd_max:
            self.rmsd_max = result
        if d > self.child_depth:
            self.child_depth = d
            self.rmsd_depth_dict[d] = result
        else:
            if result > self.rmsd_depth_dict[d]:
                self.rmsd_depth_dict[d] = result

    def MDrun(self):
        state = self.state
        pstate = self.parentNode.state
        self.parentNode.try_num += 1
        tmp = str(state) + '_tmp'
        if pstate == 0:
            os.system('gmx grompp -f md.mdp -c %s.gro -t %s.cpt -p %s.top -o %s.tpr -maxwarn 5' % (reactant, reactant, topol, tmp))
        else:
            os.system('gmx grompp -f md.mdp -t md_%d.trr -o %s.tpr -c md_%d.gro -maxwarn 5' %(pstate, tmp, pstate))
        # os.system('gmx mdrun -deffnm %s' % tmp) # pstate.trrからmdrun
        os.system('gmx mdrun -deffnm %s  -ntmpi %d  -ntomp %d -dlb auto' % (tmp,  ntmpi, ntomp))

        os.system("echo 4 4 | gmx rms -s %s.gro -f %s.trr  -o rmsd_%d.xvg -tu ns" % (target, tmp, state)) # rmsdを測定
        rmsds = np.array(read_rmsd('rmsd_%d.xvg'%state))
        # 初期RMSDを書き込み
        if pstate == 0:
            first_rmsd = rmsds[0]
            o = open('log_pats.txt','w')
            o.write(str(first_rmsd) + '\n')
            o.close()

        min_rmsd = np.min(rmsds)
        min_i = rmsds.argsort()[0]
        os.system('echo 0 | gmx trjconv -s %s.tpr -f %s.trr -o md_%s.trr -e %d' % (tmp, tmp, state, min_i)) # 最小値までのトラジェクトリーを切り出し
        os.system('echo 0 | gmx trjconv -s %s.tpr -f %s.trr -o md_%s.gro -e %d -b %d' % (tmp, tmp, state, min_i, min_i))
        for file in glob.glob("*#"):
            os.remove(file)
        for ext in ['trr', 'tpr', 'edr', 'log','gro', 'cpt']:
            for file in glob.glob('%s.%s' % (tmp,ext)):
                os.remove(file)
        if delete == 1:
            for file in glob.glob('%s*' % tmp):
                os.remove(file)
        self.rmsd = min_rmsd
        self.rmsd_depth_dict[self.depth] = min_rmsd
        return min_rmsd

    def prog_widenning(self):
        return sqrt(self.visits) < (3/2) * len(self.childNodes)


# 類似構造でかつ、rmsdが小さい構造がすでにある場合はFalse
def check_similarity(nd, rmsd_list):
    if nd.state == 1:
        return True
    os.system('echo 4 4 |gmx rms -s md_%s.gro -f all_structure.gro -o rmsd_tmp.xvg'%(nd.state))
    rmsd_tmp = np.array(read_rmsd('rmsd_tmp.xvg'))
    print(rmsd_list)
    for file in glob.glob("*#"):
        os.remove(file)
    if any((rmsd_tmp < 0.05) & (np.array(rmsd_list) < nd.rmsd)):
        return False
    else:
        return True


def UCT(rootstate):
    steps     = args.steps
    c_        = args.c
    cn        = args.continue_
    
    succeed = 0 
    if cn:
        with open('vars.pickle','rb') as f:
            var_list = pickle.load(f)
        rootnode = var_list[0]
        n_state = var_list[1]
        max_rmsd = var_list[2]
        max_node = var_list[3]
        rmsd_list = var_list[4]
            
    else:
        os.system('rm all_structure.gro')
        rootnode = Node(state = rootstate, c = c_)
        n_state = rootstate
        max_rmsd = -10000
        max_node    = rootnode
        rmsd_list = []
        o = open('log_pats.txt','w')
        o.close()
    for i in range(steps):
        o = open('log_pats.txt','a')
        n_o = open('near_count.txt', 'a')
        node = rootnode
        state = rootstate
        flag = 0
        # Select
        # while (node.untriedMoves == 0 or node.prog_widenning()) and node.childNodes != []: # progressive widennning
        while (node.untriedMoves == 0 or node.try_num >= MAX_try) and node.childNodes != []: # node is fully expanded and non-terminal
            # node = node.UCTSelectChildByDepth()
            node = node.UCTSelectChild()
            state = node.state

        # Expand
        parent_node =node
        parent_rmsd = node.rmsd
        parent_depth = node.depth
        state = n_state + 1
        depth = parent_depth + 1
        # node = node.AddChild(state) # add child and descend tree
        print('state is ' + str(state))	
        node = node.MakeChild(s = state, d = depth)
        result = node.MDrun()
        if result < parent_rmsd: # RMSDが減少した場合のみexpandする
            if check_similarity(node, rmsd_list):
                parent_node.AddChild(node)
                os.system('cat md_%s.gro >> all_structure.gro'%state) # 構造を保存
                rmsd_list.append(result) # RMSDを保存
                n_state += 1
            else:
                n_o.write(str(state) + '\n')
        # MAX_tryやっても追加されない場合は枝を刈り取る
        # else:
        #     flag = 1

        n_o.close()
        # Backpropagate
        result = -1 * result
        if result > max_rmsd:
            max_rmsd = result
            max_node = node
        while node != None:
            node.Update(result, depth)
            node = node.parentNode

        # if flag == 1:
        #     tmp_node = parent_node
        #     while tmp_node.try_num >= MAX_try and tmp_node.childNodes == []:
        #         tmp_parent = tmp_node.parentNode
        #         if tmp_parent == None:
        #             break
        #         tmp_parent.DeleteChild(tmp_node)
        #         tmp_parent.SearchMaxRmsd()
        #         tmp_node = tmp_parent

        o.write(str(-1 * max_rmsd) + '\n')
        o.close()
        if i % 100 == 0:
            G = Graph(format='svg')
            G.attr('node', shape='circle')
            G.graph_attr.update(size="1200")
            make_graph(G,rootnode)
            G.render('./tree/tree_' + str(i) + 'epoch')
        if max_rmsd > -MIN_RMSD:
            succeed = 1
            break

    # bestなノードまでのトラジェクトリを結合
    if succeed:
        trjs = ""
        node = max_node
        while True:
            print(node.state)
            trjs = "md_" + str(node.state) + ".trr " + trjs
            node = node.parentNode
            if node.parentNode == None:
                break
        o_trj = open('trjs.txt', 'w')
        o_trj.write(trjs)
        o_trj.close()
        os.system("gmx trjcat -f " + trjs + " -o merged_pats.trr -cat")
        os.system("echo 4 4 | gmx rms -s %s.gro -f merged_pats.trr -tu ns -o rmsd_pats_tmp.xvg" % target)
        modify_rmsd('rmsd_pats_tmp.xvg', 'rmsd_pats.xvg')
        os.remove('rmsd_pats_tmp.xvg')
        # for file in (glob.glob("*#") + glob.glob("md_*") + glob.glob("rmsd_[0-9]*")):
        #     os.remove(file)
        G = Graph(format='svg')
        G.attr('node', shape='circle')
        G.graph_attr.update(size="1200")
        make_graph(G,rootnode)
        G.render('./tree/pats_tree')

    # 成功しなかった場合は、途中経過をpickleに保存
    else:
        var_list = []
        var_list.append(rootnode)
        var_list.append(n_state)
        var_list.append(max_rmsd)
        var_list.append(max_node)
        var_list.append(rmsd_list)
        with open('vars.pickle', mode = 'wb') as f:
            pickle.dump(var_list, f)



def make_graph(G, nd):
    state = nd.state
    uct = nd.CalcUCT()
    G.node(str(state), str(state) + '\n' + "{:.4}".format(float(nd.rmsd))  + '\n' + str(nd.visits) + '\n' + str(nd.child_depth) + '\n' + str(uct))
    parent_node = nd.parentNode
    if parent_node != None:
        parent_state = parent_node.state
        G.edge(str(parent_state), str(state))
    for child_node in nd.childNodes:
        make_graph(G,child_node)

if __name__ == "__main__":
    UCT(0)
