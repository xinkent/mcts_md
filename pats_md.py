import numpy as np
import os, glob
from math import *
from util import *
import random
import copy
from graphviz import Graph
import argparse
import dill
MAX_child = 3
MAX_try = 3
rmsd_list = []

class Node:
    def __init__(self, move = None, parent = None, state = None, c = 0.02, depth = 0):
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
            os.system('gmx grompp -f md.mdp -c 0_320.gro -t 0_320.cpt -p topol.top -o %s.tpr -maxwarn 5' % tmp)
        else:
            os.system('gmx grompp -f md.mdp -t md_%d.trr -o %s.tpr -c md_%d.gro -maxwarn 5' %(pstate, tmp, pstate))
        # os.system('gmx mdrun -deffnm %s' % tmp) # pstate.trrからmdrun
        os.system('gmx mdrun -deffnm %s -ntmpi 1 -ntomp 6 -dlb auto -gpu_id 0' % tmp)
        os.system("echo 4 4 | gmx rms -s target_npt_320.gro -f %s.trr  -o rmsd_%d.xvg -tu ns" % (tmp, state)) # rmsdを測定
        rmsds = np.array(read_rmsd('rmsd_%d.xvg'%state))
        # 初期RMSDを書き込み
        if pstate == 0:
            first_rmsd = rmsds[0]
            o = open('log_mcts.txt','w')
            o.write(str(first_rmsd) + '\n')
            o.close()

        min_rmsd = np.min(rmsds)
        min_i = rmsds.argsort()[0]
        os.system('echo 0 | gmx trjconv -s %s.tpr -f %s.trr -o md_%s.trr -e %d' % (tmp, tmp, state, min_i)) # 最小値までのトラジェクトリーを切り出し
        os.system('echo 0 | gmx trjconv -s %s.tpr -f %s.trr -o md_%s.gro -e %d -b %d' % (tmp, tmp, state, min_i, min_i))
        for file in glob.glob("%s*" % tmp) + glob.glob("*#"):
            os.remove(file)
        self.rmsd = min_rmsd
        self.rmsd_depth_dict[self.depth] = min_rmsd
        return min_rmsd

    def prog_widenning(self):
        return sqrt(self.visits) < (3/2) * len(self.childNodes)


# 類似構造でかつ、rmsdが小さい構造がすでにある場合はFalse
def check_similarity(nd):
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
    parser = argparse.ArgumentParser()
    parser.add_argument('--steps',     '-s',  type = int,   default = 1000)
    parser.add_argument('--c',         '-c',  type = float, default = 0.05)
    parser.add_argument('--interrupt', '-ir', type = int,   default = 0)
    parser.add_argument('--continue',  '-cn', type = int,   default = 0)
    args = parser.parse_args()
    steps     = args.steps
    c         = args.c
    interrupt = args.interrupt
    cn  = args.continue
    
    if cn:
        dill.load_session('session.pkl')
    else:
        rootnode = Node(state = rootstate)
        n_state = rootstate
        max_rmsd = -10000
        max_node    = rootnode
        o = open('log_mcts.txt','w')
        o.close()
    for i in steps:
        o = open('log_mcts.txt','a')
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
        node = node.MakeChild(s = state, d = depth)
        result = node.MDrun()
        if result < parent_rmsd: # RMSDが減少した場合のみexpandする
            if check_similarity(node):
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
        if i % 10 == 0:
            G = Graph(format='svg')
            G.attr('node', shape='circle')
            G.graph_attr.update(size="1200")
            make_graph(G,rootnode)
            G.render('./tree/tree_' + str(i) + 'epoch')
        if max_rmsd > -0.1:
            break

    # bestなノードまでのトラジェクトリを結合
    if not interrupt:
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
        os.system("gmx trjcat -f " + trjs + " -o merged_mcts.trr -cat")
        os.system("echo 4 4 | gmx rms -s target_npt.gro -f merged_mcts.trr -tu ns -o rmsd_mcts_tmp.xvg")
        modify_rmsd('rmsd_mcts_tmp.xvg', 'rmsd_mcts.xvg')
        os.remove('rmsd_mcts_tmp.xvg')
        # for file in (glob.glob("*#") + glob.glob("md_*") + glob.glob("rmsd_[0-9]*")):
        #     os.remove(file)
        G = Graph(format='svg')
        G.attr('node', shape='circle')
        G.graph_attr.update(size="1200")
        make_graph(G,rootnode)
        G.render('./tree/mcts_tree')

    # 途中経過を保存(→どうやる？)
    else:
        dill.dump_session('session.pkl')




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
    os.system('rm all_structure.gro')
    UCT(0, 5000, verbose=True)