import numpy as np
import os, glob
from multiprocessing import Pool
from util import *

MIN_RMSD = 0.1

def first_cycle(n):
    # print(n)
    os.system('gmx grompp -f md.mdp -c 0.gro -t 0.cpt -p topol.top -o md_0_%d.tpr -maxwarn 5' % n)
    os.system('gmx mdrun -deffnm md_0_%d -ntmpi 1 -ntomp 10 -dlb auto' % n)
    os.system("echo 4 4 | gmx rms -s target_processed.gro -f md_0_%d.trr  -o rmsd_0_%d.xvg -tu ns" % (n, n))

def md_cycle(args):
    step = args[0]
    n = args[1]
    md_ind = args[2]
    last = args[3]
    input  = 'md_%d_%d' % (step-1, md_ind[0])
    temp   = 'res_%d_%d' % (step-1, n)
    md     = 'md_%d_%d' % (step, n)
    rmsd   = 'rmsd_%d_%d' % (step, n)
    # 全ステップで計算したトラジェクトリを指定の部分まで切り取る
    os.system('echo 0 | gmx trjconv -s %s.tpr -f %s.trr -o %s.trr -e %d' % (input, input, temp, md_ind[1]))
    os.system('echo 0 | gmx trjconv -s %s.tpr -f %s.trr -o %s.gro -e %d' % (input, input, temp, md_ind[1]))
    # 切り取った箇所から新たにMD
    if not last:
        os.system('gmx grompp -f md.mdp -t %s.trr -o %s.tpr -c %s.gro -maxwarn 5' %(temp, md, temp))
        os.system('gmx mdrun -deffnm %s -ntmpi 1 -ntomp 10 -dlb auto' % md)
        os.system("echo 4 4 | gmx rms -s target_processed.gro -f %s.trr  -o %s.xvg -tu ns" % (md, rmsd))

def pacs_md(MAX_CYCLE = 100, n_para = 5, continue_flag = True):
    dt = 1 # ps for each step
    nsteps = 101 # ショートMDのステップ数+1(時刻0を含む)
    if continue_flag:
        edge_log = np.loadtxt("edge_log.csv", delimiter = ",")
        edge_log = [list(l) for l in log]
        edge_log = log[0:len(log)-1]
        cycle_step = len(edge_log)
        cycle_num = -1 # 現実行時におけるステップ数
    else:
        cycle_step = 0
        o = open('log_pacs.txt','w')
        o.close()
        edge_log = []  # Treeの遷移関係のlog
        cycle_num = 0
    min_rmsd = 10000 # 初期値

    while cycle_num < MAX_CYCLE and min_rmsd >= MIN_RMSD:
        if not continue_flag:
            if cycle_step == 0:
                for i in range(n_para):
                    first_cycle(i)
            else:
                for i in range(n_para):
                    md_cycle([cycle_step, i, min_rmsd_idx[i],False])
        else:
            continue_flag = False
        # 最小RMSDを記録
        rmsd_list = np.zeros(n_para * nsteps)
        for i in range(n_para):
            rmsd_list[i * nsteps : (i+1) * nsteps] = read_rmsd('rmsd_%d_%d.xvg'%(cycle_step, i))
        min_rmsd = min(rmsd_list)

       # 最初の一回だけ初期値を記録
        if cycle_step == 0:
            first_rmsd = rmsd_list[0]
            write_log(first_rmsd)
        write_log(min_rmsd)
        # Treeの遷移関係を記録
        min_rmsd_idx = rmsd_list.argsort()[:n_para]
        min_rmsd_idx = [(r // nsteps, r % nsteps) for r in min_rmsd_idx]
        edge_log.append([r[0] for r in min_rmsd_idx])

        cycle_step += 1
        cycle_num  += 1

        # 不要なファイルを削除
        # for file in glob.glob('md_[0-9]*') + glob.glob("*#"):
        #     os.remove(file)

    # 最後のMDのトラジェクリについて切り取り処理をする
    for i in range(n_para):
        md_cycle([cycle_step, i, min_rmsd_idx[i],True])
    # logをファイルに保存(concat_traj)に渡す
    np.savetxt('edge_log.csv', np.array(edge_log), delimiter=',')

    if min_rmsd < MIN_RMSD:
        concat_traj()
        make_tree_pacs('edge_log.csv')

# log.csvを元にshort MDのトラジェクリを繋げる
def concat_traj():
    log = np.loadtxt("edge_log.csv", delimiter = ",")
    back_list = [0]
    step = log.shape[0] - 1
    log_index = 0
    # 最もRMSDが小さくなったものについて、backtrackする
    for i in range(step,0, -1):
        log_index = int(log[i, log_index])
        back_list.insert(0, log_index)
    # print(back_list)
    # トラジェクトリを結合
    trajs = []
    for i, j in enumerate(back_list):
        traj = "res_%d_%d"%(i,j)
        print(traj)
        trajs.append(traj)
    traj_str = ""
    for traj in trajs:
        traj_str += (traj + " ")
    # print(traj_str)
    os.system( "gmx trjcat -f " + traj_str + " -o merged_pacs.trr -cat")
    os.system("echo 4 4| gmx rms -s target_processed.gro -f merged_pacs.trr -tu ns -o rmsd_pacs_tmp.xvg")
    modify_rmsd('rmsd_pacs_tmp.xvg', 'rmsd_pacs.xvg') # ’ダブり’があるので除去
    os.remove('rmsd_pacs_tmp.xvg')

# 各ステップでのrmsdを記録
def write_log(rmsd):
    o = open('log_pacs.txt','a')
    o.write(str(rmsd) + '\n')
    o.close()


if __name__ == '__main__':
    M = 1000 # サイクル数
    k = 5 # 並列数
    pacs_md(M, k)
    # for file in (glob.glob("*#") + glob.glob("md_[0-9]*") + glob.glob("res_[0-9]*") + glob.glob("rmsd_[0-9]*")):
    #     os.remove(file)
