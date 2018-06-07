import numpy as np
import os, glob
from multiprocessing import Pool
from util import *

def first_cycle(n):
    print(n)
    os.system('gmx_mpi grompp -f md.mdp -c 0.gro -t 0.cpt -p topol.top -o md_0_%d.tpr -maxwarn 5' % n)
    os.system('gmx_mpi mdrun -deffnm md_0_%d' % n)
    os.system("echo 4 4 | gmx_mpi rms -s target_npt.gro -f md_0_%d.trr  -o rmsd_0_%d.xvg -tu ns" % (n, n))

def md_cycle(args):
    step = args[0]
    n = args[1]
    md_ind = args[2]
    input  = 'md_%d_%d' % (step-1, md_ind[0])
    temp = 'res_%d_%d' % (step-1, n)
    md    = 'md_%d_%d' % (step, n)
    os.system('echo 0 | gmx_mpi trjconv -s %s.tpr -f %s.trr -o %s.trr -e %d' % (input, input, temp, md_ind[1]))
    os.system('echo 0 | gmx_mpi trjconv -s %s.tpr -f %s.trr -o %s.gro -e %d' % (input, input, temp, md_ind[1]))
    os.system('gmx_mpi grompp -f md.mdp -t %s.trr -o %s.tpr -c %s.gro -maxwarn 5' %(temp, md, temp))
    os.system('gmx_mpi mdrun -deffnm %s' % md)
    os.system("echo 4 4 | gmx_mpi rms -s target_npt.gro -f %s.trr  -o rmsd_%d_%d.xvg -tu ns" % (md, step, n))

def pacs_md(MAX_CYCLE = 100, n_para = 5):
    dt = 1 # ps for each step
    nsteps = 101 # number of steps for each md including inital state
    pass_flag = False
    cycle_step = 0
    log = np.zeros((MAX_CYCLE - cycle_step, n_para))  # trajectory間の接続を記録
    log_i = 0
    # 途中(100回目)から再開する場合
    # pass_flag = True
    # cycle_step = 100

    while cycle_step < MAX_CYCLE:
        if not pass_flag:
            if cycle_step == 0:
                for i in range(n_para):
                    first_cycle(i)
            else:
                for i in range(n_para):
                    md_cycle([cycle_step, i, min_rmsd[i]])
        else:
            pass_flag = False
        rmsd_list = np.zeros(n_para * nsteps)
        for i in range(n_para):
            rmsd_list[i * nsteps : (i+1) * nsteps] = read_rmsd('rmsd_%d_%d.xvg'%(cycle_step, i))
        min_rmsd = rmsd_list.argsort()[:n_para]
        min_rmsd = [(r // nsteps, r % nsteps) for r in min_rmsd]
        log[log_i,:] = [r[0] for r in min_rmsd]
        cycle_step += 1
        log_i += 1

    np.savetxt('log_0_100.csv', log, delimiter=',')

def concat_traj():
    log = np.loadtxt("log_0_100.csv", delimiter = ",")
    back_list = []
    step = log.shape[0] - 1
    log_index = 0

    for i in range(step,0, -1):
        log_index = int(log[i, log_index])
        back_list.insert(0, log_index)
    # print(back_list)
    trajs = []
    for i, j in enumerate(back_list):
        traj = "res_%d_%d"%(i,j)
        print(traj)
        trajs.append(traj)

    traj_str = ""
    for traj in trajs:
        traj_str += (traj + " ")

    # print(traj_str)
    os.system( "gmx_mpi trjcat -f " + traj_str + " -o merged_pacs.trr ")
    os.system("echo 4 4| gmx_mpi rms -s target_npt.gro -f merged_pacs.trr -tu ns -o rmsd_pacs_tmp.xvg")
    modify_rmsd('rmsd_pacs_tmp.xvg', 'rmsd_pacs.xvg')
    os.remove('rmsd_pacs_tmp.xvg')

# 各ステップでのrmsdを記録
def write_log(m):
    o = open('log_pacs.txt','w')
    for i in range(m):
        rmsd = read_rmsd('rmsd_%d_%d.xvg'%(i,0))
        o.write(str(rmsd[0]) + '\n')
    o.close()


if __name__ == '__main__':
    M = 3 # サイクル数
    k = 2 # 並列数
    # pacs_md(M, k)
    concat_traj()
    write_log(M)
    for file in (glob.glob("*#") + glob.glob("md_[0-9]*") + glob.glob("res_[0-9]*") + glob.glob("rmsd_[0-9]*")):
        os.remove(file)    


