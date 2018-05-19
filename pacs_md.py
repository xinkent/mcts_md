#--------------------------------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------------------------------
import numpy as np
# from gromacs import tools
import os, glob
from multiprocessing import Pool

def initialize():
    os.system('gmx_mpi pdb2gmx -f 2lzm.pdb -o 2lzm_processed.gro -ignh -water spce -ff amber99sb')

    # Set the simulation config
    os.system('gmx_mpi editconf -f 2lzm_processed.gro -o 2lzm_newbox.gro -c -d 1.0 -bt cubic')

    # Add solvate
    os.system('gmx_mpi solvate -cp 2lzm_newbox.gro -cs spc216.gro -o 2lzm_solv.gro -p topol.top')

    # Add ions
    os.system('gmx_mpi grompp -f ions.mdp -c 2lzm_solv.gro -p topol.top -o ions.tpr -maxwarn 5')
    os.system("echo 13 | gmx_mpi genion -s ions.tpr -o 2lzm_solv_ions.gro -p topol.top -pname NA -nname CL -nn 8")

    # Energy minimization
    os.system('gmx_mpi grompp -f minim.mdp -c 2lzm_solv_ions.gro -p topol.top -o em.tpr')
    os.system('gmx_mpi mdrun -deffnm em')

    # Equilibration
    os.system('gmx_mpi grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr -r em.gro')
    os.system('gmx_mpi mdrun -deffnm nvt')

    os.system('gmx_mpi grompp -f npt.mdp -c nvt.gro -p topol.top -o npt.tpr -r nvt.gro')
    os.system('gmx_mpi mdrun -deffnm npt')

def initialize_target():
    os.system('gmx_mpi pdb2gmx -f 150l_C_K162mod.pdb -o 150l_processed.gro -ignh -water spce -ff amber99sb -p 150l_topol.top')

    # Set the simulation config
    os.system('gmx_mpi editconf -f 150l_processed.gro -o 150l_newbox.gro -c -d 1.0 -bt cubic')

    # Add solvate
    os.system('gmx_mpi solvate -cp 150l_newbox.gro -cs spc216.gro -o 150l_solv.gro -p 150l_topol.top')

    # Add ions
    os.system('gmx_mpi grompp -f ions.mdp -c 150l_solv.gro -p 150l_topol.top -o 150l_ions.tpr -maxwarn 5')
    os.system("echo 13 | gmx_mpi genion -s 150l_ions.tpr -o 150l_solv_ions.gro -p 150l_topol.top -pname NA -nname CL -nn 8")

    # Energy minimization
    os.system('gmx_mpi grompp -f minim.mdp -c 150l_solv_ions.gro -p 150l_topol.top -o 150l_em.tpr')
    os.system('gmx_mpi mdrun -deffnm 150l_em')

    # Equilibration
    os.system('gmx_mpi grompp -f nvt.mdp -c 150l_em.gro -p 150l_topol.top -o 150l_nvt.tpr -r 150l_em.gro')
    os.system('gmx_mpi mdrun -deffnm 150l_nvt')

    os.system('gmx_mpi grompp -f npt.mdp -c 150l_nvt.gro -p 150l_topol.top -o 150l_npt.tpr -r 150l_nvt.gro')
    os.system('gmx_mpi mdrun -deffnm 150l_npt')

def first_cycle(n):
    print(n)
    os.system('gmx_mpi grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_%d.tpr -maxwarn 5' % n)
    os.system('gmx_mpi mdrun -deffnm md_0_%d' % n)
    os.system("echo 4 4 | gmx_mpi rms -s 150l_npt.gro -f md_0_%d.trr  -o rmsd_0_%d.xvg -tu ns" % (n, n))

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
    os.system("echo 4 4 | gmx_mpi rms -s 150l_npt.gro -f %s.trr  -o rmsd_%d_%d.xvg -tu ns" % (md, step, n))



def read_rmsd(step, n):
    f = open('rmsd_%d_%d.xvg' % (step, n))
    rmsd = []
    while True:
        l = f.readline()
        if l[0] != '#' and l[0] != '@':
            break
    while l != '':
        rmsd.append(l.split()[1])
        l = f.readline()
    f.close()
    return rmsd


dt = 1 # ps for each step
nsteps = 101 # step for each md including inital state
n_para = 5
MAX_CYCLE = 100

cycle_step = 40
log = np.zeros((MAX_CYCLE - cycle_step + 1, n_para))
# initialize()
# initialize_target()

# 並列処理
# p = Pool(20)
# p.map(first_cycle, range(n_para))
# 逐次処理
# for i in range(n_para):
#     first_cycle(i)


# 途中から再開する場合
pass_flag = True
cycle_step -= 1
log_i = 0

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
        rmsd_list[i * nsteps : (i+1) * nsteps] = read_rmsd(cycle_step, i)
    min_rmsd = rmsd_list.argsort()[:n_para]
    min_rmsd = [(r // nsteps, r % nsteps) for r in min_rmsd]
    log[log_i,:] = [r[0] for r in min_rmsd]
    cycle_step += 1
    log_i += 1


np.savetxt('log_39_100.csv', log, delimiter=',')

