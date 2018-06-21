import numpy as np
import os

log = np.loadtxt("log_100.csv", delimiter = ",")
back_list = []
step = log.shape[0] - 1
log_index = 0

for i in range(step,0, -1):
    log_index = int(log[i, log_index])
    back_list.insert(0, log_index)

print(back_list)

times = []
trajs = []
for i, j in enumerate(back_list):
    traj = "res_%d_%d"%(i,j)
    print(traj)
    trajs.append(traj)
    """
    os.system("gmx_mpi check -f res_%d_%d.trr >& check.txt"%(i,j))
    check = open('check.txt')
    while True:
        l = check.readline()
        if l[:4] == "Time":
            time = int(l.split()[1]) - 1
            break
    check.close()
    times.append(time)
    """

traj_str = ""
for traj in trajs:
    traj_str += (traj + " ")
# traj_str += "res_" + str(step) + "_0.trr "

"""
time_str = ""
time_accum = 0
for time in times:
    time_str += (str(time_accum) + " ")
    time_accum += time
time_str += (str(time_accum) + " ")
"""

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

print(traj_str)
os.system( "gmx_mpi trjcat -f " + traj_str + " -o merged_tmp.trr ")
modify_rmsd('merged_tmp.trr', 'merged.trr')
os.remove('merged_tmp.trr')
os.system("gmx_mpi rms -s 150l_npt.gro -f merged.trr -tu ns -o rmsd_merged.xvg")
