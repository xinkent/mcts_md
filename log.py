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


o = open('min_rmsd.txt','w')
for i in range(100):
    rmsd = read_rmsd(i,0)
    o.write(rmsd[0] + '\n')
o.close()
