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


