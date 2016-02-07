#!/usr/bin/env python

#
#   read_snes_ksp_data
#
def read_snes_ksp_data(file):
    import re
    # max readin
    nmax = 10000
    # read file, parse and store
    kspiter = 0
    ksplist = []
    snesidx = []
    sneslist = []
    f = open(file, "r")
    for line in f:
        # SNES
        snesmatch = re.search(r"^\s+(\d+) SNES Function norm (\d+.\d+e[+-]\d+)", line)
        if snesmatch:
            snesidx.append(kspiter)
            sneslist.append(snesmatch.groups()[1])
        # KSP
        kspmatch = re.search(r"^\s+(\d+) KSP Residual norm (\d+.\d+e[+-]\d+)", line)
        if kspmatch:
            ksplist.append(kspmatch.groups()[1])
            kspiter = kspiter + 1
            if kspiter >= nmax:
                break
    f.close()
    return snesidx, sneslist, kspiter, ksplist

#
#   read_data
#
def read_data():
    import numpy as np
    print "# Reading data ..."
    # read all samples
    data  = []
    nsample = 100
    for isample in range(0, nsample):
        filename = "../work/metos3d.nesh-fe.sample.%02d.out.txt" % isample
#        print filename
        sampledata = read_snes_ksp_data(filename)
        # number of Newton steps, number of model years
        data.append((len(sampledata[0]), sampledata[2]))
    return np.array(data)

#
#   create_figures
#
def create_figures(data):
    import numpy as np
    print "# Creating figure ..."
    # prepare matplotlib
    import matplotlib
    matplotlib.rc("font",**{"family":"sans-serif"})
    matplotlib.rcParams.update({'font.size': 14})
    matplotlib.rc("text", usetex=True)
    matplotlib.use("PDF")
    import matplotlib.pyplot as plt

    # KSP
    plt.figure(1)
    n, bins, patches = plt.hist(data[:, 1], bins = 50, fc = "k", ec = "w")
    plt.xticks(range(300, 1001, 100), range(300, 1001, 100))
    plt.yticks(range(0, 17, 2), range(0, 17, 2))
    plt.xlabel("Model years")
    plt.ylabel("Occurrence")
    plt.savefig("parameterrange-ksp", bbox_inches = "tight")

    # SNES
    plt.figure(2)
    n, bins, patches = plt.hist(data[:, 0], bins = 50, fc = "k", ec = "w")
    plt.xticks(range(5, 46, 5), range(5, 46, 5))
    plt.yticks(range(0, 15, 2), range(0, 15, 2))
    plt.xlabel("Newton steps")
    plt.ylabel("Occurrence")
    plt.savefig("parameterrange-snes", bbox_inches = "tight")

#
#   main
#
if __name__ == "__main__":
    data = read_data()
    create_figures(data)


