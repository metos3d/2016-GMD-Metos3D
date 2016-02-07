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
#   create_figures
#
def create_figures():
    import sys
    import numpy as np
    import math as m
    # create figures
    print "# Creating figures ..."
    # ew params
    gammas = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5]
    gammas = np.array(gammas[::-1])
    alphas = [(1.0 + m.sqrt(5.0))/2.0, 1.5, 1.4, 1.3, 1.2, 1.1]
    alphas = np.array(alphas[::-1])
    ngamma = len(gammas)
    nalpha = len(alphas)
    # data to plot
    ksp2d = np.zeros(shape=(ngamma, nalpha))
    snes2d = np.zeros(shape=(ngamma, nalpha))
    # files
    sneskspdata = []
    for gidx in range(0, ngamma):
        for aidx in range(0, nalpha):
            astr = "{0:.5f}".format(alphas[aidx])
            gstr = "{0:1.1f}".format(gammas[gidx])
            file = "../work/NK-" + gstr + "-" + astr + ".out.txt"
#            print file
            snesidx, sneslist, kspiter, ksplist = read_snes_ksp_data(file)
            ksp2d[gidx, aidx] = kspiter
            snes2d[gidx, aidx] = len(snesidx)
#            print gstr, astr, len(snesidx), kspiter
    # plot
    create_plot(alphas, gammas, ksp2d, snes2d)

#
#   create_plot
#
def create_plot(alphas, gammas, ksp2d, snes2d):
    import numpy as np
    import matplotlib
    matplotlib.rc("font", **{"family" : "sans-serif"})
    matplotlib.rcParams.update({'font.size': 14})
    matplotlib.rc("text", usetex = True)
    matplotlib.use("PDF")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    # grid
    X, Y = np.meshgrid(alphas, gammas)
    ngamma = len(gammas)
    nalpha = len(alphas)
    # view
    elev = 30
    azi = 65
    
    # KSP
    fig = plt.figure(1)
    ax = Axes3D(fig)
    ax.plot_surface(X, Y, ksp2d, lw = 0.1, rstride = 1, cstride = 1, cmap = plt.cm.bone_r, alpha = 0.5)
    ax.plot_wireframe(X, Y, ksp2d, color = "k", linewidth = 1.0)
    # axes
    ax.set_xlim([1.01, 1.69])
    ax.set_ylim([0.41, 1.09])
    ax.set_zlim([401, 699])
    # labels
    ax.set_xlabel(r"Alpha ($\alpha$)")
    ax.set_ylabel(r"Gamma ($\gamma$)")
    ax.set_zlabel("Model years")
    ax.set_xticklabels(["%.1f" % alpha for alpha in alphas])
    ax.set_yticklabels(["%.1f" % gamma for gamma in gammas])
    ax.set_zticklabels(["450", "500", "550", "600", "650"])
    # view
    ax.view_init(elev, azi)
    # save
    from matplotlib.transforms import Bbox
    b = Bbox.from_bounds(0.25, 0, 7.25, 5.75)
    plt.savefig("eisenstat-walker-ksp.pdf", bbox_inches = b)

    # SNES
    fig = plt.figure(2)
    ax = Axes3D(fig)
    # plot
    ax.plot_surface(X, Y, snes2d, lw = 0.1, rstride = 1, cstride = 1, cmap = plt.cm.bone_r, alpha = 0.5)
    ax.plot_wireframe(X, Y, snes2d, color = "k", linewidth = 1.0)
    # axes
    ax.set_xlim([1.01, 1.69])
    ax.set_ylim([0.41, 1.09])
    ax.set_zlim([8.01, 21.9])
    # labels
    ax.set_xlabel(r"Alpha ($\alpha$)")
    ax.set_ylabel(r"Gamma ($\gamma$)")
    ax.set_zlabel(r"Newton steps")
    ax.set_xticklabels(["%.1f" % alpha for alpha in alphas])
    ax.set_yticklabels(["%.1f" % gamma for gamma in gammas])
    ax.set_zticklabels(["10", "12", "14", "16", "18", "20"])
    # view
    ax.view_init(elev, azi)
    # save
    from matplotlib.transforms import Bbox
    b = Bbox.from_bounds(0.25, 0, 7.25, 5.75)
    plt.savefig("eisenstat-walker-snes.pdf", bbox_inches = b)

#
#   main
#
if __name__ == "__main__":
    create_figures()




