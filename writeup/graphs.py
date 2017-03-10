import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
# from scipy import stats

# read files
datadir = "../data/"
complexityfiles = [(datadir+"serial.txt", "Serial Binning"),
              (datadir+"barneshut.txt", "Barnes-Hut")]
scalingfiles = [
                (datadir+"openmp.txt", "OpenMP - 50000 particles"),
                (datadir+"openmp5000.txt", "OpenMP - 5000 particles"),
                (datadir+"mpi.txt", "MPI - 50000 particles"),
                (datadir+"mpi5000.txt", "MPI - 5000 particles"),
                (datadir+"serial.txt", "Serial Binning"),
                (datadir+"barneshut.txt", "Serial Barnes-Hut")]
strongtrials = 7

# Add a loglog plot and best fit line of xs[i],ys[i] for i in len(xs)
def plotloglog(ax, xs, ys, ideals):
    colors = "gbmkcy"
    shapes = "^voxsd"
    regression = True
    for i in range(len(xs)):
        xs[i] = np.asarray(xs[i], dtype=float)
        ys[i] = np.asarray(ys[i], dtype=float)
        logx = np.log10(xs[i])
        logy = np.log10(ys[i])
        ax.plot(xs[i], ys[i], colors[i]+shapes[i], label=labels[i])
        regression = (xs[i][0] != xs[i][-1])
        if regression:
            [slope,const] = np.polyfit(logx, logy, 1)
            const = np.power(10,const)
            ax.plot(xs[i], const * (xs[i] ** slope), colors[i]+"-",
                    label="slope = "+str(slope))
    ax.plot(xs[0], ideals, "r-", label="Ideal performance, slope = "+
            str(np.polyfit(np.log10(xs[0]), np.log10(ideals), 1)[0]))

def makegraph(title, outfile, base, xlabel, ylabel, xs, ys, idealys, plotfunc):
    fig, ax = plt.subplots()
    plt.margins(x=0.1, y=0.1)
    plt.xlabel(xlabel)
    ax.set_xscale("log", basex=base)
    plt.ylabel(ylabel)
    ax.set_yscale("log")
    plotfunc(ax, xs, ys, idealys)
    plt.title(title)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
    plt.savefig(outfile+'.png', bbox_inches='tight')


# read in data files
def read_input_files(list_of_files):
    serialt = []
    ns = []
    ps = []
    ts = []
    labels = []
    for filename, label in list_of_files:
        with open(filename) as f:
            lines = f.readlines()
            # parse serial line first
            if len(lines[0].split()) == 2:
                serialn = lines[0].split()[0]
                serialt.append(lines[0].split()[1])
                starti = 1
            else:
                serialn = lines[0].split()[0]
                serialt.append(lines[0].split()[2])
                starti = 0

            n = [line.split()[0] for line in lines[starti:]]
            if len(lines[1].split()) == 2: # it's a serial file
                p = [1 for line in lines]
                n = [line.split()[0] for line in lines]
                t = [line.split()[1] for line in lines]
            else:
                n = [line.split()[0] for line in lines[starti:]]
                p = [line.split()[1] for line in lines[starti:]]
                t = [line.split()[2] for line in lines[starti:]]

            print "n for", label, ":",n
            ns.append(n)
            ps.append(p)
            ts.append(t)
            labels.append(label)
    serialt = np.asarray(serialt, dtype=float)
    return serialn, serialt, ns, ps, ts, labels

def prepend(obj, lst):
    out = [obj]
    out.extend(deepcopy(lst))
    return out

# serial complexity
serialn, serialt, ns, ps, ts, labels = read_input_files(complexityfiles)
makegraph("Serial Performance", 'serial',
          10,
          "Number of particles",
          "Simulation time (seconds)",
          # [ns[i][:7] for i in range(len(ns))],
          # [ts[i][:7] for i in range(len(ts))],
          # (ns[-1]/250)[:7],
          ns,
          ts,
          [float(i)/1500 for i in ns[0]],
          plotloglog)
serialn, serialt, ns, ps, ts, labels = read_input_files(scalingfiles)
# strong
ps1 = [ps[i][:strongtrials] for i in range(len(ps) - 2)]
ts1 = [ts[i][:strongtrials] for i in range(len(ts) - 2)]
# # exclude outliers from best fit curve
# ps1[1] = ps1[1][:-2]
# ts1[1] = ts1[1][:-2]
ys1 = np.mean(serialt)/np.asarray(ps1[0], dtype=float)
makegraph("Strong Scaling: 5000 Particles",
          "strong",
          2,
          "Number of processors",
          "Simulation time (seconds)",
          ps1, ts1, ys1, plotloglog)
# weak
ns2 = [prepend(ns[i][0],ns[i][strongtrials:]) for i in range(len(ns) - 2)]
ns2.append(ns[-2][4:])
ns2.append(ns[-1][3:])
ps2 = [prepend(ps[i][0],ps[i][strongtrials:]) for i in range(len(ps) - 2)]
ps2.append(ps[-2][4:])
ps2.append(ps[-1][3:])
ts2 = [serialt[i]/np.asarray(prepend(ts[i][0],ts[i][strongtrials:]),
                                        dtype=float) for i in
                  range(len(ts)-2)]
ts2.append(serialt[i]/np.asarray(ts[-2], dtype=float)[4:])
ts2.append(serialt[i]/np.asarray(ts[-1], dtype=float)[3:])
ys2 = np.asarray([1.0]*len(ns2[0]))
print "ns2[0]",ns2[0],"ys2",ys2
makegraph("Weak Scaling",
          "weak",
          2,
          "Number of particles",
          "Efficiency (serial time / parallel time)",
          ns2, ts2, ys2, plotloglog)
