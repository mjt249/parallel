import numpy as np
import matplotlib.pyplot as plt
# from scipy import stats

# read files
datadir = "../data/"
complexityfiles = [(datadir+"serial.txt", "Serial Binning"),
              (datadir+"barneshut.txt", "Barnes-Hut")]
scalingfiles = [(datadir+"openmp.txt", "OpenMP"),
              (datadir+"mpi.txt", "MPI")]
strongtrials = 7

# Add a loglog plot and best fit line of xs[i],ys[i] for i in len(xs)
def plotloglog(ax, xs, ys, ideals):
    colors = "gb"
    shapes = "ov"
    for i in range(len(xs)):
        # print i
        # print xs[i]
        logx = np.log10(xs[i])
        logy = np.log10(ys[i])
        [slope,const] = np.polyfit(logx, logy, 1)
        const = np.power(10,const)
        ax.plot(xs[i], ys[i], colors[i]+shapes[i], label=labels[i])
        ax.plot(xs[i], const * (xs[i] ** slope), colors[i]+"-",
                label="slope = "+str(slope))
    ax.plot(xs[-1], ideals, "r-", label="Ideal performance, slope = "+
            str(np.polyfit(np.log10(xs[-1]), np.log10(ideals), 1)[0]))

def makegraph(title, outfile, base, xlabel, ylabel, xs, ys, idealys, plotfunc):
    fig, ax = plt.subplots()
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
                p = [1 for line in lines[starti:]]
                t = [line.split()[1] for line in lines[starti:]]
            else:
                p = [line.split()[1] for line in lines[starti:]]
                t = [line.split()[2] for line in lines[starti:]]

            n = np.asarray(n, dtype=float)
            p = np.asarray(p, dtype=float)
            t = np.asarray(t, dtype=float)

            ns.append(n)
            ps.append(p)
            ts.append(t)
            labels.append(label)
    serialt = np.asarray(serialt, dtype=float)
    return serialn, serialt, ns, ps, ts, labels

def prepend(obj, lst):
    out = [obj]
    out.extend(lst)
    return out

serialn, serialt, ns, ps, ts, labels = read_input_files(complexityfiles)
makegraph("Serial Performance", 'serial',
          10,
          "Number of particles",
          "Simulation time (seconds)",
          ns, ts, ns[-1]/250, plotloglog)
serialn, serialt, ns, ps, ts, labels = read_input_files(scalingfiles)
# strong
ps1 = [ps[i][:strongtrials] for i in range(len(ps))]
ts1 = [ts[i][:strongtrials] for i in range(len(ts))]
ys1 = np.mean(serialt)/ps1[0]
# weak
ps2 = [prepend(ps[i][0],ps[i][strongtrials:]) for i in range(len(ps))]
ts2 = np.asarray([serialt[i]/prepend(ts[i][0],ts[i][strongtrials:]) for i in
                  range(len(ts))])
ys2 = np.asarray([1.0]*len(ps2[i]))
makegraph("Strong Scaling: "+str(serialn)+" Particles",
          "strong",
          2,
          "Number of processors",
          "Simulation time (seconds)",
          ps1, ts1, ys1, plotloglog)
makegraph("Weak Scaling: "+str(serialn)+" Particles Per Processor",
          "weak",
          2,
          "Number of processors",
          "Efficiency (serial time / parallel time)",
          ps2, ts2, ys2, plotloglog)

