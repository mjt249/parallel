import numpy as np
import matplotlib.pyplot as plt
# from scipy import stats

# read files
datadir = "../data/"
list_of_files = [(datadir+"openmp.txt", "OpenMP"),
              (datadir+"mpi.txt", "MPI")]
strongtrials = 7
weakslope = 0
weakconst = 1
strongslope = 0
strongconst = 1
serialt = []
ns = []
ps = []
ts = []
labels = []

# Add a loglog plot and best fit line of xs[i],ys[i] for i in len(xs)
def plotloglog(ax, xs, ys):
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

def makegraph(title, outfile, xlabel, ylabel, xs, ys, idealys):
    fig, ax = plt.subplots()
    plt.xlabel(xlabel)
    ax.set_xscale("log", basex=2)
    plt.ylabel(ylabel)
    ax.set_yscale("log")
    plotloglog(ax, xs, ys)
    plt.title(title)
    print xs[0]
    # print ys
    print idealys
    ax.plot(xs[-1], idealys, "r-", label="Ideal performance, slope = "+
            str(np.polyfit(np.log10(xs[-1]), np.log10(idealys), 1)[0])) # serial T = "+str(serialt[0])+"/p seconds")
    plt.legend()
    plt.savefig(outfile+'.png', bbox_inches='tight')


# read in data files
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
        p = [line.split()[1] for line in lines[starti:]]
        t = [line.split()[2] for line in lines[starti:]]

        n = np.asarray(n, dtype=float)
        p = np.asarray(p, dtype=float)
        t = np.asarray(t, dtype=float)

        ps.append(p)
        ts.append(t)
        labels.append(label)


def prepend(obj, lst):
    out = [obj]
    out.extend(lst)
    return out

serialt = np.asarray(serialt, dtype=float)
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
          "Number of processors",
          "Simulation time (seconds)",
          ps1, ts1, ys1)
makegraph("Weak Scaling: "+str(serialn)+" Particles Per Processor",
          "weak",
          "Number of processors",
          "Efficiency (serial time / parallel time)",
          ps2, ts2, ys2)

