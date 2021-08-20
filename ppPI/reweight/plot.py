from sys import argv
from math import sqrt,log10
import matplotlib.pyplot as plt
import matplotlib

TEXT_WIDTH = 6      # about 6 in text width


matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    'axes.labelsize': 'small',
    'xtick.labelsize': 'small',
    'ytick.labelsize': 'small',
    'legend.fontsize': 'small',
    'legend.handlelength': 0.5,
    'legend.borderaxespad': 0,
    'legend.frameon': False,
    'savefig.pad_inches': 0.0,
    'savefig.transparent': True,
    'savefig.bbox': 'tight',
});

def savefig(path):
    plt.savefig(path + '.pgf')


def setfigsize(fig, width, height):
    fig.set_size_inches(w=width, h=height)


def make_ratio_fig(rows, cols, h_ratios=[3,1]):
    fig = plt.figure()
    gs = fig.add_gridspec(nrows=rows, ncols=cols, height_ratios=h_ratios, hspace=0, wspace=0)
    axes = gs.subplots(sharex=True, sharey=False)
    return fig, axes


DATA_PATH = "../data"
THEORY_PATH = "../theory"
EXPERIMENTS_PATH = "../reweight/plot_experiments.txt"

cutoff = 5

def arr_sum(a, b):
    return [a[i] + b[i] for i in range(len(a))]

def arr_dif(a, b):
    return [a[i] - b[i] for i in range(len(a))]

def arr_ratio(a, b):
    return [a[i] / b[i] for i in range(len(a))]

def get_theory_path(experiment, prefix, replica):
    formatted_replica = "%03d" % replica
    return "../theory/" + experiment + "/results/" + prefix + formatted_replica + ".res"

def get_data_path(experiment):
    return "../data/" + experiment + "/" + experiment + ".txt"

def read_datafile(experiment):
    data_x = []
    data_y = []
    data_err = []
    i = 0
    with open(get_data_path(experiment), 'r') as infile:
        for line in infile:
            if i < 5:
                i += 1
            else:
                items = [float(x) for x in line.split()]
                if items[1] < cutoff: continue
                data_x.append(items[1])
                data_y.append(items[2])
                data_err.append(sqrt(items[3]*items[3] + items[4]*items[4]))
    return data_x, data_y, data_err
 
def compute_avg_and_error(ys, weights=[]):
    if (len(weights) == 0):
        weights = [1. for y in ys]

    #for i in range(len(weights)):
    #    for y in ys[i]:
    #        if y < 0:
    #            weights[i] = 0
    #            print("Weight %d set to 0" % i)

    # Compute average and 2nd moment
    theory_y_avg = [0 for i in range(len(ys[0]))]
    theory_y_2mt = [0 for i in range(len(ys[0]))]
    theory_y_err = [0 for i in range(len(ys[0]))]
    for i in range(len(theory_y_avg)):
        N = 0
        for j in range(len(ys)):
            y = ys[j][i]
            weight = weights[j]
            # If the theory point is negative / NaN, don't include it
            if y < 0:
                weight = 0
                print("y[%d][%d] = %f" % (j, i, y))
            theory_y_avg[i] += weight * y
            theory_y_2mt[i] += weight * y * y
            N += weight
        theory_y_avg[i] /= N
        theory_y_2mt[i] /= N
        # compute std dev from average and second moment
        theory_y_err[i] = sqrt(theory_y_2mt[i] - theory_y_avg[i]*theory_y_avg[i])

    return theory_y_avg, theory_y_err


def read_theoryfile(experiment, prefix, replicas):
    theory_x = []
    theory_y = [[] for i in range(1, replicas+1)]
    for i in range(1, replicas+1):
        with open(get_theory_path(experiment, prefix, i), 'r') as infile:
            for line in infile:
                items = [float(x) for x in line.split()]
                if items[0] < cutoff: continue
                if (i == 1):
                    theory_x.append(items[0])
                theory_y[i-1].append(items[1])

    return theory_x, theory_y

def read_centraltheory(experiment, prefix):
    theory_x = []
    theory_y = []
    with open(get_theory_path(experiment, prefix, 0), 'r') as infile:
        for line in infile:
            items = [float(x) for x in line.split()]
            if items[0] < cutoff: continue
            theory_x.append(items[0])
            theory_y.append(items[1])

    return theory_x, theory_y

def read_unweightfile(filename):
    data = []
    with open(filename, 'r') as infile:
        for line in infile:
            data.append(float(line))
    return data

def draw_fill_between(axis, x, y, err, label):
    axis.plot(x, y, label=label)
    axis.fill_between(x, arr_sum(y, err), arr_dif(y, err), alpha=0.6)


def read_post_process(experiment):
    obs = ""
    conv = 1
    name = ''
    exp = ''
    with open("../theory/" + experiment + "/post_process.txt", 'r') as infile:
        for line in infile:
            if line.startswith("obs: "):
                obs = line.strip()[5:]
            elif line.startswith("conv: "):
                conv = float(line.strip().split()[1])
            elif line.startswith("name: "):
                name = line.strip()[6:]
            elif line.startswith("exp: "):
                exp = line.strip()[5:]

    return (obs, conv, name, exp)

def read_chi2(path, experiment):
    with open(path + "/result_chi2s.txt", 'r') as infile:
        for line in infile:
            if line.startswith(experiment):
                items = [float(x) for x in line.strip().split()[2:]]
                return (items[0], items[1])
    return (-1, -1)

def plot(experiment, prefix, path, replicas):
    data_x, data_y, data_y_err = read_datafile(experiment)                  # Experimental data
    #theory_x, theory_y = read_centraltheory(experiment, prefix)            # Central data
    theory_x, theory_ys = read_theoryfile(experiment, prefix, replicas)     # Theory w/ error data
    unweights = read_unweightfile(path + "/reweights.data")                         # Reweight data
    observable, conversion, name, _ = read_post_process(experiment)                  # Post-process data
    original_chi2, weighted_chi2 = read_chi2(path, experiment)

    # Apply unit conversion
    for i in range(len(theory_ys)):
        for j in range(len(theory_ys[i])):
            theory_ys[i][j] *= conversion

    # Compute central and average from replicas
    theory_y, theory_y_err = compute_avg_and_error(theory_ys)

    if (len(theory_x) == 0):
        print("No theory predictions found above the cutoff!!")
        return

    print("Num unweights: " + str(len(unweights)))
    print("Num theory_ys: " + str(len(theory_ys)))

    #print("Unweights: " + str(unweights))

    weighted_y, weighted_err = compute_avg_and_error(theory_ys, unweights)

    #for i in range(len(data_x)):
    #    print("data  [%6.4e] = %6.4e" % (data_x[i], data_y[i]))
    #    print("theory[%6.4e] = %6.4e" % (theory_x[i], theory_y[i]))
    #    print("weight[%6.4e] = %6.4e" % (theory_x[i], weighted_y[i]))

    fig, axes = make_ratio_fig(2, 1)
    axes[0].set_yscale('log')

    epsilon = 0.1
    #theory_xmid = [(theory_x[i+1] + theory_x[i])/2.0 for i in range(len(theory_x)-1)]      # Bar midpoint
    #theory_xlow = [theory_x[0] - (theory_x[1] - theory_x[0])/2.0] + theory_xmid            # Bar low
    #theory_xhigh = theory_xmid + [theory_x[-1] + (theory_x[-1]-theory_x[-2])/2.0]          # Bar high
    #axes[0].bar(theory_xlow, theory_y_err, arr_dif(theory_xhigh, theory_xlow), arr_dif(theory_y, theory_y_err), align='edge', color='red', alpha=0.6, label="theory")              # Bar plots
    draw_fill_between(axes[0], theory_x, theory_y, theory_y_err, "Original")              # Fill between theory w/ error
    #axes[0].plot(theory_x, theory_y, '-', label="theory")                              # Plot central
    draw_fill_between(axes[0], theory_x, weighted_y, weighted_err, "Reweighted")       # Fill between reweighted
    axes[0].errorbar(data_x, data_y, yerr=data_y_err, fmt='.', label=name)      # Plot experiment
    #axes[0].set_xlabel(r'$p_T$ [GeV]')
    axes[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    axes[0].set_ylabel(r'$' + observable + '$')
    #axes[0].set_xlim(0.1, 40)
    #axes[0].set_xlim(min(theory_xlow)-epsilon, max(theory_xhigh)+epsilon)
    #a = int(log10(min(theory_y)))-2
    #b = int(log10(max(theory_y)))+1
    #axes[0].set_ylim(10**a, 10**b)

    norm_data = arr_ratio(data_y, data_y)                       # Experimental ratio
    norm_data_err = arr_ratio(data_y_err, data_y)               # Experimental error ratio
    norm_theory = arr_ratio(theory_y, data_y)                   # Theory ratio
    norm_theory_err = arr_ratio(theory_y_err, data_y)           # Theory error ratio
    norm_weighted = arr_ratio(weighted_y, data_y)              # Reweighted ratio
    norm_weighted_err = arr_ratio(weighted_err, data_y)        # Reweighted error ratio

    #axes[1].bar(theory_xlow, norm_theory_err, arr_dif(theory_xhigh, theory_xlow), arr_dif(norm_theory, norm_theory_err), align='edge', color='red', alpha=0.6, label="theory")     # Bar plots
    draw_fill_between(axes[1], theory_x, norm_theory, norm_theory_err, "Original")                # Fill between theory w/ error
    draw_fill_between(axes[1], theory_x, norm_weighted, norm_weighted_err, "Reweighted")       # Fill between reweighted w/ error
    #axes[1].plot(theory_x, norm_theory, '-', label="theory")                                   # Plot central
    axes[1].errorbar(data_x, norm_data, yerr=norm_data_err, fmt='.', label=name)        # Plot experiment
    axes[1].set_xlabel(r'$p_T$ [GeV]')
    axes[1].set_ylabel("Ratio")
    #axes[1].set_xlim(0.1, 40)
    #axes[1].set_xlim(min(theory_xlow)-epsilon, max(theory_xhigh)+epsilon)                      # Enforce theory x-limit

    axes[0].legend()


    #setfigsize(fig, 3, 2.5)
    if experiment.startswith("ALICE"):
        setfigsize(fig, 5, 2.5)
    else:
        setfigsize(fig, 3, 2.5)
    savefig('%s/%s_plot' % (path, experiment))
    #plt.savefig(path + "/" + experiment + "_plot.png")
    #plt.show()
    plt.clf()
    plt.close(fig)


def plot_chi2s(path, exp_names):
    ns = []
    chi2s = [[], []]
    labels = []
    with open(path + '/result_chi2s.txt') as infile:
        for line in infile:
            items = line.strip().split()
            name = items[0]
            if name in exp_names:
                ns.append(int(items[1]))
                chi2s[0].append(float(items[2]))
                chi2s[1].append(float(items[3]))
                labels.append(name)

    for i in range(len(labels)):
        _, _, name, exp = read_post_process(labels[i])
        labels[i] = '%s %s' % (exp, name)

    ns.append(sum(ns))
    print("total points: %d" % ns[-1])
    for i in range(len(chi2s)):
        chi2s[i].append(sum(chi2s[i]))
        for j in range(len(chi2s[i])):
            chi2s[i][j] /= float(ns[j])
    labels.append('Total')

    width = 0.2
    indexes = [i for i in range(len(chi2s[0]))]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar([x         for x in indexes], chi2s[0], width, label='Original')
    ax.bar([x + width for x in indexes], chi2s[1], width, label='Reweighted')

    ax.set_ylabel(r'$\chi^2$/point')
    ax.set_xticks(indexes)
    ax.set_xticklabels(labels, rotation=30, ha='right')
    ax.legend()

    setfigsize(fig, 6, 2.5)
    #plt.savefig(path + '/chi2_dist.png')
    #savefig('%s/chi2_dist' % path)
    plt.show()

def main():
    for path in argv[1:]:
        exps = []
        with open(path + "/plot_experiments.txt") as experiment_file:
            i = 0
            N_replica = 0
            cutoff = 0
            for line in experiment_file:
                if line.startswith("#"):
                    continue
                elif len(line.strip()) == 0:
                    continue
                elif i == 0:
                    N_replica = int(line)
                    i += 1
                    continue
                elif i == 1:
                    cutoff = float(line)
                    i += 1
                    continue
                else:
                    items = line.strip().split()
                    experiment_name = items[0]
                    experiment_prefix = items[1]
                    exps.append(experiment_name)
                    print("Experiment: " + experiment_name)
                    plot(experiment_name, experiment_prefix, path, N_replica)
        plot_chi2s(path, exps)


if __name__ == '__main__':
    main()
