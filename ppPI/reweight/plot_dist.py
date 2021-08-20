from math import sqrt
import sys
import numpy as np
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
    'savefig.pad_inches': 0,
    'savefig.transparent': True,
    'savefig.bbox': 'tight',
});

def savefig(path):
    plt.savefig(path + '.pgf', pad_inches=0)


def setfigsize(fig, width, height):
    fig.set_size_inches(w=width, h=height)


def make_ratio_fig(rows, cols, h_ratios=[3,1]):
    fig = plt.figure()
    gs = fig.add_gridspec(nrows=rows, ncols=cols, height_ratios=h_ratios, hspace=0, wspace=0)
    axes = gs.subplots(sharex=True, sharey=False)
    return fig, axes

def average(data):
    return sum(data) / float(len(data))

def stddev(data):
    a = average(data)
    return sqrt(sum([(x-a)*(x-a) for x in data]) / float(len(data)-1))


def readfile(filename):
    data = []
    with open(filename, 'r') as infile:
        for line in infile:
            data.append(float(line.strip()))
    return data

def plot_hist(axis, xs, bin_count, weights, label, alpha=1, color='green', style='bar'):
    _, bins = np.histogram(xs, bins=bin_count, weights=weights)
    axis.hist(xs, bins=bins, range=(0,8), density=False, weights=weights, color=color, alpha=alpha, histtype=style, label=label)

def calc_Neff(reweights):
    N = float(len(reweights))
    return int(np.exp(sum([w * np.log(N / w) if w > 1e-10 else 0.0 for w in reweights]) / N))

path = sys.argv[1]

chi2_original = readfile(path + "/replica_chi2_original.txt")
weights = readfile(path + "/reweights.data")

original_average = average(chi2_original)
original_stddev = stddev(chi2_original)
chi2_central = [1 if abs(x-original_average) < original_stddev else 0 for x in chi2_original]

print(len(chi2_original))

# ----- CHI 2 DIST ----- #
fig, axes = make_ratio_fig(1, 1, [1])

bin_count = 50
plot_hist(axes, chi2_original, bin_count, weights, "Reweighted", alpha=1, color='orange')
plot_hist(axes, chi2_original, bin_count, [1 for i in range(len(chi2_original))], "Original", alpha=1, color='black', style='step')
axes.legend()

axes.set_xlabel(r'$\chi^2$ distribution')
axes.set_ylabel(r'Count')

setfigsize(fig, 3.3, 2.5)
savefig('%s/chi_replica_dist' % path)
plt.clf()
plt.close(fig)

# ----- WEIGHT DIST ----- #
fig, axes = make_ratio_fig(1, 1, [1])
_, bins = np.histogram(weights, bins=400)
log_bins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
axes.hist(weights, bins=log_bins, density=False, color='green', label='All')

hist, _ = np.histogram(weights, bins=log_bins)
Neff = calc_Neff(weights)
print("Neff: %d" % Neff)

bin_eff = len(hist)-1
for i in range(bin_eff):
    bin_eff = len(hist)-1-i
    Neff -= hist[bin_eff]
    print(Neff)
    if Neff < 0:
        break
print(hist)
print(bin_eff)
print(hist[bin_eff:])
print(log_bins[bin_eff:])
axes.hist(weights, bins=log_bins[bin_eff:], density=False, color='black', histtype='step')
axes.legend()

axes.set_xlabel(r'Weight distribution')
axes.set_xscale('log')
axes.set_xlim(np.power(10.0, -10.0), np.power(10.0, 2.0))
axes.set_ylabel('Count')

setfigsize(fig, 3.3, 2.5)
savefig('%s/weight_dist' % path)
plt.clf()
plt.close(fig)

# ----- PROB DIST ----- #
fig, axes = make_ratio_fig(1, 1, [1])
prob_data = [[], []]
with open(path + '/probability_density.txt', 'r') as infile:
    for line in infile:
        items = [float(x) for x in line.strip().split()]
        prob_data[0].append(items[0])
        prob_data[1].append(items[1])

axes.plot(prob_data[0], prob_data[1], '-')
axes.set_xlabel(r'$\alpha$')
axes.set_ylabel(r'$\mathcal{P}(\alpha)$')

setfigsize(fig, 3.3, 2.5)
savefig(path + '/probability_density')
