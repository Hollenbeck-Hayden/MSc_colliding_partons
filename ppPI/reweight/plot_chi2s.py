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


path = argv[1]
s_cutoffs = [x for x in argv[2:]]

ns = []
chi2s_original = []
chi2s_weighted = []
for s_cutoff in s_cutoffs:
    cutoff = float(s_cutoff)
    chi2_original = 0
    chi2_weighted = 0
    n = 0
    with open(path + '_' + s_cutoff + '/result_chi2s.txt', 'r') as infile:
        for line in infile:
            items = line.strip().split()
            n += int(items[1])
            chi2_original += float(items[2])
            chi2_weighted += float(items[3])
    chi2s_original.append(chi2_original / float(n))
    chi2s_weighted.append(chi2_weighted / float(n))
    ns.append(n)

cutoffs = [float(x) for x in s_cutoffs]

fig = plt.figure()
gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[3,1], hspace=0)
axes = gs.subplots(sharex=True, sharey=False)

axes[0].plot(cutoffs, chi2s_original, 'o-', label='original')
axes[0].plot(cutoffs, chi2s_weighted, 'o-', label='weighted')

print(ns)

axes[1].plot(cutoffs, ns, 'o-', color='black')
axes[1].set_xlabel('Cutoff [GeV]')
axes[1].set_ylim(0, 140)

axes[0].set_ylabel(r'$\chi^2$/point')
axes[1].set_ylabel(r'Total $N_{dat}$')

axes[0].set_title(r'$\chi^2$ distribution for %s' % path)
axes[0].legend()

#plt.show()
#plt.savefig('%s_chi_dist.png' % path)

setfigsize(fig, 4, 4)
savefig('%s_chi_dist' % path)
