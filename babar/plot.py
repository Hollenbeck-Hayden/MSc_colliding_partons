import matplotlib.pyplot as plt
from math import log
from matplotlib import gridspec
import matplotlib

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
    'axes.titlesize': 'medium',
});

def savefig(path):
    plt.savefig(path + '.pgf', pad_inches=0)
    #plt.savefig(path + '.png')


def setfigsize(fig, width, height):
    fig.set_size_inches(w=width, h=height)
    fig.tight_layout()


def make_ratio_fig(rows, cols, h_ratios=[3,1]):
    fig = plt.figure()
    gs = fig.add_gridspec(nrows=rows, ncols=cols, height_ratios=h_ratios, hspace=0, wspace=0)
    axes = gs.subplots(sharex=True, sharey=False)
    return fig, axes


z_min = 0.075
z_max = 0.9

def arr_sum(a, b):
    return [a[i] + b[i] for i in range(len(a))]

def arr_dif(a, b):
    return [a[i] - b[i] for i in range(len(a))]

def accept_z(z):
    return z > z_min and z < z_max

def z_cut(arr, zs):
    return [arr[i] for i in range(len(arr)) if accept_z(zs[i])]

def unwrap_data(data):
    midpoint = [(data[0][i] + data[1][i]) / 2 for i in range(len(data[0]))]
    return (midpoint, data[2], data[3], data[0], data[1], data[4], data[5])
#    return (data[0], data[1], data[2],
#            z_cut(data[0], data[0]),
#            z_cut(data[3], data[0]),
#            z_cut(data[4], data[0]))

def readfile(filename):
    data = [[], [], [], [], [], []]
    with open(filename, "r") as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            items = line.split()
            for i in range(len(data)):
                data[i].append(float(items[i]))
    return data


data = readfile("table2.data")

ex_zs, ex_ys, ex_er, th_zl, th_zh, th_ys, th_er = unwrap_data(data)

fig, axes = make_ratio_fig(2, 1)

axes[0].set_yscale("log")
#axes[0].set_xlabel("z")
axes[0].set_xlim(0, 1)
axes[0].set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
axes[0].set_ylabel(r"$\frac{1}{\sigma_{tot}} \frac{d\sigma}{dz}$")
axes[0].set_ylim(0.001, 10)
#axes[0].set_title("BaBar Differential Cross-Section for the Charged Pion")

axes[0].fill_between(ex_zs, arr_sum(th_ys, th_er), arr_dif(th_ys, th_er), alpha=0.8, color='red', label='MAPFF')
axes[0].errorbar(ex_zs, ex_ys, yerr=ex_er, fmt='.', label='BABAR')

axes[1].set_xlabel("z")
axes[1].set_xlim(0, 1)
axes[1].set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
axes[1].set_ylim(0.95, 1.05)
axes[1].set_yticks([0.95, 1.0])
axes[1].set_ylabel("Ratio")

for i in range(len(data[0])):
    data[5][i] /= data[2][i]
    data[4][i] /= data[2][i]
    data[3][i] /= data[2][i]
    data[2][i] /= data[2][i]

ex_zs, ex_yz, ex_er, th_zl, th_zh, th_ys, th_er = unwrap_data(data)

#axes[1].fill_between(th_zs, arr_sum(th_ys, th_er), arr_dif(th_ys, th_er), alpha=0.5, color='red')
axes[1].bar(th_zl, arr_sum(th_er, th_er), arr_dif(th_zh, th_zl), arr_dif(th_ys, th_er), align='edge', color='red', alpha=0.8, label='MAPFF')
axes[1].errorbar(ex_zs, ex_ys, yerr=ex_er, fmt='.', label='BABAR')
#axes[1].errorbar(th_zs, th_ys, yerr=th_er, fmt='-')

axes[0].legend()

setfigsize(fig, 5, 2.5)
savefig('babar')
