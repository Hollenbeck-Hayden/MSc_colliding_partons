import matplotlib.pyplot as plt
import matplotlib
from math import sqrt

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
    gs = fig.add_gridspec(nrows=rows, ncols=cols, height_ratios=h_ratios, hspace=0, wspace=0.1)
    axes = gs.subplots(sharex=True, sharey=True)
    return fig, axes

def read_bin(infile):
    return [float(x) for x in infile.readline().split()]

def read_file(filename):
    NO_STATE = 0
    READ_BIN = 1
    READ_DATA = 2

    with open(filename, 'r') as infile:
        data = {}
        state = NO_STATE
        target = ''
        header = []

        for line in infile:
            if line.startswith("# x"):
                state = READ_BIN
                target = 'x_bin'
                continue

            if line.startswith("# y"):
                state = READ_BIN
                target = 'y_bin'
                continue

            if line.startswith("# z"):
                state = READ_BIN
                target = 'z_bin'
                continue

            if line.startswith("# data"):
                state = READ_DATA
                target = 'data'
                data[target] = {}
                header = line[6:].split()
                for h in header:
                    data[target][h] = []
                continue

            if state == READ_BIN:
                items = sorted([float(x) for x in line.split()])
                data[target] = [(items[i], items[i+1]) for i in range(len(items)-1)]

            if state == READ_DATA:
                items = [float(x) for x in line.split()]
                for i in range(len(items)):
                    data[target][header[i]].append(items[i])
    return data

def in_bin(x, bounds):
    return x >= bounds[0] and x <= bounds[1]

def select(data, x_bin, y_bin):

    data_bin = {}
    data_bin['z'] = []
    data_bin['M'] = []
    data_bin['Mt'] = []
    data_bin['Merr'] = []
    for i in range(len(data['data']['z'])):
        if in_bin(data['data']['x'][i], x_bin) and in_bin(data['data']['y'][i], y_bin):
            data_bin['z'].append(data['data']['z'][i])
            data_bin['M'].append(data['data']['M'][i])
            data_bin['Mt'].append(data['data']['Mt'][i])
            if 'Merr' in data['data'].keys():
                data_bin['Merr'].append(data['data']['Merr'][i])

    data_bin['dz'] = [(data_bin['z'][i] + data_bin['z'][i+1])/2 for i in range(len(data_bin['z'])-1)]
    data_bin['dM'] = [(data_bin['Mt'][i+1] - data_bin['Mt'][i])/(data_bin['z'][i+1]-data_bin['z'][i]) for i in range(len(data_bin['z'])-1)]

    #print('x: [%5.3f, %5.3f]\ty: [%5.3f, %5.3f]' % (x_bin[0], x_bin[1], y_bin[0], y_bin[1]))
    #print('\tz: ' + str(data_bin['z']))
    #print('\tM: ' + str(data_bin['M']))
    #print('\tMt: ' + str(data_bin['Mt']))
    return data_bin

def arr_sum(a, b):
    return [a[i] + b[i] for i in range(len(a))]

def arr_diff(a, b):
    return [a[i] - b[i] for i in range(len(a))]



data = read_file("results0.data")

print(data.keys())
print(data['data'].keys())

bins = [(data['x_bin'][i], data['y_bin'][3]) for i in range(0, 2)]

data_avg = [[0 for x in select(data, b[0], b[1])['z']] for b in bins]
data_mmt = [[0 for x in select(data, b[0], b[1])['z']] for b in bins]

for r in range(1,201):
    replica_data = read_file("results" + str(r) + ".data")
    for i in range(2):
        replica_binned = select(replica_data, bins[i][0], bins[i][1])
        for j in range(len(replica_binned['Mt'])):
            data_avg[i][j] += replica_binned['Mt'][j]                           / 200.0
            data_mmt[i][j] += replica_binned['Mt'][j] * replica_binned['Mt'][j] / 200.0

data_var = [[sqrt(data_mmt[i][j] - (data_avg[i][j] * data_avg[i][j])) for j in range(len(data_avg[i]))] for i in range(2)]

fig, axes = make_ratio_fig(1, 2, h_ratios=[1])

axes[0].set_ylabel(r'$M$')
for i in range(2):
    data_bin = select(data, bins[i][0], bins[i][1])

    print('bins: ' + str(bins[i]))
    print(data_bin['Merr'])

    axes[i].set_xlim(0.3, 0.7)
    axes[i].set_ylim(0, 1)

    axes[i].set_xlabel(r'$z$')

    axes[i].plot(data_bin['z'], data_avg[i], '-', label='MAPFF')
    axes[i].fill_between(data_bin['z'], arr_sum(data_avg[i], data_var[i]), arr_diff(data_avg[i], data_var[i]), alpha=0.6)
    axes[i].errorbar(data_bin['z'], data_bin['M'], yerr=data_bin['Merr'], fmt='.', label='COMPASS')
    axes[i].legend()

setfigsize(fig, 6, 3)
savefig('compass')
