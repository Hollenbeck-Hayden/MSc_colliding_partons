import matplotlib.pyplot as plt
import matplotlib

#matplotlib.use("pgf")
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
    #plt.savefig(path + '.pgf', pad_inches=0)
    plt.savefig(path + '.png')


def setfigsize(fig, width, height):
    fig.set_size_inches(w=width, h=height)
    fig.tight_layout()


def make_ratio_fig(rows, cols, h_ratios=[3,1]):
    fig = plt.figure()
    gs = fig.add_gridspec(nrows=rows, ncols=cols, height_ratios=h_ratios, hspace=0, wspace=0)
    axes = gs.subplots(sharex=True, sharey=False)
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
    for i in range(len(data['data']['z'])):
        if in_bin(data['data']['x'][i], x_bin) and in_bin(data['data']['y'][i], y_bin):
            data_bin['z'].append(data['data']['z'][i])
            data_bin['M'].append(data['data']['M'][i])
            data_bin['Mt'].append(data['data']['Mt'][i])

    data_bin['dz'] = [(data_bin['z'][i] + data_bin['z'][i+1])/2 for i in range(len(data_bin['z'])-1)]
    data_bin['dM'] = [(data_bin['Mt'][i+1] - data_bin['Mt'][i])/(data_bin['z'][i+1]-data_bin['z'][i]) for i in range(len(data_bin['z'])-1)]

    print('x: [%5.3f, %5.3f]\ty: [%5.3f, %5.3f]' % (x_bin[0], x_bin[1], y_bin[0], y_bin[1]))
    print('\tz: ' + str(data_bin['z']))
    print('\tM: ' + str(data_bin['M']))
    print('\tMt: ' + str(data_bin['Mt']))
    return data_bin


data = read_file("results.data")

n_x_bin = len(data['x_bin'])
n_y_bin = len(data['y_bin'])

fig = plt.figure()
gs = fig.add_gridspec(nrows=n_y_bin+1, ncols=n_x_bin+1, wspace=0, hspace=0,
                        width_ratios =[1 if in_bin(i, [1,n_x_bin  ]) else 0.5 for i in range(n_x_bin+1)],
                        height_ratios=[1 if in_bin(i, [0,n_y_bin-1]) else 0.5 for i in range(n_y_bin+1)])
axs = gs.subplots(sharex=False, sharey=False)

axs[n_y_bin][0].set_xlim(0, data['x_bin'][0][0])
axs[n_y_bin][0].set_ylim(0, data['y_bin'][0][0])
axs[n_y_bin][0].set_xticks([data['x_bin'][0][0]])
axs[n_y_bin][0].set_yticks([data['y_bin'][0][0]])

axs[0][int((n_x_bin+1)/2)].set_title("COMPASS pi+ multiplicities")
axs[int(n_y_bin/2)][0].set_ylabel("y")
axs[n_y_bin][int((n_x_bin+1)/2)].set_xlabel("x")

for i in range(n_x_bin):
    x_bin = data['x_bin'][i]
    ax = axs[n_y_bin][i+1]
    ax.set_xlim(x_bin[0], x_bin[1])
    ax.set_xticks([x_bin[1]])
    ax.set_yticks([])

for i in range(n_y_bin):
    y_bin = data['y_bin'][i]
    ax = axs[n_y_bin-i-1][0]
    ax.set_ylim(y_bin[0], y_bin[1])
    ax.set_yticks([y_bin[1]])
    ax.set_xticks([])


for i in range(n_y_bin):
    for j in reversed(range(n_x_bin)):
        x_bin = data['x_bin'][j]
        y_bin = data['y_bin'][i]

        ax = axs[n_y_bin-i-1][j+1]
        data_bin = select(data, x_bin, y_bin)


        ax.set_xlim(0.2, 0.8)
        #ax.set_ylim(0, 2)

        if (i == 0):
            ax.set_xticks([0.3, 0.5, 0.7])
            ax.set_xlabel("z")
        else:
            ax.set_xticks([0.5])

        if (j == 0):
            ax.set_yticks([0.5, 1.0, 1.5])
            #ax.set_ylabel("M")
        else:
            ax.set_yticks([])
        
#        if len(data_bin['z']) == 0:
#            ax.set_xticks([])
#            ax.set_yticks([])
#            continue

        ax.plot(data_bin['z'], data_bin['M'], '.')
        ax.plot(data_bin['z'], data_bin['Mt'], '-')
#        ax.plot(data_bin['dz'], data_bin['dM'], '-')


plt.show()
#plt.savefig("compass.png")
