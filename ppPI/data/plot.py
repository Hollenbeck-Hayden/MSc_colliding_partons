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
    'axes.titlesize': 'medium',
});

def savefig(path):
    plt.savefig(path + '.pgf', pad_inches=0)


def setfigsize(fig, width, height):
    fig.set_size_inches(w=width, h=height)
    fig.tight_layout()


def make_ratio_fig(rows, cols, h_ratios=[3,1]):
    fig = plt.figure()
    gs = fig.add_gridspec(nrows=rows, ncols=cols, height_ratios=h_ratios, hspace=0, wspace=0)
    axes = gs.subplots(sharex=True, sharey=False)
    return fig, axes

cutoff = 5

experiments = []
with open("all_experiments.txt", 'r') as infile:
    for line in infile:
        if line.startswith("#"): continue
        experiments.append(line.strip())


fig, axes = make_ratio_fig(1, 1, h_ratios=[1])

data = {}
for exp in experiments:
    data[exp] = ([], [])
    with open(exp + "/" + exp + "_unfiltered.txt", 'r') as infile:
        cme = 0
        i = 0
        for line in infile:
            if i < 5:
                if i == 2:
                    cme = float(line.strip())
                i += 1
            else:
                #if (float(line.split()[1]) < cutoff): continue
                data[exp][0].append(float(line.split()[1]))
                data[exp][1].append(cme)

    i = 0
    style = "."
    color = "red"
    exp_set = ""

    if exp.startswith("ALICE"):
        style = "."
        color = "red"
        exp_set = "ALICE"
    elif exp.startswith("CMS"):
        style = "x"
        color = "blue"
        exp_set = "CMS"
    elif exp.startswith("PHENIX"):
        style = "^"
        color = "green"
        exp_set = "PHENIX"
    elif exp.startswith("STAR"):
        style = "s"
        color = "orange"
        exp_set = "STAR"

    cut = [ [data[exp][0][i] for i in range(len(data[exp][0])) if data[exp][0][i] < cutoff],
            [data[exp][1][i] for i in range(len(data[exp][0])) if data[exp][0][i] < cutoff]]

    print("exp: %s, cut: %d" % (exp, len(cut[0])))


    uncut=[ [data[exp][0][i] for i in range(len(data[exp][0])) if data[exp][0][i] >= cutoff],
            [data[exp][1][i] for i in range(len(data[exp][0])) if data[exp][0][i] >= cutoff]]



    axes.scatter(  cut[0],   cut[1], color='grey', marker=style, label=exp_set)
    axes.scatter(uncut[0], uncut[1], color=color, marker=style, label=exp_set)

    #if exp.endswith("PIp") or exp.endswith("PIm"): i = 1
    #elif exp.endswith("PIm"): continue
    #elif exp.endswith("PI0"): i = 2
    #else: i = 3

    #axes[i].scatter(  cut[0],   cut[1], color='grey', marker=style, label=exp)
    #axes[i].scatter(uncut[0], uncut[1], color=color, marker=style, label=exp)


axes.set_xlabel(r"$p_T$ [GeV]")
axes.set_ylabel(r"$\sqrt{s}$ [GeV]")

setfigsize(fig, 5, 2.5)

#for ax in axes:
#    ax.set_xlabel("pT [GeV]")
#    ax.set_ylabel("CME [GeV]")
#    #ax.legend()

#axes[0].set_title("All")
#axes[1].set_title("P P --> PI+/- X")
#axes[2].set_title("P P --> PI0 X")
#axes[3].set_title("P P --> (PI+ + PI-) X")

savefig('kinematic_coverage')
#plt.show()
#plt.savefig('kinematic_coverage.png')
