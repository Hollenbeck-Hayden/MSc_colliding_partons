
def normalization(exp):
    pT = []
    obs = []

    with open(exp + "/results/1_3_3_1.0_1000.res", 'r') as theory_file:
        for line in theory_file:
            item = [float(x) for x in line.strip().split()]
            pT.append(item[0])
            obs.append(item[1])

    intervals = [pT[i] - pT[i-1] for i in range(1, len(pT))]
    midpoints = [(obs[i] + obs[i-1])/2.0 for i in range(1, len(obs))]

    total = 0
    for i in range(len(intervals)):
        total += intervals[i] * midpoints[i]
    return total




mb_to_mub = 1e3
mb_to_pb = 1e9

#gev_to_mb = 0.3894
#gev_to_mub = gev_to_mb * mb_to_mub
#gev_to_pb  = gev_to_mb * mb_to_pb

experiments = []
with open("../data/successful_experiments.txt", 'r') as exp_file:
    for line in exp_file:
        if len(line) == 0:
            continue
        if line.startswith("#"):
            continue
        experiments.append(line.strip())

for exp in experiments:
    observable = 0
    with open(exp + "/" + exp + ".txt", 'r') as config_file:
        i = 0
        for line in config_file:
            if i == 6:
                observable = int(line.strip().split()[0])
            i += 1

    conversion = 1
    units = "mb"
    if exp == "ALICE_7000_PI0":
        conversion = mb_to_mub
        units = r"\mu b"
        conversion /= 2
    elif exp == "ALICE_8000_PI0" or exp == "ALICE_2760a_PI0" or exp == "ALICE_2760b_PI0":
        conversion = mb_to_pb
        units = r"pb"
        conversion /= 2

    with open(exp + "/post_process.txt", 'w') as post_file:
        post_file.write("obs: ")
        if observable == 1:
            post_file.write(r"\frac{d^2 \sigma}{dy \, dp_T^2} \quad [" + units + r"/GeV^2]" + "\n")
        if observable == 2:
            post_file.write(r"E\frac{d^3 \sigma}{dp^3} \quad [" + units + r"/GeV^2]" + "\n")
        if observable == 3:
            post_file.write(r"\frac{d^2 \sigma}{dy \, dp_T} \quad [" + units + r"/Gev^2]" + "\n")

        post_file.write("conv: %f\n" % conversion)

        items = exp.split('_')
        post_file.write('exp: %s\n' % items[0])
        post_file.write('name: ')
        if items[2] == 'PI0':
            post_file.write(r'$\pi^0$')
        elif items[2] == 'PIsum':
            post_file.write(r'$\pi^+ + \pi^-$')
        if items[1][-1] == 'a' or items[1][-1] == 'b':
            items[1] = items[1][:-1]
        post_file.write(' %s GeV\n' % items[1])
