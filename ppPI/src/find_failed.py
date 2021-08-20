
def replica_filename(exp, replica):
    replica_str = "%03d" % replica
    return "../theory/" + exp + "/results/1_3_3_1.0_1" + replica_str + ".res"


def search_file(exp, replica):
    with open(replica_filename(exp, replica)) as result_file:
        for line in result_file:
            #if line.split()[1] == "NaN":
            try:
                if float(line.split()[1]) < 0:
                    return True
            except ValueError:
                continue

    return False

def replace_file(exp, replica):
    data = []
    found = False
    with open(replica_filename(exp, replica), 'r') as result_file:
        for line in result_file:
            items = line.strip().split()
            if items[1] == "NaN":
                found = True
                data.append((float(items[0]), -1))
                print("%f : %f" % (float(items[0]), -1))
            else:
                data.append((float(items[0]), float(items[1])))

    if found:
        with open(replica_filename(exp, replica), 'w') as result_file:
            for line in data:
                result_file.write("%f\t%f\n" % (line[0], line[1]))

    return found




def main():
    ff_count = 201

    experiments = []
    with open("../data/successful_experiments.txt", 'r') as exp_file:
        for line in exp_file:
            experiments.append(line.strip())


    with open("failed_replicas.txt", 'w') as failed_file:
        for exp in experiments:
            pdf_replicas = [0 for i in range(ff_count)]
            with open("../theory/" + exp + "/pdf_replicas.txt", 'r') as replica_file:
                for line in replica_file:
                    items = [int(x) for x in line.strip().split()]
                    pdf_replicas[items[0]] = items[1]

            for i in range(201):
                #if search_file(exp, i):
                if replace_file(exp, i):
                    print("Failed: %s - %3d %3d" % (exp, i, pdf_replicas[i]))
                    failed_file.write("%s %d %d\n" % (exp, i, pdf_replicas[i]))


if __name__ == '__main__':
    main()
