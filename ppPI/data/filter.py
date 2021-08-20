# Collection of filters to transfor raw data into a common format

import pandas as pd

#mub_to_mb = 1e-3
#pb_to_mb = 1e-9

#mb_to_gev = 2.56819
#mub_to_gev = mb_to_gev * mub_to_mb
#pb_to_gev = mub_to_gev * pb_to_mb

cutoff = 3

successful = []

# Define data set filters

class FilterSet:
    def __init__(self, name):
        self.nameexp = name
        self.ndata = 0
        self.nsys = 0
        self.cme = 0
        self.ymin = 0
        self.ymax = 0
        self.sysnor = 0
        self.data = None

    def readcsv(self, filename, skip, cols, footer=0):
        infile = self.nameexp + '/' + filename
        self.data = pd.read_csv(
            infile,
            dtype={"user_ld": float},
            skiprows=skip,
            sep=',',
            header=None,
            usecols=cols,
            skipfooter=footer,
            na_values="-",
            engine='python').dropna()

    def write(self, stat_mult=False, sys_mult=False):
        print("experiment: " + self.nameexp)
        self.nsys = len(self.data.columns) + 1 - 3
        self.data[""] = self.data[self.data.columns[1]] * self.sysnor / 100.
        self.data.columns = ['pT', 'obs', 'stat'] + ['sys' + str(i+1) for i in range(self.nsys)]
        #self.data['obs'] *= conversion
        if (stat_mult):
            if '%' in str(self.data['stat'][0]):
                self.data['stat'] = self.data['stat'].str.rstrip('%').astype(float) / 100.
            self.data['stat'] = self.data['obs'] * self.data['stat']
        #else:
        #    self.data['stat'] *= conversion
        if (sys_mult):
            if '%' in str(self.data['sys1'][0]):
                self.data['sys1'] = self.data['sys1'].str.rstrip('%').astype(float) / 100.
            self.data['sys1'] = self.data['obs'] * self.data['sys1']
        #else:
        #    self.data['sys1'] *= conversion
        self.ndata = len(self.data['pT'])

        self.write_file("_unfiltered")
        self.data = self.data[self.data['pT'] >= cutoff].dropna().reset_index(drop=True)
        if self.data.empty:
            return False
        else:
            self.write_file("")
            return True

    def write_file(self, addon):
        with open(self.nameexp + '/' + self.nameexp + addon + '.txt', 'w') as f:
            print(self.nameexp, file=f)
            print(len(self.data), self.nsys, file=f)
            print(self.cme, file=f)
            print(self.ymin, self.ymax, file=f)
            pd.set_option("display.max_rows", None, "display.max_columns", None)
            print(self.data, file=f)

def absolute(y):
    return (-y, +y)


def filter_ALICE_2760a_PI0():
    f = FilterSet("ALICE_2760a_PI0")
    f.readcsv('HEPData-ins1296306-v1-Table_2.csv', 13, [0,3,4,6])
    f.cme = 2760.0
    #f.ymin, f.ymax = absolute(0.6)
    f.ymin, f.ymax = absolute(0)
    f.sysnor = 3.9    # ???
    if f.write():
        successful.append(f.nameexp)

def filter_ALICE_2760b_PI0():
    f = FilterSet("ALICE_2760b_PI0")
    f.readcsv('HEPData-ins1512110-v1-Table_3.csv', 13, [0,3,4,6])
    f.cme = 2760.0
    f.ymin, f.ymax = absolute(0.6)
    f.sysnor = 5.7
    if f.write():
        successful.append(f.nameexp)

def filter_ALICE_2760_PIsum():
    f = FilterSet("ALICE_2760_PIsum")
    f.readcsv('HEPData-ins1276299-v1-Table_1.csv', 152, [0,3,4,6])
    f.cme = 2760.0
    f.ymin, f.ymax = absolute(0.8)
    f.sysnor = 7.6
    if f.write():
        successful.append(f.nameexp)

def filter_ALICE_5020_PIsum():
    f = FilterSet("ALICE_5020_PIsum")
    f.readcsv('HEPData-ins1759506-v1-Table_2.csv', 12, [0,3,4,6])
    f.cme = 5020.0
    f.ymin, f.ymax = absolute(0.5)
    f.sysnor = 6.4
    if f.write():
        successful.append(f.nameexp)

def filter_ALICE_7000_PI0():
    f = FilterSet("ALICE_7000_PI0")
    f.readcsv('HEPData-ins1116147-v1-Table_1.csv', 13, [0,3,4,6])
    f.cme = 7000.0
    f.ymin = 0
    f.ymax = 0
    f.sysnor = 6.25
    #f.write(conversion=mub_to_mb)
    if f.write():
        successful.append(f.nameexp)

def filter_ALICE_7000_PIsum():
    f = FilterSet("ALICE_7000_PIsum")
    f.readcsv('ALICE_7000_PI.csv', 13, [0,3,4,6], footer=183-67)
    f.cme = 7000.0
    f.ymin, f.ymax = absolute(0.5)
    f.sysnor = 5.5
    if f.write():
        successful.append(f.nameexp)

def filter_ALICE_8000_PI0():
    f = FilterSet("ALICE_8000_PI0")
    f.readcsv('HEPData-ins1620477-v2-Table_1.csv', 13, [0,3,4,6])
    f.cme = 8000.0
    f.ymin = 0
    f.ymax = 0      # 0.8?
    f.sysnor = 2.6
    #f.write(conversion=pb_to_mb)
    if f.write():
        successful.append(f.nameexp)

def filter_ALICE_900_PIp():
    f = FilterSet("ALICE_900_PIp")
    f.readcsv('HEPData-ins885104-v1-Table_1.csv', 12, [0,3,4,6], footer=83-46)
    f.cme = 900.0
    f.ymin, f.ymax = absolute(0.5)
    f.sysnor = 3.6
    if f.write():
        successful.append(f.nameexp)

def filter_ALICE_900_PIm():
    f = FilterSet("ALICE_900_PIm")
    f.readcsv('HEPData-ins885104-v1-Table_1.csv', 49, [0,3,4,6])
    f.cme = 900.0
    f.ymin, f.ymax = absolute(0.5)
    f.sysnor = 3.6
    if f.write():
        successful.append(f.nameexp)

def filter_CMS_900_PIp():
    f = FilterSet("CMS_900_PIp")
    f.readcsv('HEPData-ins1123117-v1-Table_1.csv', 13, [0,1,2,4], footer=120-46)
    f.cme = 900.0
    f.ymin, f.ymax = absolute(1)
    f.sysnor = 3.0
    if f.write():
        successful.append(f.nameexp)

def filter_CMS_900_PIm():
    f = FilterSet("CMS_900_PIm")
    f.readcsv('HEPData-ins1123117-v1-Table_2.csv', 13, [0,1,2,4], footer=120-46)
    f.cme = 900.0
    f.ymin, f.ymax = absolute(1)
    f.sysnor = 3.0
    if f.write():
        successful.append(f.nameexp)

def filter_CMS_2760_PIp():
    f = FilterSet("CMS_2760_PIp")
    f.readcsv('HEPData-ins1123117-v1-Table_3.csv', 13, [0,1,2,4], footer=120-46)
    f.cme = 2760.0
    f.ymin, f.ymax = absolute(1)
    f.sysnor = 3.0
    if f.write():
        successful.append(f.nameexp)

def filter_CMS_2760_PIm():
    f = FilterSet("CMS_2760_PIm")
    f.readcsv('HEPData-ins1123117-v1-Table_4.csv', 13, [0,1,2,4], footer=120-46)
    f.cme = 2760.0
    f.ymin, f.ymax = absolute(1)
    f.sysnor = 3.0
    if f.write():
        successful.append(f.nameexp)

def filter_CMS_7000_PIp():
    f = FilterSet("CMS_7000_PIp")
    f.readcsv('HEPData-ins1123117-v1-Table_5.csv', 13, [0,1,2,4], footer=120-46)
    f.cme = 7000.0
    f.ymin, f.ymax = absolute(1)
    f.sysnor = 3.0
    if f.write():
        successful.append(f.nameexp)

def filter_CMS_7000_PIm():
    f = FilterSet("CMS_7000_PIm")
    f.readcsv('HEPData-ins1123117-v1-Table_6.csv', 13, [0,1,2,4], footer=120-46)
    f.cme = 7000.0
    f.ymin, f.ymax = absolute(1)
    f.sysnor = 3.0
    if f.write():
        successful.append(f.nameexp)

def filter_PHENIX_200a_PIp():
    f = FilterSet("PHENIX_200a_PIp")
    f.readcsv('HEPData-ins886590-v1-Table_1.csv', 13, [0,3,4,6], footer=72-41)
    f.cme = 200
    f.ymin, f.ymax = absolute(0.6)
    f.sysnor = 9.2
    if f.write():
        successful.append(f.nameexp)

def filter_PHENIX_200a_PIm():
    f = FilterSet("PHENIX_200a_PIm")
    f.readcsv('HEPData-ins886590-v1-Table_1.csv', 44, [0,3,4,6])
    f.cme = 200
    f.ymin, f.ymax = absolute(0.6)
    f.sysnor = 9.2
    if f.write():
        successful.append(f.nameexp)

def filter_PHENIX_200b_PIp():
    f = FilterSet("PHENIX_200b_PIp")
    f.readcsv('HEPData-ins1315330-v1-Table_1.csv', 13, [0,3,4,6], footer=31-20)
    f.cme = 200
    f.ymin, f.ymax = absolute(0.35)
    f.sysnor = 9.6
    if f.write():
        successful.append(f.nameexp)

def filter_PHENIX_200b_PIm():
    f = FilterSet("PHENIX_200b_PIm")
    f.readcsv('HEPData-ins1315330-v1-Table_1.csv', 24, [0,3,4,6])
    f.cme = 200
    f.ymin, f.ymax = absolute(0.35)
    f.sysnor = 9.6
    if f.write():
        successful.append(f.nameexp)

def filter_PHENIX_200_PI0():
    f = FilterSet('PHENIX_200_PI0')
    f.readcsv('HEPData-ins617784-v1-Table_1.csv', 13, [0,3,4,6])
    f.cme = 200.0
    f.ymin, f.ymax = absolute(0.35)
    f.sysnor = 9.6  #%
    #f.write(conversion=mb_to_gev, stat_mult=True, sys_mult=True)
    if f.write(stat_mult=True, sys_mult=True):
        successful.append(f.nameexp)

def filter_PHENIX_62_PIp():
    f = FilterSet('PHENIX_62_PIp')
    f.readcsv('HEPData-ins886590-v1-Table_5.csv', 13, [0,3,4,6], footer=70-40)
    f.cme = 62
    f.ymin, f.ymax = absolute(0.6)
    f.sysnor = 11
    if f.write():
        successful.append(f.nameexp)

def filter_PHENIX_62_PIm():
    f = FilterSet('PHENIX_62_PIm')
    f.readcsv('HEPData-ins886590-v1-Table_5.csv', 43, [0,3,4,6])
    f.cme = 62
    f.ymin, f.ymax = absolute(0.6)
    f.sysnor = 11
    if f.write():
        successful.append(f.nameexp)

def filter_STAR_200a_PI0():
    f = FilterSet('STAR_200a_PI0')
    f.readcsv('HEPData-ins836952-v1-Figure_1a.csv', 9, [0,3,4,6])
    f.cme = 200.0
    f.ymin = 0.
    f.ymax = 1.0
    f.sysnor = 11.7
    if f.write():
        successful.append(f.nameexp)

def filter_STAR_200b_PI0():
    f = FilterSet('STAR_200b_PI0')
    f.readcsv('HEPData-ins840766-v1-Figure_22a.csv', 9, [0,3,4,6], footer=39-23)
    f.cme = 200.0
    f.ymin = 0.8
    f.ymax = 2.0
    f.sysnor = 7.7
    if f.write():
        successful.append(f.nameexp)


#def filter_PHENIX_200b_PI0():
#    f = FilterSet('PHENIX_200b_PI0')
#    f.readcsv('HEPData-ins617784-v1-Table_1.csv', 13, [0,3,4,6])
#    f.cme = 200.0
#    f.ymin, f.ymax = absolute(0.35)
#    f.sysnor = 9.6
#    if f.write(stat_mult=True, sys_mult=True):
#        successful.append(f.nameexp)

def filter_PHENIX_510_PI0():
    f = FilterSet('PHENIX_510_PI0')
    f.readcsv('data.csv', 1, [0,1,2,3])
    f.cme = 510.0
    f.ymin, f.ymax = absolute(0.35)
    f.sysnor = 10
    if f.write():
        successful.append(f.nameexp)


filter_ALICE_2760a_PI0()
filter_ALICE_2760b_PI0()
filter_ALICE_2760_PIsum()
filter_ALICE_5020_PIsum()
filter_ALICE_7000_PI0()
filter_ALICE_7000_PIsum()
filter_ALICE_8000_PI0()
filter_ALICE_900_PIp()
filter_ALICE_900_PIm()

filter_CMS_2760_PIp()
filter_CMS_2760_PIm()
filter_CMS_7000_PIp()
filter_CMS_7000_PIm()
filter_CMS_900_PIp()
filter_CMS_900_PIm()

filter_PHENIX_200a_PIp()
filter_PHENIX_200a_PIm()
filter_PHENIX_200b_PIp()
filter_PHENIX_200b_PIm()
filter_PHENIX_200_PI0()
filter_PHENIX_62_PIp()
filter_PHENIX_62_PIm()

filter_STAR_200a_PI0()
filter_STAR_200b_PI0()

filter_PHENIX_510_PI0()

with open('successful_experiments.txt', 'w') as outfile:
    for exp in successful:
        outfile.write(exp + "\n")
