"""
@author: Seth

for plotting .csv convergence files
"""

import matplotlib.pyplot as plt
import csv
from math import sqrt

def get_csv_data(filename, col):
    data = []
    times = []
    file = open(filename, 'rU')
    reader = csv.reader(file)
    get_data = False
    for row in reader:
        if get_data:
            times.append(float(row[0]))
            data.append(float(row[col]))
        elif len(row) > 0:
            if row[0] == 'time':
                get_data = True
    file.close()
    return [times, data]

def get_csv_data_skip(filename, col, skip):
    data = []
    times = []
    file = open(filename, 'rU')
    reader = csv.reader(file)
    get_data = 0
    for row in reader:
        get_data += 1
        if get_data > skip:
            times.append(float(row[0]))
            data.append(float(row[col]))
    file.close()
    return [times, data]

def get_csv_data_c1c2skip(filename, c1, c2, skip):
    data = []
    times = []
    file = open(filename, 'rU')
    reader = csv.reader(file)
    get_data = 0
    for row in reader:
        get_data += 1
        if get_data > skip:
            times.append(float(row[c1]))
            data.append(float(row[c2]))
    file.close()
    return [times, data]

def plot_Qdata(data, names):
    for i in range(len(data)):
        plt.plot(data[i][0], data[i][1], label=names[i])
    plt.xlabel('time')
    plt.ylabel('Q(t)')
    plt.legend(loc='best')
    return

def plot_Mdata(data, names):
    for i in range(len(data)):
        plt.plot(data[i][0], data[i][1], label=names[i])
    plt.xlabel('time')
    plt.ylabel('M(t, R)')
    plt.legend(loc='best')
    return

def plot_Mnames(names, c1, c2, skip):
    data = []
    for name in names:
        data = [get_csv_data_c1c2skip(name, c1, c2, skip)]
        plt.plot(data[0], data[1], label=name)
    plt.xlabel('time')
    plt.ylabel('M(t, R)')
    plt.legend(loc='best')
    return


""" input: full filename (with extension) as string """
def get_i1data(filename, i1=1):
    file = open(filename, 'rU')
    reader = csv.reader(file)
    data = [row[i1] for row in reader]
    file.close()
    return data

def get_i1i2data(filename, i1=1, i2=3):
    file = open(filename, 'rU')
    reader = csv.reader(file)
    data = [[row[i1], row[i2]] for row in reader]
    file.close()
    return data

""" ([[str,str]]) data = [[m(t), m_{r>r0}(t)] for all t] """
def get_mi(data, i0=1):
    return float(data[i0][0]) - float(data[len(data)/2][1])

def get_mf(data):
    return float(data[len(data)-1][0]) - float(data[len(data)-1][1])

"""d0 = val(4h), d1 = val(2h), d2 = val(h)"""
def get_err(d0, d1, d2):
    e01 = (d0 - d1) / 12.0
    e12 = (d1 - d2) / 3.0
    return max(abs(e01), abs(e12))

def make_mfname(fname, resn):
    return "mass-" + str(resn) + "-" + fname + ".csv"

def make_ricci_fname(fname, resn):
    return "maxRicci-" + str(resn) + "-" + fname + ".csv"

""" give fname as just the outfile variable from the simulation
    and for fout give full file name or '' to prevent writing 
    returns [mass(t(i0)), mass_err(t(i0)), mass(tf), mass_err(tf)]"""
def mass_err(fname, resn0=8, fout='', i0=1):
    filename0 = make_mfname(fname, resn0)
    m0 = get_i1data(filename0, 1)
    
    filename1 = make_mfname(fname, 2*resn0)
    m1 = get_i1i2data(filename1, 1, 3)
        
    filename2 = make_mfname(fname, 4*resn0)
    m2 = get_i1data(filename2, 1)

    mi = get_mi(m1, i0)
    ei = get_err(float(m0[i0]), float(m1[i0][0]), float(m2[i0]))
    mf = get_mf(m1)
    k = len(m1) - 1
    ef = get_err(float(m0[k]), float(m1[k][0]), float(m2[k]))
    mlist = [mi, ei, mf, ef]

    if fout:
        times = get_i1data(filename0, 0)
        fwr = open(fout, 'w')
        writer = csv.writer(fwr)
        writer.writerow(["time","(4h-2h)/12","(2h-h)/3","max(|err|)"])
        for i in range(i0, len(m2)):
            wrow = [float(times[i])]
            wrow.append((float(m0[i]) - float(m1[i][0])) / 12.0)
            wrow.append((float(m1[i][0]) - float(m2[i])) / 3.0)
            wrow.append(max(abs(wrow[1]), abs(wrow[2])))
            writer.writerow(wrow)
        fwr.close()
    
    return mlist

def absorbed(mlist):
    mi = mlist[0]
    ei = mlist[1]
    mf = mlist[2]
    ef = mlist[3]
    absorb = 1 - (mf / mi)
    lob = 1 - ((mf + ef) / (mi - ei))
    elo = absorb - lob
    upb = 1 - ((mf - ef) / (mi + ei))
    eup = upb - absorb
    return [absorb, elo, eup]
    
def merrplot(fnames, dsqs, resn0=8, fout='', i0=1):
    x = [sqrt(float(dsq)) for dsq in dsqs]
    y = []
    errlo = []
    errup = []
    for fname in fnames:
        data = absorbed(mass_err(fname, resn0, fout, i0))
        y.append(100*data[0])
        errlo.append(400*data[1])
        errup.append(400*data[2])
    plt.errorbar(x, y, yerr=[errlo, errup], lw=1, elinewidth=2, capsize=2)
    plt.xlabel('initial pulse width / black hole mass')
    plt.ylabel('mass absorption percent')
    plt.title('scalar field incident on static black hole:\n' \
              'mass absorption vs. pulse width')
    return "enter plt.show() to view plot"

def ricci_max(fname, resn):
    filename = make_ricci_fname(fname, resn)
    data = get_csv_data_skip(filename, 2, 1)
    return max(data[1])

def get_ricci_max_data(fnames, amps, resn=16):
    max_data = []
    amp_data = []
    for k in range(len(fnames)):
        max_data.append(ricci_max(fnames[k],resn))
        amp_data.append(amps[k])
    return [amp_data, max_data]

def ricci_max_plot(fnames, amps, resn=8, dsq=25):
    plot_data = get_ricci_max_data(fnames, amps, resn)
    plt.plot(plot_data[0], plot_data[1])
    plt.xlabel('Initial Amplitude')
    plt.ylabel('Ricci Scalar Maximum Value')
    titletext = 'max_Ricci vs. ic_Amplitude for ic_Dsq = ' + str(dsq) + ', resn = ' + str(resn)
    plt.title(titletext)
    return

use_names = ["_0100", "_0102", "_0104", "_0106", "_0108", "_0110", "_0111", "_0112", "_0113", "_0114", "_0115", "_0116", "_0117", "_0118", "_0119", "_01200", "_01202", "_01204", "_01206", "_01208", "_01210", "_01212", "_01214", "_01216", "_01218", "_01220", "_01222", "_01224", "_01226", "_01228"]

use_amps = [0.0100, 0.0102, 0.0104, 0.0106, 0.0108, 0.0110, 0.0111, 0.0112, 0.0113, 0.0114, 0.0115, 0.0116, 0.0117, 0.0118, 0.0119, 0.01200, 0.01202, 0.01204, 0.01206, 0.01208, 0.01210, 0.01212, 0.01214, 0.01216, 0.01218, 0.01220, 0.01222, 0.01224, 0.01226, 0.01228]

run_default = raw_input("run default plot? (y/n)\n")
if run_default == "y":
    ricci_max_plot(use_names, use_amps)


