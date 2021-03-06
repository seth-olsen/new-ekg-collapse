"""
@author: Seth

for plotting .csv convergence files
"""

import matplotlib.pyplot as plt
import csv
from math import sqrt
from math import log
import numpy as np
from scipy.optimize import curve_fit

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

# **************************************************************************
# **************************************************************************


# **************************************************************************
# **************************************************************************
# ************************* SCALAR COLLLAPSE *******************************
# **************************************************************************
# **************************************************************************


# **************************************************************************
# **************************************************************************

def make_ricci_fname(fname, resn):
    return "maxRicci-" + str(resn) + "-" + fname + ".csv"

def ricci_max(fname, resn):
    filename = make_ricci_fname(fname, resn)
    data = get_csv_data_skip(filename, 2, 1)
    return max(data[1])

def get_ricci_max_data(fnames, amps, resn=8):
    max_data = []
    amp_data = []
    for k in range(len(fnames)):
        max_data.append(ricci_max(fnames[k],resn))
        amp_data.append(amps[k])
    return [amp_data, max_data]

def ricci_max_plot(fnames, amps, resn=8, dsq=25):
    plot_data = np.array(get_ricci_max_data(fnames, amps, resn))
    plt.plot(plot_data[0], plot_data[1], 'g^')
    plt.xlabel('Initial Amplitude')
    plt.ylabel('Ricci Scalar Maximum Value')
    titletext = 'max_Ricci vs. ic_Amplitude for ic_Dsq = ' + str(dsq) + ', resn = ' + str(resn)
    plt.title(titletext)
    return plot_data

def ricci_max_crit_plot(fnames, amps, pcrit, resn=8, dsq=25, logplot=False, marker='bo'):
    plot_data = np.array(get_ricci_max_data(fnames, amps, resn))
    plot_data[0] = pcrit - plot_data[0]
    if (logplot):
        plot_data = np.log(plot_data)
    plt.plot(plot_data[0], plot_data[1], marker)
    if (logplot):
        plt.xlabel('log(Critical Amplitude - Initial Amplitude)')
        plt.ylabel('log(Ricci Scalar Maximum)')
    else:
        plt.xlabel('Critical Amplitude - Initial Amplitude')
        plt.ylabel('Ricci Scalar Maximum')
    titletext = 'max_Ricci vs. ic_Amplitude for ic_Dsq = ' + str(dsq) + ', resn = ' + str(resn)
    plt.title(titletext)
    return plot_data

def get_ricci_max_crit_data(fnames, amps, pcrit, emax=[0,0,0], resn=8, nparr=False):
    max_lob = []
    max_upb = []
    amp_lob = []
    amp_upb = []
    max_data = []
    amp_data = []
    for k in range(len(fnames)):
        amp = amps[k]
        # ***** ADJUST THESE*********************************************
        amp_data.append(pcrit - amp)
        amp_lob.append(0.012285 - amp)
        amp_upb.append(0.012291 - amp)
        max_ric = ricci_max(fnames[k],resn)
        max_data.append(max_ric)
        if (amp < 0.0122):
            max_lob.append(max_ric - emax[0])
            max_upb.append(max_ric + emax[0])
        elif (amp < 0.01227):
            max_lob.append(max_ric - emax[1])
            max_upb.append(max_ric + emax[1])
        else:
            max_lob.append(max_ric - emax[2])
            max_upb.append(max_ric + emax[2])
    if (nparr):
        return np.array([amp_data, amp_lob, amp_upb, max_data, max_lob, max_upb])
    else:
        return [amp_data, amp_lob, amp_upb, max_data, max_lob, max_upb]

def get_ricci_max_log_crit_data(fnames, amps, pcrit, emax=[0,0,0], resn=8):
    return np.log(get_ricci_max_crit_data(fnames,amps,pcrit,emax,resn,True))

def ricci_max_errplot(logplot, fnames, amps, emax=[0,0,0], pcrit=0.0122885, resn=8, dsq=25, marker='g^'):
    plot_data = None
    if (logplot):
        plot_data = get_ricci_max_log_crit_data(fnames, amps, emax, pcrit, resn)
    else:
        plot_data = get_ricci_max_crit_data(fnames, amps, emax, pcrit, resn)
    x = plot_data[0]
    y = plot_data[3]
    xerrlo = []
    xerrup = []
    yerrlo = []
    yerrup = []
    for k in range(len(x)):
        xerrlo.append(x[k] - plot_data[1][k])
        xerrup.append(plot_data[2][k] - x[k])
        yerrlo.append(y[k] - plot_data[4][k])
        yerrup.append(plot_data[5][k] - y[k])
    plt.errorbar(x, y, marker, xerr=[xerrlo, xerrup], yerr=[yerrlo, yerrup], lw=1, elinewidth=2, capsize=2)
    if (logplot):
        plt.xlabel('Log(Critical Amplitude - Initial Amplitude)')
        plt.ylabel('Log(Ricci Scalar Maximum Value)')
    else:
        plt.xlabel('Critical Amplitude - Initial Amplitude')
        plt.ylabel('Ricci Scalar Maximum Value')
    titletext = 'max_Ricci vs. ic_Amplitude for ic_Dsq = ' + str(dsq) + ', resn = ' + str(resn)
    plt.title(titletext)
    return

def pwr_fit_fn(xdata, gam, c):
    return c * np.power(xdata, -2*gam)

def log_fit_fn(xdata, gam, c):
    return -2*gam*xdata + c

def pcrit_fit_fn(xdata, pcrit, gam, c):
    return c * np.power(pcrit - xdata, -2*gam)

def pcrit_log_fit_fn(xdata, pcrit, gam, c):
    return c - 2*gam*np.log(pcrit - xdata)

def fit_ricci_data(fnames, amps, pcrit=0.012288, emax=[0,0,0], resn=8, pwr_fit=False, pcrit_fit=False, guess=None):
    popt = None
    pcov = None
    if (pcrit_fit):
        if (pwr_fit):
            plot_data = get_ricci_max_data(fnames, amps, resn)
            xdata = np.array(plot_data[0])
            ydata = np.array(plot_data[1])
            popt, pcov = curve_fit(pcrit_fit_fn, xdata, ydata, p0=guess)
        else:
            plot_data = get_ricci_max_data(fnames, amps, resn)
            xdata = np.array(plot_data[0])
            ydata = np.log(np.array(plot_data[1]))
            popt, pcov = curve_fit(pcrit_log_fit_fn, xdata, ydata, p0=guess)
    else:
        if (pwr_fit):
            plot_data = get_ricci_max_crit_data(fnames, amps, pcrit, emax, resn)
            xdata = np.array(plot_data[0])
            ydata = np.array(plot_data[3])
            popt, pcov = curve_fit(pwr_fit_fn, xdata, ydata, p0=guess)
        else:
            plot_data = get_ricci_max_log_crit_data(fnames, amps, pcrit, emax, resn)
            xdata = np.array(plot_data[0])
            ydata = np.array(plot_data[3])
            popt, pcov = curve_fit(log_fit_fn, xdata, ydata, p0=guess)
    return [popt, pcov]


def make_amp_file_name(digits):
    dstr = str(digits)
    namestr = "_01" + dstr
    numstr = "0.01" + dstr
    return [namestr, float(numstr)]

def append_amp_file_names(namelist, amplist, dstart, dstop, dstep=1):
    for d in range(dstart, dstop, dstep):
        entry = make_amp_file_name(d)
        namelist.append(entry[0])
        amplist.append(entry[1])
    return 0

def generate_amp_file_names(ranges=[[200,225,2], [225,227,1], [2255,2270,10], \
                                    [2270,2289,2], [2275,2289,2], [22775,22890,10]]):
    namelist = []
    amplist = []
    for r in ranges:
        for d in range(r[0], r[1], r[2]):
            entry = make_amp_file_name(d)
            namelist.append(entry[0])
            amplist.append(entry[1])
    return [namelist, amplist]


all_names = ["_0100", "_0102", "_0104", "_0106", "_0108", "_0110", \
             "_0111", "_0112", "_0113", "_0114", "_0115", \
             "_0116", "_0117", "_0118", "_0119","_01200", \
             "_01202", "_01204", "_01206", "_01208", "_01210", \
             "_01212", "_01214", "_01216", "_01218", "_01220", \
             "_01222", "_01224", "_01225", "_012255", "_01226", \
             "_012265", "_012268", "_012270", "_012272", "_012274", "_012275", \
             "_012276", "_012277", "_0122775", "_012278", "_0122785", "_012279", \
             "_0122795", "_012280", "_0122805", "_012281", "_0122815", "_012282", \
             "_0122825", "_012283", "_0122835", "_012284", "_0122845", "_012285", \
             "_0122855", "_012286", "_0122865", "_012287", "_0122875", "_012288"]

all_amps = [0.0100, 0.0102, 0.0104, 0.0106, 0.0108, 0.0110, \
            0.0111, 0.0112, 0.0113, 0.0114, 0.0115, \
            0.0116, 0.0117, 0.0118, 0.0119, 0.01200, \
            0.01202, 0.01204, 0.01206, 0.01208, 0.01210, \
            0.01212, 0.01214, 0.01216, 0.01218, 0.01220, \
            0.01222, 0.01224, 0.01225, 0.012255, 0.01226, \
            0.012265, 0.012268, 0.012270, 0.012272, 0.012274, 0.012275, \
            0.012276, 0.012277, 0.0122775, 0.012278, 0.0122785, \
            0.012279, 0.0122795, 0.012280, 0.0122805, \
            0.012281, 0.0122815, 0.012282, 0.0122825, 0.012283, 0.0122835, \
            0.012284, 0.0122845, 0.012285, 0.0122855, \
            0.012286, 0.0122865, 0.012287, 0.0122875, 0.012288]

nms = ["_01224", "_01225", "_012255", "_01226", \
       "_012265", "_012268", "_012270", "_012272", "_012274", "_012275", \
       "_012276", "_012277", "_0122775", "_012278", "_0122785", "_012279", \
       "_0122795", "_012280", "_0122805", "_012281", "_0122815", "_012282", \
       "_0122825", "_012283", "_0122835", "_012284", "_0122845", "_012285", \
       "_0122855", "_012286"]

ams = [0.01224, 0.01225, 0.012255, 0.01226, \
       0.012265, 0.012268, 0.012270, 0.012272, 0.012274, 0.012275, \
       0.012276, 0.012277, 0.0122775, 0.012278, 0.0122785, \
       0.012279, 0.0122795, 0.012280, 0.0122805, \
       0.012281, 0.0122815, 0.012282, 0.0122825, 0.012283, 0.0122835, \
       0.012284, 0.0122845, 0.012285, 0.0122855, 0.012286]

sm_nms = ["_012277", "_0122775", "_012278", "_0122785", \
          "_012279", "_0122795", "_012280", "_0122805", "_012281", "_0122815", \
          "_012282", "_0122825", "_012283"]

sm_ams = [0.012277, 0.0122775, 0.012278, 0.0122785, \
          0.012279, 0.0122795, 0.012280, 0.0122805, 0.012281, 0.0122815, \
          0.012282, 0.0122825, 0.012283]

pcrit_arr = 0.0000001 * np.array(range(122835, 122935, 5))
iend = len(nms)
if (iend != len(ams)):
    print("ERROR: names and amps differ")

gamma_results = []
best_results = []

for i in range(18):
    print("\n" + str(i) + "\n") 
    for j in range(10):
        print("******** off front = " + str(i) + "********")
        print("******** off back = " + str(j) + "********\n")
        for pcr in pcrit_arr:
            if (pcr > max(ams[i:(iend - j)])):
                params = fit_ricci_data(nms[i:(iend - j)], ams[i:(iend - j)], \
                                        pcr, emax=[0,0,0], resn=8, pwr_fit=False, \
                                        pcrit_fit=False, guess=[0.37, -3])
                gam = params[0][0]
                if ((gam < 0.375) and (gam > 0.365)):
                    print("gamma_fit = " + str(gam))
                    print("FOR pcrit = " + str(pcr) + "\n")
                    gamma_results.append([gam, pcr, i, j, params])
                    if (j < 5):
                        best_results.append([gam, pcr, i, j, params])


#for i in range(8):
#    print("\n" + str(i) + "\n") 
#    for j in range(4):
#        print("******** off front = " + str(i) + "********")
#        print("******** off back = " + str(j) + "********\n")
#        for pcr in pcrit_arr:
#            params = fit_ricci_data(sm_nms[i:(iend - j)], sm_ams[i:(iend - j)], \
#                                    pcr, emax=[0,0,0], resn=8, pwr_fit=False, \
#                                    pcrit_fit=False, guess=[0.37, -3])
#            gam = params[0][0]
#            gamma_results.append([gam, pcr, i, j])
#            if ((gam < 0.42) and (gam > 0.32)):
#                print("gamma_fit = " + str(gam))
#                print("FOR pcrit = " + str(pcr) + "\n")

#for i in range(10):
    #print(i)
    #all_names.pop(len(all_names)-1)
    #all_amps.pop(len(all_amps)-1)
    #print("fit pcrit:")
    #print(fit_ricci_data(all_names,all_amps,pcrit=0.012288, emax=[0,0,0], resn=8, pwr_fit=False, pcrit_fit=True, guess=[0.012288, 0.37, 1]))
    #print("given pcrit:")
#print(fit_ricci_data(all_names, all_amps, pcrit=0.012288, emax=[0,0,0], resn=8, pwr_fit=False, pcrit_fit=False, guess=[0.37, 1]))






