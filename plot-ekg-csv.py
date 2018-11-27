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

def get_ricci_max_crit_data(fnames, amps, pcrit, emax=[0.02,0.8,3.1], resn=8, nparr=False):
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
        amp_lob.append(pcrit - 0.0000005 - amp)
        amp_upb.append(pcrit + 0.0000005 - amp)
        max_ric = ricci_max(fnames[k],resn)
        max_data.append(max_ric)
        if (amp < 0.0122):
            max_lob.append(max_ric - emax[0])
            max_upb.append(max_ric + emax[0])
        elif (amp < 0.01225):
            max_lob.append(max_ric - emax[1])
            max_upb.append(max_ric + emax[1])
        else:
            max_lob.append(max_ric - emax[2])
            max_upb.append(max_ric + emax[2])
    if (nparr):
        return np.array([amp_data, amp_lob, amp_upb, max_data, max_lob, max_upb])
    else:
        return [amp_data, amp_lob, amp_upb, max_data, max_lob, max_upb]

def get_ricci_max_log_crit_data(fnames, amps, pcrit, emax=[0.02,0.8,3.1], resn=8):
    return np.log(get_ricci_max_crit_data(fnames,amps,pcrit,emax,resn,True))

def fit_label(params, pcrit):
    s = "$log(R_{max}) = C - 2 \gamma \,log(A_{crit} - A)$"
    s += "\n  $A_{crit} = " + str(pcrit)[:9] + " \pm 0.000005$"
    gam = params[0][0]
    c = params[0][1]
    err_gam = sqrt(params[1][0][0])
    err_c = sqrt(params[1][1][1])
    s += "\n  $\gamma = " + str(gam)[:8] + " \pm " + str(2*err_gam)[:8] + "$"
    s += "\n  $C = " + str(c)[:8] + " \pm " + str(2*err_c)[:8] + "$"
    return s

def ricci_max_errplot(fnames, amps, pcrit, emax=[0.02,0.8,3.1], resn=8, dsq=25, dmarker='bo', logplot=True):
    plot_data = None
    fit_p = None
    if (logplot):
        plot_data = get_ricci_max_log_crit_data(fnames, amps, pcrit, emax, resn)
        xdata = np.array(plot_data[0])
        ydata = np.array(plot_data[3])
        popt, pcov = curve_fit(log_fit_fn, xdata, ydata, p0=[0.37, -3])
        fit_p = [popt, pcov]
    else:
        plot_data = get_ricci_max_crit_data(fnames, amps, pcrit, emax, resn)
    x = np.array(plot_data[0])
    y = np.array(plot_data[3])
    xerrlo = []
    xerrup = []
    yerrlo = []
    yerrup = []
    for k in range(len(x)):
        xerrlo.append(x[k] - plot_data[1][k])
        xerrup.append(plot_data[2][k] - x[k])
        yerrlo.append(y[k] - plot_data[4][k])
        yerrup.append(plot_data[5][k] - y[k])
    plt.errorbar(x, y, yerr=[yerrlo, yerrup], xerr=[xerrlo, xerrup], fmt=dmarker, markersize=2, lw=0.5, elinewidth=1, capsize=1, label="data")
    plt.legend()
    if (logplot):
        plt.plot(x, log_fit_fn(x, fit_p[0][0], fit_p[0][1]), 'r', label=fit_label(fit_p, pcrit))
        plt.legend()
        plt.xlabel('Log(Critical Amplitude - Initial Amplitude)')
        plt.ylabel('Log(Ricci Scalar Maximum Value)')
    else:
        plt.xlabel('Critical Amplitude - Initial Amplitude')
        plt.ylabel('Ricci Scalar Maximum Value')
    titletext = 'Max Ricci Scalar vs. Initial Amplitude for $\Delta^{2}$ = ' + str(dsq) + ', dr = ' + str(0.1 / float(resn)) + ' ($r_{max} = 50.0$)'
    plt.title(titletext)
    return plot_data, fit_p

def pwr_fit_fn(xdata, gam, c):
    return c * np.power(xdata, -2*gam)

def log_fit_fn(xdata, gam, c):
    return -2*gam*xdata + c

def pcrit_fit_fn(xdata, pcrit, gam, c):
    return c * np.power(pcrit - xdata, -2*gam)

def pcrit_log_fit_fn(xdata, pcrit, gam, c):
    return c - 2*gam*np.log(pcrit - xdata)

def fit_ricci_data(fnames, amps, pcrit, emax=[0.02,0.8,3.1], resn=8, pwr_fit=False, pcrit_fit=False, guess=None):
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

def get_titletext(pcrit, offstart, offend):
    return "pcrit = " + str(pcrit) + "\nremoved: smallest " + str(offstart) + " & largest " + str(offend)

def logplot_vs_fit(titletext, fnames, amps, pcrit, emax=[0.02,0.8,3.1], resn=8, guess=[0.37, -3]):
    plot_data = get_ricci_max_log_crit_data(fnames, amps, pcrit, emax, resn)
    x = plot_data[0]
    y = plot_data[3]
    popt, pcov = curve_fit(log_fit_fn, x, y, p0=guess)
    plt.plot(x, y, 'bo', markersize=2)
    plt.plot(x, log_fit_fn(x, popt[0], popt[1]), 'r')
    plt.title(titletext)
    return 0

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

nms = ["_01224", "_01225", "_012255", "_01226", \
       "_012265", "_012268", "_012270", "_012272", "_012274", "_012275", \
       "_012276", "_012277", "_0122775", "_012278", "_0122785"]

ams = [0.01224, 0.01225, 0.012255, 0.01226, \
       0.012265, 0.012268, 0.012270, 0.012272, 0.012274, 0.012275, \
       0.012276, 0.012277, 0.0122775, 0.012278, 0.0122785]

append_amp_file_names(nms, ams, 22790, 22826, dstep=1)

nms += ["_012283", "_0122835", "_012284", "_0122845", \
        "_012285", "_0122855", "_012286"]
       
ams += [0.012283, 0.0122835, 0.012284, 0.0122845, 0.012285, 0.0122855, 0.012286]

iend = len(nms)
if (iend != len(ams)):
    print("ERROR: names and amps differ")


n83_84 = nms[11:(iend - 14)]
a83_84 = ams[11:(iend - 14)]
# *******************************
n82_83 = nms[11:(iend - 15)]
a82_83 = ams[11:(iend - 15)]

n3_4 = nms[12:(iend - 16)]
a3_4 = ams[12:(iend - 16)]
# *******************************
n2_3 = nms[12:(iend - 17)]
a2_3 = ams[12:(iend - 17)]

n3_4t = nms[13:(iend - 18)]
a3_4t = ams[13:(iend - 18)]
# *********************************
n2_3t = nms[13:(iend - 19)]
a2_3t = ams[13:(iend - 19)]



#plt.figure()
#logd90, p90 = ricci_max_errplot(n90, a90)
#plt.figure()
#logd90t, p90t = ricci_max_errplot(n90t, a90t)
# plt.figure()
# logplot_vs_fit(get_titletext(0.012283, 11, 14), n83_84, a83_84, 0.012283)
# plt.figure()
# logplot_vs_fit(get_titletext(0.012284, 11, 14), n83_84, a83_84, 0.012284)
# plt.figure()
# logplot_vs_fit(get_titletext(0.012283, 12, 16), n3_4, a3_4, 0.012283)
# plt.figure()
# logplot_vs_fit(get_titletext(0.012284, 12, 16), n3_4, a3_4, 0.012284)
# plt.figure()
# logplot_vs_fit(get_titletext(0.012283, 13, 18), n3_4t, a3_4t, 0.012283)
# plt.figure()
# logplot_vs_fit(get_titletext(0.012284, 13, 18), n3_4t, a3_4t, 0.012284)


# plt.figure()
# logplot_vs_fit(get_titletext(0.012282, 11, 15), n82_83, a82_83, 0.012282)
# plt.figure()
# logplot_vs_fit(get_titletext(0.012283, 11, 15), n82_83, a82_83, 0.012283)
# plt.figure()
# logplot_vs_fit(get_titletext(0.012282, 12, 17), n2_3, a2_3, 0.012282)
# plt.figure()
# logplot_vs_fit(get_titletext(0.012283, 12, 17), n2_3, a2_3, 0.012283)
# plt.figure()
# logplot_vs_fit(get_titletext(0.012282, 13, 19), n2_3t, a2_3t, 0.012282)
# plt.figure()
# logplot_vs_fit(get_titletext(0.012283, 13, 19), n2_3t, a2_3t, 0.012283)


#plt.show(block=False)


pcrit_arr = 0.0000001 * np.array(range(122880, 122901, 1))
more_results = []
more_tight = []
for i in range(25):
    #print("\n" + str(i) + "\n") 
    for j in range(25):
    #    print("******** off front = " + str(i) + "********")
    #    print("******** off back = " + str(j) + "********\n")
        for pcr in pcrit_arr:
            if (pcr > max(ams[i:(iend - j)])):
                params = fit_ricci_data(nms[i:(iend - j)], ams[i:(iend - j)], \
                                        pcr, emax=[0,0,0], resn=8, pwr_fit=False, \
                                        pcrit_fit=False, guess=[0.37, -3])
                gam = params[0][0]
                if ((gam < 0.38) and (gam > 0.36)):
    #                print("gamma_fit = " + str(gam))
    #                print("FOR pcrit = " + str(pcr) + "\n")
                    more_results.append([gam, pcr, i, j, params])
                    if ((i > 10) and (j > 5)):
                        more_tight.append([gam, pcr, i, j, params])



# n89 = nms[5:]
# a89 = ams[5:]

# n895 = nms[6:]
# a895 = ams[6:]

# n915 = nms[11:]
# a915 = ams[11:]

# n90 = nms[18:(iend - 8)]
# a90 = ams[18:(iend - 8)]

# n90t = n90[9:(iend - 9)]
# a90t = a90[9:(iend - 9)]

# n88_89 = nms[11:(iend - 9)]
# a88_89 = ams[11:(iend - 9)]

# n87_88 = nms[11:(iend - 10)]
# a87_88 = ams[11:(iend - 10)]

# n86_87 = nms[11:(iend - 11)]
# a86_87 = ams[11:(iend - 11)]

# n85_86 = nms[11:(iend - 12)]
# a85_86 = ams[11:(iend - 12)]

# n84_85 = nms[11:(iend - 13)]
# a84_85 = ams[11:(iend - 13)]

# n83_84 = nms[11:(iend - 14)]
# a83_84 = ams[11:(iend - 14)]
# # *******************************
# n82_83 = nms[11:(iend - 15)]
# a82_83 = ams[11:(iend - 15)]

# n81_82 = nms[11:(iend - 16)]
# a81_82 = ams[11:(iend - 16)]

# n8_9 = nms[12:(iend - 11)]
# a8_9 = ams[12:(iend - 11)]

# n7_8 = nms[12:(iend - 12)]
# a7_8 = ams[12:(iend - 12)]

# n6_7 = nms[12:(iend - 13)]
# a6_7 = ams[12:(iend - 13)]

# n5_6 = nms[12:(iend - 14)]
# a5_6 = ams[12:(iend - 14)]

# n4_5 = nms[12:(iend - 15)]
# a4_5 = ams[12:(iend - 15)]
# # *******************************
# n3_4 = nms[12:(iend - 16)]
# a3_4 = ams[12:(iend - 16)]

# n2_3 = nms[12:(iend - 17)]
# a2_3 = ams[12:(iend - 17)]

# n1_2 = nms[12:(iend - 18)]
# a1_2 = ams[12:(iend - 18)]

# n8_9t = nms[13:(iend - 13)]
# a8_9t = ams[13:(iend - 13)]

# n7_8t = nms[13:(iend - 14)]
# a7_8t = ams[13:(iend - 14)]

# n6_7t = nms[13:(iend - 15)]
# a6_7t = ams[13:(iend - 15)]

# n5_6t = nms[13:(iend - 16)]
# a5_6t = ams[13:(iend - 16)]

# n4_5t = nms[13:(iend - 17)]
# a4_5t = ams[13:(iend - 17)]

# n3_4t = nms[13:(iend - 18)]
# a3_4t = ams[13:(iend - 18)]
# # *********************************
# n2_3t = nms[13:(iend - 19)]
# a2_3t = ams[13:(iend - 19)]

# n8_9tt = nms[14:(iend - 16)]
# a8_9tt = ams[14:(iend - 16)]

# n7_8_9 = nms[14:(iend - 17)]
# a7_8_9 = ams[14:(iend - 17)]

# n6_7_8 = nms[14:(iend - 18)]
# a6_7_8 = ams[14:(iend - 18)]

# n6_7tt = nms[14:(iend - 19)]
# a6_7tt = ams[14:(iend - 19)]






