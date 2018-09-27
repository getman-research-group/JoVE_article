#!/usr/bin/env python
# This code use lammps dump file to calculate HB lifetime distribution
# based on TJ's HB calculation code
# Xiaohong Zhang
# Getman Research Group
# 06/25/2018

import string,sys
import math
import os
from inspect import currentframe, getframeinfo
import re
import fileinput
import numpy as np
from operator import itemgetter
import scipy as sp
import scipy.stats
import warnings



actualStart = 0         # if dump file not start at 0
timestep = 10           # in the original data, the minimum time span between two frames
N_first = 0             # the beginning time which computation starts, absolute time 
N_last  = 5000000       # the end time which computation finishes
nevery = 1              # print every this many snapshots 
width = 200000          # the time span which computation uses, test 50ps lifetime
frameLine = 203         # the number of lines between two "ITEM:TIMESTEP", i.e. atom# + 9
warnings.simplefilter("ignore")

############## sepcify atom type from data file ########################################                                                                                     ################                                                                      
rC = [2]                                                                                                           ################                                                       
rHS = [3,4]                                                                                                            ################                                                      
rOS = [4]                                                                                                              ################                                                    
rHW = [5]                                                                                                                ################                                                  
rOW = [6]                                                                                                                ################      
#################################################################################################
################################################################################################################

def dist(X,Y): # square of distance
    x1 = X[1]
    x2 = Y[1]
    y1 = X[2]
    y2 = Y[2]
    z1 = X[3]
    z2 = Y[3]
    return 2*(x1-x2)*(y1-y2)*a*b*cosgamma + 2*(y1-y2)*(z1-z2)*b*c*cosalpha + 2*(x1-x2)*(z1-z2)*a*c*cosbeta + (x1-x2)*(x1-x2)*a*a + (y1-y2)*(y1-y2)*b*b + (z1-z2)*(z1-z2)*c*c

def dipo(X,Y):  # input is dipole vector rather than point coords
    x1 = X[1]*(float(words[1])-float(words[0])-float(words[2])) + X[2]*float(words[2])
    x2 = Y[1]*(float(words[1])-float(words[0])-float(words[2])) + Y[2]*float(words[2])
    y1 = X[2]*(float(words[4])-float(words[3]))
    y2 = Y[2]*(float(words[4])-float(words[3]))
    z1 = X[3]*(float(words[7])-float(words[6]))
    z2 = Y[3]*(float(words[7])-float(words[6]))
    if (x2*x2+y2*y2+z2*z2 == 0 or x1*x1+y1*y1+z1*z1 == 0):
        return 0
    else:
        return (x1*x2+y1*y2+z1*z2)/math.sqrt((x2*x2+y2*y2+z2*z2)*(x1*x1+y1*y1+z1*z1))

def init(items):
    temp = []
    x = []
    tmp = []
    OW = []
    h1 = h2 = o1 = o2 = [None] * 100000   # at least 8 times larger than one box h2o atoms.
    for j in range(len(items)):
        temp = [float(n) for n in items[j]]
        x.append(temp)

    for i in range(len(x)):
        x[i][0] = int(x[i][0])
        x[i][1] = int(x[i][1])
        x[i][2] = x[i][2] #+ 1 - x[0][2] 
        x[i][3] = x[i][3] #- x[0][3]
        x[i][4] = x[i][4] #- x[0][4]
    h1 = findNonperiodic (rHS,x)   # item# and x, y, z for each h1 or h2 or o1 or o2  
    h2 = find (rHW,x)   # item# and x, y, z for each h1 or h2 or o1 or o2
    o1 = findNonperiodic (rOS,x)   # item# and x, y, z for each h1 or h2 or o1 or o2
    o2 = find (rOW,x)   # item# and x, y, z for each h1 or h2 or o1 or o2
    h1 = [k for k in h1 if k != None]  #h1
    h2 = [k for k in h2 if k != None]  #h2
    o1 = [k for k in o1 if k != None]  #o1
    o2 = [k for k in o2 if k != None]  #o2


    for i in range(len(o1)):
       for j in range(len(o2)):
            do1o2 = dist(o1[i] , o2[j])
            if do1o2 <= float(12.25):
               for l in range(len(h2)):
                   do2h2 = dist(o2[j],h2[l])
                   if do2h2<=1.2:
                        do1h2 = dist(o1[i] , h2[l])
                        A1 = ang( do1o2 , do2h2, do1h2)
                        A2 = ang( do2h2 , do1h2, do1o2)
                        # if A2>=120:  
                        if A1<=30 and A2>=120 and do1h2<=6.25:
                            for n in range(len(h2)):
                                if dist(o2[j], h2[n])<=1.2: 
                                    tmp.append((o2[j][4],h2[n][4]))
                                      
               for l in range(len(h1)):
                   do1h1 = dist(o1[i],h1[l])
                   if do1h1<=1.2:
                        do2h1 = dist(o2[j] , h1[l])
                        # print do1h1,do1o2,do2h1
                        A1 = ang( do1o2 , do1h1, do2h1 )
                        A2 = ang( do1h1 ,do2h1, do1o2 )
                        # if  A2>=120:
                        if  A1<=30 and A2>=120 and do2h1<=6.25:
                            for m in range(len(h2)):
                                if dist(o2[j],h2[m])<=1.2:  # output HBed Ow and corresponed Hw. At least two pairs, Ow-Hw1 AND Ow-Hw2
                                    tmp.append((o2[j][4],h2[m][4])) # o2[j][4] is j line with item# ([4], already reordered
                               
    tmp = list(set(tmp))   # convert 2 lines(2 pairs in 1 h2o mlc: Ow-Hw1 AND Ow-Hw2) to 1 lines: Ow Hw1 Hw2, thus deletes repeating Ow
    for p in range(len(tmp)):
        for q in range(p+1,len(tmp)):  # set the repeated value to -1
            if(tmp[p][0] == tmp[q][0] and tmp[p][1] == tmp[q][1]):
                tmp[q][0]= -1
        for t in range(p+1, len(tmp)):
            if(tmp[t][0]>0 and tmp[p][0] == tmp[t][0]):
                OW.append((tmp[p][0],tmp[p][1],tmp[t][1]))
                #print tmp[p][0],tmp[p][1],tmp[t][1]
    return OW

def ang(d1,d2,d3):# Note: square of the distances
    return math.degrees(math.acos((d1 + d2 - d3)/(2.0 * math.sqrt(d1) * math.sqrt(d2))))

def findNonperiodic(y,c):
    x = grab = [None] *100000
    for i in range(len(c)):
        grab = ((c[i][1],c[i][2],c[i][3],c[i][4],c[i][0]))  # read the whole line and reorder each string
        for j in range (len(y)):
            if grab[0] == y[j]:
                x[i]=grab   # x[i] is the re-ordered whole line with item# and x, y, z
    return (x)

#below is periodic boundary 
def find(y,c):
    x=grab= grabxp = grabxn = grabyp = grabyn = grabxpyp = grabxnyn = grabxnyp = grabxpyn = [None] * 100000
    for i in range(len(c)):
        grab = ((c[i][1],c[i][2],c[i][3],c[i][4],c[i][0]))      
        grabxp = ((c[i][1],c[i][2]+1,c[i][3],c[i][4],c[i][0]))
        grabxn = ((c[i][1],c[i][2]-1,c[i][3],c[i][4],c[i][0]))
        grabyp = ((c[i][1],c[i][2],c[i][3]+1,c[i][4],c[i][0]))
        grabyn = ((c[i][1],c[i][2],c[i][3]-1,c[i][4],c[i][0]))
        grabxpyp = ((c[i][1],c[i][2]+1,c[i][3]+1,c[i][4],c[i][0]))
        grabxpyn = ((c[i][1],c[i][2]+1,c[i][3]-1,c[i][4],c[i][0]))
        grabxnyp = ((c[i][1],c[i][2]-1,c[i][3]+1,c[i][4],c[i][0]))
        grabxnyn = ((c[i][1],c[i][2]-1,c[i][3]-1,c[i][4],c[i][0]))
        for j in range(len(y)):
            if grab[0] == y[j]:
                x[9*i+0]=grab
                x[9*i+1]=grabxp
                x[9*i+2]=grabxn
                x[9*i+3]=grabyp
                x[9*i+4]=grabyn
                x[9*i+5]=grabxpyp
                x[9*i+6]=grabxpyn
                x[9*i+7]=grabxnyp
                x[9*i+8]=grabxnyn
    return (x)

def coor(words):  # calculate a,b,c,cosalpha,cosbeta,cosgamma
    data = [None] * 6
    xlo_bound = float(words[0])
    xhi_bound = float(words[1])
    xy =  float(words[2])
    ylo_bound =  float(words[3])
    yhi_bound =  float(words[4])
    xz =  float(words[5])
    zlo =  float(words[6])
    zhi =  float(words[7])
    yz =  float(words[8])
    xlo = xlo_bound - min(0.0,xy,xz,xy+xz)
    xhi = xhi_bound - max(0.0,xy,xz,xy+xz)
    ylo = ylo_bound - min(0.0,yz)
    yhi = yhi_bound - max(0.0,yz)
    lx = xhi-xlo
    ly = yhi-ylo
    lz = zhi-zlo
    data[0] = lx   #a
    data[1] = math.sqrt(ly*ly+xy*xy)   #b
    data[2] = math.sqrt(lz*lz+xz*xz+yz*yz)  #c
    data[3] = (xy*xz+ly*yz)/data[1]/data[2]   #cosalpha
    data[4] = xz/data[2]    #cosbeta
    data[5] = xy/data[1]   #cosgamma
    return (data)


def compute_avg(avgVal, last_avgVal):
    output = {}
    avg = []
    for x in range(len(avgVal)):
        if avgVal[x][0] in output.keys():
            output[avgVal[x][0]][0] += avgVal[x][1]
            output[avgVal[x][0]][1] += avgVal[x][2]
        else:
            output[avgVal[x][0]] = [avgVal[x][1], avgVal[x][2]]
    for y in range(len(last_avgVal)):
        if last_avgVal[y][0] in output.keys():
            output[last_avgVal[y][0]][0] += last_avgVal[y][1]
            output[last_avgVal[y][0]][1] += last_avgVal[y][2]
        else:
            output[last_avgVal[y][0]] = [last_avgVal[y][1],last_avgVal[y][2]]
    for key in output.keys():
        if output[key][1] == 0:
            avg.append((key, 0))
        else:
            avg.append((key, output[key][0]/output[key][1]))
    avg.sort(key=itemgetter(0))
    return avg

# reference: 
# https://stackoverflow.com/questions/15033511/compute-a-confidence-interval-from-sample-data
def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
#    if n < 5:        
#        m = str()
#        h = str()
#    else:
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    m = format(m, '.15f')
    h = format(h, '.15f')
    return m, h   

def get_dist_parameter(myFrame):
    words = [None]*9
    for l in range(len(myFrame[5])):
        words[l] = myFrame[5][l]
    if (len(myFrame[5]) == 2):
        words[2] = "0.0"
    for m in range(len(myFrame[6])):
        words[3+m] = myFrame[6][m]
    if (len(myFrame[6]) == 2):
        words[5] = "0.0"
    for n in range(len(myFrame[7])):
        words[6+n] = myFrame[7][n]
    if (len(myFrame[7]) == 2):
        words[8] = "0.0"
    result = coor(words) 
    return result
                                            
def compute_hb_dict(one_frame):  

    frameTime = int(one_frame[1][0])
    print frameTime
    N = int(one_frame[3][0])
    items = [ None] * N
    one_hb_dict = []
    xv = yv = zv = 0
    for t in range(N):
        items[t] = one_frame[t+9]   # saved coordinate for each frame. Already delete the head 9 lines, only coordinates, e.g. 1 1 0.01 1.23 2.31
    one_hb_dict = init(items)   ## hydrogen bonded water list, each line is Ow Hw1 Hw2, e.g. 41 76 77s 
    #print(one_hb_dict)
    # o_list is Oxygen in each frame
    o_list = [item[0] for item in one_hb_dict]  # change [(30, 59, 56), (31, 58, 55)...] to [30, 31 ...]
    for key in hb_dict.keys():  # first comprasion, remove the key when hb breaks
            if key not in o_list:
                if hb_dict[key] > 0: # if value > 0
                    hb_distribution.append(hb_dict[key])
                hb_dict.pop(key)  # delete this element
    if len(o_list) > 0: # if not empty
        for i in range(len(o_list)):  # second comprasion, add new or update existing hb
            if o_list[i] in hb_dict.keys():
                hb_dict[o_list[i]] += 1
            else:
                hb_dict[o_list[i]] = 0
    #print(hb_dict)
    #print(hb_distribution)
    #print("===========")

############################################################################################### 
if __name__ == '__main__':
    frame = []
    cell_param = []
    hb_dict = {}    # global var, b/c updating each 
    hb_distribution = []  # global var, b/c already used append in compute_hb_dict() function
    rootDir = '.'
    line_counter = 0

    for lammpfile in os.listdir(rootDir):
        if lammpfile.endswith("lammpstrj"):
            myfile = lammpfile
    print(myfile)

    for line in fileinput.input("%s" %(str(myfile))):
        if (fileinput.lineno() >= frameLine*(N_first-actualStart)/timestep+1 and fileinput.lineno() <= frameLine*((N_last-actualStart)/timestep+1)):
            if ((fileinput.lineno()-1)/frameLine%nevery == 0):
                wholeline = line.split()   # split the line on runs of whitespace
                words = [s for s in wholeline if s != None]
                frame.append(words)
                line_counter += 1
            if line_counter == frameLine:  # already read one frame
                if len(cell_param) == 0:  # only execute this when compute the first frame, only execute one time
                    cell_param = get_dist_parameter(frame)
                    a = cell_param[0]
                    b = cell_param[1]
                    c = cell_param[2]
                    cosalpha = cell_param[3]
                    cosbeta = cell_param[4]
                    cosgamma = cell_param[5]
                compute_hb_dict(frame)
                frame = []
                line_counter = 0
    fileinput.close()

    m_distribution, h_distribution = mean_confidence_interval(hb_distribution)
    #print(m_distribution, h_distribution)
    # keep the output finally
    np.savetxt("distribution_HB_lifetime.dat", hb_distribution, fmt="%d")
    f = open('avg_HB_lifetime.dat', 'w')
    f.write('   '.join((m_distribution, h_distribution)))
#########################################################################################################################################################



