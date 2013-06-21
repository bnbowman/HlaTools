#! /usr/bin/env python

import os, sys, glob

import matplotlib.pyplot as plt
import numpy as np
from pylab import xlim, ylim, xlabel

def plot_score(file_name, r_count):
    score_data = []
    with open(file_name) as f:
        for l in f:
            l = l.strip().split()
            score_data.append( (int(l[0]), l[1], int(l[2]), int(l[3]), int(l[4]), int(l[5]), float(l[6]) ) )
        coord, bases, c1, c2, c3, c4, p = zip(*score_data)
        p = np.array(p)
        p2 = 1.0*np.array(c2)/np.array(c4)
        p3 = 1.0*np.array(c3)/np.array(c4)
                                                         
        c1 = np.array(c1)
        p[c1 < 10] = 0
                                                                                   
        fig = plt.figure( figsize=(10, 24) )
        i = 1
        xspan = 1000
        for x in range(0, len(coord), xspan):
            ax = fig.add_subplot( 20, 1, i )
            #ax.plot(coord[x:x+xspan], p[x:x+xspan], "g-", 
            #        coord[x:x+xspan], p2[x:x+xspan], "r-",  
            #        coord[x:x+xspan], p3[x:x+xspan], "b-")
            ax.plot(coord[x:x+xspan], p2[x:x+xspan], "r-", 
                    coord[x:x+xspan], p3[x:x+xspan], "k-")

            ylim( -0.05, 0.75)
            xlim(x, x+xspan)
            i += 1
        xlabel("position | %s (%s reads)" % (file_name, r_count) )
                
        png_file = '.'.join( file_name.split('.')[:-1] ) + '.png'
        plt.savefig( png_file )

def plot_data(wd):
    g_rn = []
    with open(wd+"/summary.txt") as sf:
        for l in sf:
            l = l.strip().split()
            if l[0] == "total":
                l[0] = "group_root"
            g_rn.append( (l[0], l[1]) )
                                                                                      
    for g, c in g_rn:
        fn  = wd+"/%s.score" % g
        plot_score(fn, c)

if __name__ == '__main__':
    plot_data( sys.argv[1] )
