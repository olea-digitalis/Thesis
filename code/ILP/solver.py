#!/usr/bin/python

## Stand-alone script to run the Shortest B-Hyperpaths problem.
## Modified from VT: /data/annaritz/hypergraph/2014-05-15-acm-bcb-hypergraphs/master-script.py

import os
import os.path
import subprocess
import sys
from optparse import OptionParser
import numpy as np
import glob


from ilp import *  ## import functions from ilp.py


## import HALP library.  https://github.com/Murali-group/halp
from halp.directed_hypergraph import DirectedHypergraph
from halp.utilities import directed_graph_transformations,directed_statistics
from halp.algorithms import directed_paths

################################
def main(args):
    opts = parseOptions(args)
    H = read_hypergraph(opts)
    if opts.verbose:
        print '%d nodes and %d hyperedges' % (directed_statistics.number_of_nodes(H),directed_statistics.number_of_hyperedges(H))
        for hedge in H.hyperedge_id_iterator():
            print '%s: %s --> %s' % (hedge,opts.nodedelim.join(H.get_hyperedge_tail(hedge)),opts.nodedelim.join(H.get_hyperedge_head(hedge)))
    pairs = []
    for i in range(1,61):
        if i % 2 == 0:
            pairs.append(tuple([i,i+1]))

    allvars, times = compile_shortest_hpaths(H,opts.source,opts.target,opts.outprefix,opts.numsols,opts.subopt, opts.verbose,path_directions = pairs)
    
    #allvars = compute_shortest_b_hyperpath(H,opts.source,opts.target,\
        #opts.outprefix,opts.numsols,opts.subopt,opts.verbose)
    return 
   

################################
def parseOptions(args):
    """
    Parses arguments passed to function.
    """

    desc = 'Computes Shortest B-Hyperpaths between a source and a target in a directed hypergraph consisting of nodes.\n\npython acm_bcb_shortest_b_hyperpaths.py [options] HEDGEFILE SOURCE TARGET\n\n\tHEDGEFILE: two column delimited file of hyperedges in the form of <TAIL> <HEAD>, where \n\t  <TAIL> and <HEAD> may be delimited sets of nodes. Note: Underscore "_" is not allowed\n\t  in the node names.\n\tSOURCE: Source node.\n\tTARGET: Target node.'
    parser = OptionParser(usage=desc)

    parser.add_option('','--coldelim',type='string',metavar='STR',default='\t',
        help='Column delimiter in HEDGEFILE. Default="\\t".')
    parser.add_option('','--nodedelim',type='string',metavar='STR',default='|',
        help='Node delimiter in HEDGEFILE. Default="|".')
    parser.add_option('','--numsols',type='int',default=1,metavar='INT',
        help='Return the top k solutions. Default=1.')
    parser.add_option('','--subopt',action='store_true',
        help='Allow suboptimal solutions (when k>1). Default=False.')
    parser.add_option('-v','--verbose',action='store_true',
        help='Print additional information to the screen. Default=False')
    parser.add_option('-o','--outprefix',type='string',metavar='STR',default='out',
        help='Output file prefix. Default is "out".')

    opts,args = parser.parse_args()

    if len(args) != 3:
        parser.print_help()
        sys.exit('\nERROR: three arguments required HEDGEFILE SOURCE TARGET.\n')

    ## add required arguments to opts object:
    opts.hedgefile = args[0]
    opts.source = args[1]
    opts.target = args[2]

    if opts.verbose:
        print 'Inputs:',opts

    return opts

def read_hypergraph(opts):
    H = DirectedHypergraph()
    H.read(opts.hedgefile,opts.nodedelim,opts.coldelim)

    nodes = H.get_node_set()
    if opts.source not in nodes:
        sys.exit('ERROR: source "%s" is not in node set. Exiting.' % (opts.source))
    if opts.target not in nodes:
        sys.exit('ERROR: target "%s" is not in node set. Exiting.' % (opts.target))
    return H

def compute_shortest_b_hyperpath(H,source,target,outprefix,numsols,subopt,verbose):
    
    ## Get induced sub-hypergraph on b-connected nodes 
    bconnected,ignore1,ignore2,ignore3 = directed_paths.b_visit(H,source)
    H = H.get_induced_subhypergraph(bconnected)

    if verbose:
        print '%d nodes are B-connected to source "%s"' % (len(bconnected),source)
        print 'Hypergraph of B-connected nodes now has %d nodes and %d hyperedges.' % (directed_statistics.number_of_nodes(H),directed_statistics.number_of_hyperedges(H))

    ## Build the ILP.
    lpfile = '%s.lp' % (outprefix)
    make_shortesthyperpath_ilp(H,source,target,lpfile)

    ## get set of nodes (for parsing output)
    nodeset = H.get_node_set()

    ## Run the ILP
    numsols,numoptobjective,allvars,times = solveILP(H,nodeset,lpfile,outprefix,numsols,subopt,verbose)
    if verbose:
        print 'Done.\n'
        
    if numsols == 0:
        print 'INFEASIBLE SOLUTION.'
        return None
  
    ## return variables (first solution indexed at 0)
    print '%d solutions returned (%d optimal)' % (numsols,numoptobjective)
    return allvars, times


def compile_shortest_hpaths(H,source,target,outprefix,numsols,subopt,verbose,path_directions=None):
    if path_directions == None:
        allvars, times = compute_shortest_b_hyperpath(H,source,target,outprefix,numsols,subopt,verbose)
        return allvars
    else:
        compile_times = []
        for pair in path_directions:
            allvars, times = compute_shortest_b_hyperpath(H,str(pair[0]),str(pair[1]),outprefix,numsols,subopt,verbose)
            compile_times.append(times)
        with open (outprefix + '_times.csv','a') as time_file:
            time_file.write(get_time_stats(compile_times,string_form=True))
            time_file.write('\n' + str(compile_times) + '\n\n')

        return allvars, times

def get_time_stats(compile_times,string_form=False):
    #returns the min time, max time, median, and mean
    all_times = []
    for times in compile_times:
        for t in times:
            all_times.append(t)

    max_time = max(all_times)
    min_time = min(all_times)
    mean = get_mean(all_times)
    median = get_median(all_times)

    if not string_form:
        return min_time,max_time,median,mean

    else:
        s = "min: " + str(min_time) + "\t max: " + str(max_time) + "\t median: " + str(median) + "\t mean: " + str(mean)
        return s
    
def get_mean(ls):
    total = 0
    for n in ls:
        total += n

    mean = total / float(len(ls))
    return mean

def get_median(ls):
    s_ls = sorted(ls)
    if len(ls) % 2 == 1:
        median = s_ls[len(ls)/2]
    else:
        m1 = s_ls[len(ls)/2]
        m2 = s_ls[len(ls)/2 -1]
        median = (m1 + m2)/float(2)
    return median
        

    
    
###########
if __name__=='__main__':
    main(sys.argv)
