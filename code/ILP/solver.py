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
import random
import time
import copy

from ilp import *  ## import functions from ilp.py
import mod_ilp

## import HALP library.  https://github.com/Murali-group/halp
from halp.directed_hypergraph import DirectedHypergraph
from halp.utilities import directed_graph_transformations,directed_statistics
from halp.algorithms import directed_paths

################################
def main(args):
    opts = parseOptions(args)
    H = read_hypergraph(opts,mod=True)
    if opts.verbose:
        print '%d nodes and %d hyperedges' % (directed_statistics.number_of_nodes(H),directed_statistics.number_of_hyperedges(H))
        """
        for hedge in H.hyperedge_id_iterator():
            print '%s: %s --> %s' % (hedge,opts.nodedelim.join(H.get_hyperedge_tail(hedge)),opts.nodedelim.join(H.get_hyperedge_head(hedge)))
        """
    """
    H_nodes = H.get_node_set()
    pairs = []
    for i in range(0,3000):
        pairs.append(random.sample(H_nodes,2))
    """

    """
    #testing cheat hedges code
    print(H.get_node_set())
    print(H.get_hyperedge_id_set())
    print_all_hedges(H)

    add_cheat_hedges(H)
    print("##############################\nAdded cheat hedges\n#####################\n")
    
    print(H.get_node_set())
    print(H.get_hyperedge_id_set())
    print_all_hedges(H)

    """

    #find_b_fragments(H)
    
    #allvars, times = compile_shortest_hpaths(H,opts.source,opts.target,opts.outprefix,opts.numsols,opts.subopt, opts.verbose,path_directions = pairs)
    
    #allvars = compute_shortest_b_hyperpath(H,opts.source,opts.target,\
    #    opts.outprefix,opts.numsols,opts.subopt,opts.verbose)
    """
    cH = copy.deepcopy(H)
    cheatset = add_cheat_hedges(cH) #add_cheating_hedges() modifies the given object and returns a list of all cheating hedge id's
    k = 1

    

    allvars = compute_cheating_hyperpath(cH,cheatset,k,opts.source,opts.target,\
        opts.outprefix,opts.numsols,opts.subopt,opts.verbose)
    """
    iterate_cheat_ILP(H,opts.source,opts.target,opts.outprefix,opts.numsols,opts.subopt,opts.verbose)

    """
    print("\n\nHere's the hgraph:")
    print_all_hedges(H)
    print_all_hedges(cH)
    print('\n')
    print("Printing allvars:")
    print(allvars)
    print(type(allvars[0][0]))
    """
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

def read_hypergraph(opts, mod=False):
    H = DirectedHypergraph()
    H.read(opts.hedgefile,opts.nodedelim,opts.coldelim)

    nodes = H.get_node_set()
    if opts.source not in nodes and mod==False:
        sys.exit('ERROR: source "%s" is not in node set. Exiting.' % (opts.source))
    if opts.target not in nodes and mod==False:
        sys.exit('ERROR: target "%s" is not in node set. Exiting.' % (opts.target))
    return H


def get_frag_pairs(H, n):
    #UNFINISHED
    i = 0
    nodes = H.get_node_set()
    while i < n:
        current = random.choice(nodes)
        bconnected, ignore1,ignore2,ignore3 = directed_paths.b_visit(H,current)
        if len(bconnected) == 1:
            pass
        else:
            pass
        



def print_hedge(H,hedge):
    print hedge, "\ttail:", H.get_hyperedge_tail(hedge),"\thead:", H.get_hyperedge_head(hedge), "\tattributes:", H.get_hyperedge_attributes(hedge)

def print_all_hedges(H):
    hedges = list(H.get_hyperedge_id_set())
    for hedge in sorted(hedges):
        print_hedge(H,hedge)



################################
#cheating hpath setup functions#
################################

def get_restrictive_hedges(H,n):
    #the restrictive hedges of a node is a subset of the forward star
    #any hedge in the forward star of n with |tail_set| > 1 (i.e. contains more nodes than just n) is a restrictive hedge
    #the forward star of a node is the set of hyperedges such that the node is in the tail of each hyperedge in that set.
    hedges = H.get_forward_star(n)
    r_hedges = []
    for h in hedges:
        if len(H.get_hyperedge_tail(h)) > 1:
            r_hedges.append(h)
    return r_hedges


def add_cheat_hedges(H):
    #adds cheat hedges to every node in H, labeled with attribute: 'cheat' = True
    #for every restrictive hedge (see function: get_restrictive_hedges() ) of each node, we add a cheat hedge
    #cheat hedges model the restrictive hedges but they remove the restriction
    #that is, the tail of a cheat hedge contains only the given node, and the head contains the head set of the given restrictive hedge
    #each cheat edge keeps track of which hedge it is modeled after via the 'original_hedge' attribute.
    #returns a list of the added cheat hedge ID's

    
    added_cheats = []
    for n in H.get_node_set():
        r_hedges = get_restrictive_hedges(H,n)
        for h in r_hedges:
            if not H.has_hyperedge([n],H.get_hyperedge_head(h)): #this conditional is necessary because the halp library won't let you have two hyperedges with the same head and tail. If there is already with the head and tail of the candidate cheat edge, we simply do not make it.
                c = H.add_hyperedge([n], H.get_hyperedge_head(h), attr_dict={'cheat': True, 'original_hedge': h})
                added_cheats.append(c)
    return added_cheats







def find_b_fragments(H,dictionary=False):
    if dictionary:
        d = {}
    ls = []
    for n in H.get_node_set():
        bconnected,ignore1,ignore2,ignore3 = directed_paths.b_visit(H,n)
        ls.append(len(bconnected))
        if dictionary:
            d[n]=bconnected


    ls.sort(reverse=True)

    if dictionary:
        return d
    print(count(ls))
    return ls


def count(ls):
    d = {}
    for i in ls:
        if i in d:
            d[i] += 1
        else:
            d[i] = 1
    return d
    


def compute_shortest_b_hyperpath(H_orig,source,target,outprefix,numsols,subopt,verbose):
    
    ## Get induced sub-hypergraph on b-connected nodes 
    bconnected,ignore1,ignore2,ignore3 = directed_paths.b_visit(H_orig,source)
    if target not in bconnected:
        print 'TARGET NOT IN INPUT.'
        return None, None
    H = H_orig.get_induced_subhypergraph(bconnected)

    if verbose:
        #print '%d nodes are B-connected to source "%s"' % (len(bconnected),source)
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
        return None, None
  
    ## return variables (first solution indexed at 0)
    print '%d solutions returned (%d optimal)' % (numsols,numoptobjective)
    return allvars, times






###############################
#MY ILP FOR CHEATING HYPERPATH#
###############################


def compute_cheating_hyperpath(cH,cheatset,k,source,target,outprefix,numsols,subopt,verbose):
    ## Get induced sub-hypergraph on b-connected nodes 
    bconnected,ignore1,ignore2,ignore3 = directed_paths.b_visit(cH,source)
    if target not in bconnected:
        print 'TARGET NOT IN INPUT.'
        print(bconnected)
        return None, None
    #print("Printing cH now:\n")
    #print_all_hedges(cH)
    
    #H = cH.get_induced_subhypergraph(bconnected)
    #print("\n\nPrinting H now:\n")
    #print_all_hedges(H)

    
    ## Build the ILP.
    lpfile = '%s.lp' % (outprefix)
    mod_ilp.make_cheatinghyperpath_ilp(cH,k,cheatset,source,target,lpfile)

    ## get set of nodes (for parsing output)
    nodeset = cH.get_node_set()

    ## Run the ILP
    numsols,numoptobjective,allvars,times = mod_ilp.solveCheatILP(cH,nodeset,lpfile,outprefix,numsols,subopt,verbose)

    if numsols == 0:
        print 'INFEASIBLE SOLUTION.'
        return None, None

    ## return variables (first solution indexed at 0)
    print '%d solutions returned (%d optimal)' % (numsols,numoptobjective)

    return allvars, times





def iterate_cheat_ILP(H_orig,source,target,outprefix,numsols,subopt,verbose):
    results_file = outprefix + "_icILP_results.txt"
    with open(results_file,'w') as rf:
        rf.write("-------------------------------------\n")
        rf.write("iterate_cheat_ILP results\n")
        rf.write("-------------------------------------\n\n\n")
    cH = copy.deepcopy(H_orig)
    cheatset = add_cheat_hedges(cH) #add_cheating_hedges() modifies the given object and returns a list of all cheating hedge id's
    k = len(cheatset)+1

    while k > 0:
        allvars, times = compute_cheating_hyperpath(cH,cheatset,k,source,target,outprefix,numsols,subopt,verbose)
        allvars = allvars[0] #for some reason this is a list containing only a dict
        print("#####################################################")
        print("#####################################################")
        print("#####################################################")
        print(allvars)
        print(type(allvars))
        print(allvars.keys())
        print(cheatset)
        print_all_hedges(cH)
        print("#####################################################")
        print("#####################################################")
        print("#####################################################")
        next_k = icILP_recorder(k,allvars,cheatset,outprefix)
        k = next_k

    return







def icILP_recorder(k,allvars,cheatset,outprefix):
    #this function will record a summary of each iteration of the icILP in a txt file
    #it will also compute and return the k for the next iteration of the icILP

    results_file = outprefix + "_icILP_results.txt"
    
    if allvars == None:
        with open(results_file,'a') as rf:
            rf.write('\n\nRESULTS: k=' + str(k) + '\n')
            rf.write('NO SOLUTION')
        return -1
    
    path_nodes = []
    path_hedges = []
    path_cheats = []
    #for some reason allvars is returned as a list containing only a dictionary
    #so that's why it's allvars[0] in the next line
    for v in allvars.keys():
        if var_is_alpha(v) and allvars[v]==1:
            if var_is_edge(v):
                path_hedges.append(get_var_id(v))
                if var_is_cheat_edge(v, cheatset):
                    path_cheats.append(get_var_id(v))
            else:
                path_nodes.append(get_var_id(v))

    num_nodes = len(path_nodes)
    num_hedges = len(path_hedges)
    num_cheats = len(path_cheats)
    next_k = num_cheats

    with open(results_file,'a') as rf:
        rf.write('\n\nRESULTS: k=' + str(k) + '\n')
        rf.write('# nodes:\t# hedges:\t# cheats:\tnodes+hedges:\n')
        rf.write(str(num_nodes) + '\t' + str(num_hedges) + '\t' + str(num_cheats) + '\t' + str(num_nodes + num_hedges) + '\n')
        rf.write('\n\nNodes in Path: ' + str(path_nodes) + '\n\n')
        rf.write('Hyperedges in Path: ' + str(path_hedges) + '\n\n')
        rf.write('Cheat Hyperedges in Path: ' + str(path_cheats) + '\n\n')
    
    return next_k


def write_ref_sheet(H, cH, cheatset,allvars=None):
    #this function will construct a reference sheet so you can quickly figure out what a returned path should actually look like.
    pass




def get_var_id(v):
    #gets the id from a variable in the ILP output
    return v[2:]

def var_is_alpha(v):
    #determines whether a variable from the ILP output is an alpha variable
    if v[0] == 'a':
        return True
    else:
        return False

def var_is_edge(v):
    #this helper function determines whether an alpha variable from the ILP
    #output belongs to an edge. On the toy graphs, edges are labeled e<number>
    #but it's concievable that naming conventions for somebody else's graphs
    #might complicate matters. So for now this is a placeholder, to make fixing
    #things easier if such a problem arises.
    var_id = get_var_id(v)
    if var_id[0] == 'e':
        try:
            int(var_id[1])
            return True
        except:
            return False
    else:
        return False

def var_is_cheat_edge(v, cheatset):
    #given that v is an alpha variable of an edge, determines whether that
    #edge is a cheat edge or not
    var_id = get_var_id(v)
    if var_id in cheatset:
        return True
    else:
        return False







def compile_shortest_hpaths(H,source,target,outprefix,numsols,subopt,verbose,path_directions=None):
    if path_directions == None:
        allvars, times = compute_shortest_b_hyperpath(H,source,target,outprefix,numsols,subopt,verbose)
        return allvars
    else:
        compile_times = []
        i = 0
        for pair in path_directions:
            with open (outprefix + '_working_record.txt','a') as wr:
                wr.write(str(pair) + time.asctime() + '\n')
            allvars, times = compute_shortest_b_hyperpath(H,str(pair[0]),str(pair[1]),outprefix,numsols,subopt,verbose)

            if allvars:
                with open (outprefix + '_working_record.txt','a') as wr:
                    wr.write("Solution found!\n")
                compile_times.append([pair,times])
                i += 1

            if i == 29:
                break
        with open (outprefix + '_times.csv','a') as time_file:
            time_file.write(get_time_stats(compile_times,string_form=True) + '\n')
            
            for thing in compile_times:
                time_file.write(str(thing) + '\n')

        return allvars, times

def get_time_stats(compile_times,string_form=False):
    #returns the min time, max time, median, and mean
    all_times = []
    for thing in compile_times:
        for t in thing[1]:
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
