#!/usr/bin/python

## Script to run the Cheating Hyperpath Algorithm
## Modified from Anna's original Shortest s-t B-Hyperpath Algorithm script
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

#path to CPLEX package on VT local server:
sys.path.append('/data/annaritz/tools/CPLEX/2014/CPLEX_Studio126/cplex/python/x86-64_linux/')

#from ilp import *  #Anna's original ilp script
import mod_ilp     #cheating hyperpath ilp script

## import HALP library.  https://github.com/Murali-group/halp
from halp.directed_hypergraph import DirectedHypergraph
from halp.utilities import directed_graph_transformations,directed_statistics
from halp.algorithms import directed_paths





################################
################################
################################
def main(args):
    opts = parseOptions(args) #parses command line arguments into an object
    H = read_hypergraph(opts) #parses the infile into a hypergraph object
    if opts.verbose:
        print '%d nodes and %d hyperedges' % (directed_statistics.number_of_nodes(H),directed_statistics.number_of_hyperedges(H))
        for hedge in H.hyperedge_id_iterator():
            print '%s: %s --> %s' % (hedge,opts.nodedelim.join(H.get_hyperedge_tail(hedge)),opts.nodedelim.join(H.get_hyperedge_head(hedge)))


    #runs the ilp
    cH = iterate_cheat_ILP(H,opts.source,opts.target,opts.outprefix,opts.numsols,opts.subopt,opts.verbose)

    #summary linking hedge id's to nodes
    if opts.verbose:
        for hedge in cH.hyperedge_id_iterator():
            print(str(hedge) + ': ' + str(opts.nodedelim.join(cH.get_hyperedge_tail(hedge))) + ' --> ' + str(opts.nodedelim.join(cH.get_hyperedge_head(hedge))))


    return
#################################
#################################
#################################


##########################
#PARSING AND INITIALIZING#
##########################

def parseOptions(args):
    """
    Parses command line arguments.
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

def read_hypergraph(opts, debug=False):
    H = DirectedHypergraph()
    H.read(opts.hedgefile,opts.nodedelim,opts.coldelim)

    nodes = H.get_node_set()
    if opts.source not in nodes and debug==False:
        sys.exit('ERROR: source "%s" is not in node set. Exiting.' % (opts.source))
    if opts.target not in nodes and debug==False:
        sys.exit('ERROR: target "%s" is not in node set. Exiting.' % (opts.target))
    return H

        
        
        


################################
#cheating hpath setup functions#
################################

def get_restrictive_hedges(H,n):
    #returns the set of restrictive hyperedges for a node n as a list

    #details:
    #the restrictive hedges of a node is a subset of the forward star
    #any hedge in the forward star of n with |tail_set| > 1 is a restrictive hedge
    #(i.e. any hedge whose tail contains more than n itself is a restrictive hyperedge)
    #the forward star of a node is the set of hyperedges such that the node is in the tail of each hyperedge in that set.

    hedges = H.get_forward_star(n)
    r_hedges = []
    for h in hedges:
        if len(H.get_hyperedge_tail(h)) > 1:
            r_hedges.append(h)
    return r_hedges


def add_cheat_hedges(H):
    #adds cheat hedges to every node in the hypergraph H, labeled with attribute: 'cheat' = True
    #returns a list of the added cheat hedge ID's

    #details:
    #for every restrictive hedge (see function: get_restrictive_hedges() ) of each node, we add a cheat hedge
    #cheat hedges model the restrictive hedges but they remove the restriction
    #that is, the tail of a cheat hedge contains only the given node, and the head contains the head set of the given restrictive hedge
    #each cheat edge keeps track of which hedge it is modeled after via the 'original_hedge' attribute.
    
    
    added_cheats = []
    for n in H.get_node_set():
        r_hedges = get_restrictive_hedges(H,n)
        
        for h in r_hedges:
            if not H.has_hyperedge([n],H.get_hyperedge_head(h)):
                #if there is an existing hedge with the same head and tail of the candidate cheat hedge, we don't need to make it. (fixes an issue with halp)
                c = H.add_hyperedge([n], H.get_hyperedge_head(h), attr_dict={'cheat': True, 'original_hedge': h})
                #c is a hyperedge whose tail is only n and whose head is the same as that of the restrictive hyperedge it's modeled after
                added_cheats.append(c)
    return added_cheats










###############################
#MY ILP FOR CHEATING HYPERPATH#
###############################


def compute_cheating_hyperpath(cH,cheatset,k,source,target,outprefix,numsols,subopt,verbose):
    #constructs and executes ilp files to find the cheating hyperpath for a given cheat parameter k
    
    ## Check whether there is a feasible solution 
    bconnected,ignore1,ignore2,ignore3 = directed_paths.b_visit(cH,source)
    if target not in bconnected:
        print 'TARGET NOT IN INPUT.'
        print(bconnected)
        return None, None

    ## get set of nodes (for parsing output)
    nodeset = cH.get_node_set()

    ## consider only induced subhypergraph to reduce the solution space
    cH = cH.get_induced_subhypergraph(bconnected)
    
    ## Build the ILP.
    lpfile = '%s.lp' % (outprefix)
    mod_ilp.make_cheatinghyperpath_ilp(cH,k,cheatset,source,target,lpfile)

    ## Run the ILP
    numsols,numoptobjective,allvars,times = mod_ilp.solveCheatILP(cH,nodeset,lpfile,outprefix,numsols,subopt,verbose)

    if numsols == 0:
        print 'INFEASIBLE SOLUTION.'
        return None, None

    ## return variables (first solution indexed at 0)
    print '%d solutions returned (%d optimal)' % (numsols,numoptobjective)

    return allvars, times





def iterate_cheat_ILP(H_orig,source,target,outprefix,numsols,subopt,verbose, debug=False):
    #finds the cheating hyperpath between source and target for all interesting values of the cheat parameter k
    #constructs a report of the solutions as a txt file:
    #outprefix + "_icILP_results.txt"
    #returns the cheat transform of the given hypergraph c(H_orig)
    
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
        if allvars:
            allvars = allvars[0] #for some reason allvars is a list containing only a dict
        else:
            #if the cheating hyperpath for some k is infeasible, then the cheating hyperpath for all lower k's is also infeasible, so we can stop iterating when we get to an infeasible k.
            return None

        if debug:
            #helpful printout for debugging, does not show by default
            print("#####################################################\n"*3)
            print(allvars)
            print(type(allvars))
            print(allvars.keys())
            print(cheatset)
            print_all_hedges(cH)
            print("#####################################################\n"*3)
        next_k = icILP_recorder(cH,k,allvars,cheatset,outprefix)
        k = next_k

    return cH







def icILP_recorder(cH,k,allvars,cheatset,outprefix):
    #this function will record a summary of each iteration of the icILP in a txt file:
    #outprefix + "_icILP_results.txt"
    #it will also compute and return the k for the next iteration of the icILP
    
    #you probably don't need to know any of the details below

    
    #details:
    #the next value of k is determined by observing the number of cheats used in the previous iteration's path
    #k is set to that number of cheats (such that the next iteration is allowed one fewer cheat than was used in the previous /solution/)
    #emphasis on /solution/ since k is an upper bound. If solution uses fewer than k cheats and we simply count down by one each time, every k value between the present k and that number of cheats will redundantly find the same solution.
    #setting k in this way skips over the redundant solutions.
    #the reason the next k is computed in this particular function is because this is the only cheating-hyperpath-relevant function that communicates with the outputs of the ilp.
    #having the recorder function compute the next k is weird, but writing a separate function to go find the results of the recorder function seemed more cumbersome

    results_file = outprefix + "_icILP_results.txt"
    
    if allvars == None:
        with open(results_file,'a') as rf:
            rf.write('\n\nRESULTS: k=' + str(k) + '\n')
            rf.write('NO SOLUTION')
        return -1
    
    path_nodes = []
    path_hedges = []
    path_cheats = []
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
        rf.write('Hyperedge Index:\n')
        rf.write('hedge_id\ttail_nodes\thead_nodes\tis_cheat\n')
        sort_ids(path_hedges)
        for e in path_hedges:
            d = cH.get_hyperedge_attributes(e)
            rf.write(e + '\t' + set_string(d['tail']) + '\t' + set_string(d['head']) + '\t' + str(safe_lookup(d, 'cheat')) + '\n')

        rf.write('\n\n')
        rf.write('Cheat Hyperedge Index:\n')
        rf.write('hedge_id\toriginal_hedge\toriginal_tail\toriginal_head\n')

        sort_ids(path_cheats)
        for e in path_cheats:
            d = cH.get_hyperedge_attributes(e)
            d_orig = cH.get_hyperedge_attributes(d['original_hedge'])
            rf.write(e+'\t'+str(d['original_hedge'])+'\t'+set_string(d_orig['tail'])+'\t'+set_string(d_orig['head'])+'\n')
    
    return next_k




###helper functions for icILP_recorder()###

def sort_ids(ls):
    #sorts a list of hedge ids numerically
    f = lambda x: int(x[1:])
    ls.sort(key=f)

def set_string(s,delim=';'):
    out = ""
    for element in s:
        out += str(element) + delim
    out = out[:-len(delim)]
    return out

def safe_lookup(d, key):
    if key in d:
        return True
    else:
        return False

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
    #helper function for parsing ILP output files.
    #determines whether an alpha variable belongs to an edge.
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











###various diagnostic functions that are irrelevant to the cheating hyperpath algorithm###



def print_hedge(H,hedge):
    print hedge, "\ttail:", H.get_hyperedge_tail(hedge),"\thead:", H.get_hyperedge_head(hedge), "\tattributes:", H.get_hyperedge_attributes(hedge)

def print_all_hedges(H):
    hedges = list(H.get_hyperedge_id_set())
    for hedge in sorted(hedges):
        print_hedge(H,hedge)



def find_b_fragments(H,dictionary=False):
    #counts the number of nodes B-connected to each node in H
    #if dictionary=True, returns a dictionary specifying the connections
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
    #helper function, basically makes a histogram in dictionary form
    d = {}
    for i in ls:
        if i in d:
            d[i] += 1
        else:
            d[i] = 1
    return d


def val_len(d, key):
    #accession helper function
    return len(d[key])

def dict_max(d, obj_function):
    #finds the key of a dictionary whose value is maximized according to the objective function
    max_key = d.keys()[0]
    max_val = obj_function(d,d.keys()[0])
    for k in d.keys():
        if obj_function(d,k) >= max_val:
            max_key = k
            max_val = obj_function(d,k)
    return max_key




def compile_shortest_hpaths(H,root,outprefix,subopt,verbose):
    #this times finding the shortest hyperpath between a root node and each
    #node in its bconnected set. 
    bconnected, ignore1, ignore2, ignore3 = directed_paths.b_visit(H,root)
    time_ls = []
    i = 0
    t = len(bconnected)

    with open(outprefix + '_times.csv','w') as time_file:
        time_file.write("root: " + str(root) + '\n' + 'size of bconnected set: ' + str(t) + '\n\n')
        time_file.write("Target\tTime (seconds)\n")
    
    for n in bconnected:
        print("Round "+str(i)+'/'+str(t)+'\t'+str(n))
        allvars, duration = compute_shortest_b_hyperpath(H,root,n,outprefix,1,subopt,verbose)
        if duration != None:
            time_ls.append(duration)
        with open (outprefix + '_times.csv','a') as time_file:
            time_file.write(str(n) + '\t' + str(duration) + '\n')
        i+=1

    with open(outprefix + '_times.csv','a') as time_file:
        time_file.write('\n\n###########################################')
        time_file.write('\nSUMMARY\n')
        time_file.write('###########################################\n')
        time_file.write(get_time_stats(time_ls,string_form=True) + '\n')

    return

def get_time_stats(compile_times,string_form=False):
    #returns the min time, max time, median, and mean
    all_times = []
    for t_ls in compile_times:
        for t in t_ls:
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
        

    







##############################################
#ANNA'S ORIGINAL SHORTEST HYPERPATH ALGORITHM#
##############################################

#for posterity/reference

"""
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

"""





###########
if __name__=='__main__':
    main(sys.argv)
