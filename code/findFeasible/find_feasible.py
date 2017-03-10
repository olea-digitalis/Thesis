from parseCount import *
from gridObj import *
from parseNodes import *
from fragCount import *


#this list contains all of the wnt entries in ncipid-elements.txt:
wnt_ls = ['pid_154','pid_168','pid_185','pid_206','pid_246','pid_250','pid_302','pid_319','pid_330','pid_334','pid_350','pid_361','pid_368','pid_400','pid_454','pid_467','pid_492','pid_532','pid_791','pid_82477']



out = "feasible_out.txt"


node_file = "signaling-hypergraph-hyperedges.txt"
hedge_file = "signaling-hypergraph-hyperedges.txt"
elements_file = "ncipid-elements.txt"

def find_entry(id, elements):
    #given a pid_number, searches through ncipid-elements.txt to find the entry associated with the id
    #returns that entry
    with open(elements,'r') as f:
        for line in f:
            if line.split('\t')[0] == id:
                return line
        s = id + " not found.\n"
        return s
                
def find_feasible(source, hedgefile, nodefile, elements, outfile):
    #reads in node and hedge files and creates a standard graph representation
    #runs BFS (as implemented in my hypergraph-stats scripts from summer 2016) starting from the given source node
    #writes down the entry from ncipid-elements.txt that corresponds with each of the reachable nodes in an output file
    nodes = parse_nodes(nodefile)
    hedges = parse_hedges(hedgefile)
    populate_nodes(nodes, hedges, True)
    node_names = [n.name for n in nodes]
    if source not in node_names:
        print("Warning! " + source + " not found in node file!")
        return
    feas_ls = frag_BFS(nodes,source, ignore_direction = False, out_ls = True)
    with open(outfile,'a') as o:
        for id in feas_ls:
            o.write(find_entry(id, elements))



def iterate_find_feasible(ls, hedgefile, nodefile, elements, outfile):
    #wrapper function for find_feasible
    #runs find_feasible on a super source represented as a list of nodes
    #gives an itemized list of which node can reach which targets
    with open(outfile,'w') as o:
        o.write("#####################################\n")
        o.write("FEASIBLE TARGETS FOR " + str(ls) + "\n")
        o.write("#####################################\n")
        o.write("\n\n")
        for source in ls:
            o.write("##################\n")
            o.write("TARGETS FOR " + source + ":\n")
            o.write(find_entry(source, elements))
            o.write("##################\n")
            find_feasible(source, hedgefile, nodefile, elements, outfile)
            o.write("\n\n")

iterate_find_feasible(wnt_ls,hedge_file,node_file, elements_file, out)