from parseCount import *
from gridObj import *
from parseNodes import *
from fragCount import *


#this list contains all of the wnt entries in ncipid-elements.txt:
wnt_ls = ['pid_154','pid_168','pid_185','pid_206','pid_246','pid_250','pid_302','pid_319','pid_330','pid_334','pid_350','pid_361','pid_368','pid_400','pid_454','pid_467','pid_492','pid_532','pid_791','pid_82477']

notch_ls = ['pid_132413','pid_132418','pid_132530','pid_93907','pid_93910','pid_94101','pid_94106','pid_94109','pid_94137','pid_94148','pid_94154','pid_94344','pid_94435','pid_94551']

#Hedgehog
Hh_ls = "pid_24145 pid_24217 pid_24221 pid_24225 pid_24226 pid_24286 pid_24310 pid_24312 pid_24322 pid_24333 pid_36242 pid_87931".split(' ')

#retinoic acid receptor
RAR_ls = "pid_62912 pid_62915 pid_62929 pid_62948 pid_62952 pid_62984 pid_63052 pid_62984 pid_63066 pid_63139 pid_63211 pid_63227 pid_63282 pid_63368".split(' ')

dopamine_ls = "pid_67037 pid_67039".split(' ')

#breast cancer gene
BRCA1_ls = "pid_40897 pid_7404 pid_90826 pid_92007 pid_92209 pid_92256".split(' ')

EGF_ls = "pid_14480 pid_21562 pid_60514".split(' ')

out = "feasible_out.txt"


node_file = "signaling-hypergraph-hypernodes.txt"
hedge_file = "signaling-hypergraph-hyperedges.txt"
elements_file = "ncipid-elements.txt"
reactions_file = "ncipid-reactions.txt"
subpathways_file = "ncipid-subpathways.txt"
controls_file = "ncipid-controls.txt"
complexes_file = "ncipid-complexes.txt"

ref_ls = [elements_file, complexes_file, subpathways_file, controls_file, reactions_file]

def find_entry(id, reference_files):
    #given a pid_number, searches through ncipid-elements.txt to find the entry associated with the id
    #returns that entry
    for filename in reference_files:
        with open(filename,'r') as f:
            for line in f:
                if line.split('\t')[0] == id:
                    return line
    s = id + " not found.\n"
    return s
                
def find_feasible(source, hedgefile, nodefile, outfile, reference_files):
    #reads in node and hedge files and creates a standard graph representation
    #runs BFS (as implemented in my hypergraph-stats scripts from summer 2016) starting from the given source node
    #writes down the entry from ncipid-elements.txt that corresponds with each of the reachable nodes in an output file
    nodes = parse_nodes(nodefile)
    hedges = parse_hedges(hedgefile)
    populate_nodes(nodes, hedges, True)

    source_node = None
    for n in nodes:
        if n.name == source:
            source_node = n
    if not source_node:
        print("Warning! " + source + " not found in node file!")
        return
    
    feas_ls = frag_BFS(nodes,source_node, ignore_direction = False, out_ls = True)
    found_entries = ""
    with open(outfile,'a') as o:
        for id in feas_ls:
            entry = find_entry(id, reference_files)
            if 'not found.' not in entry:
                found_entries += entry
            o.write(entry)
    return found_entries



def iterate_find_feasible(ls, hedgefile, nodefile, outfile, reference_files):
    #wrapper function for find_feasible
    #runs find_feasible on a super source represented as a list of nodes
    #gives an itemized list of which node can reach which targets
    with open(outfile,'w') as o:
        o.write("#####################################\n")
        o.write("FEASIBLE TARGETS FOR " + str(ls) + "\n")
        o.write("#####################################\n")
        o.write("\n\n")

    summary_str = "\n\nListing all successfully identified outputs:\n"
    for source in ls:
        with open(outfile,'a') as o:
            o.write("##################\n")
            o.write("TARGETS FOR " + source + ":\n")
            o.write(find_entry(source,reference_files))
            o.write("##################\n")
        summary_str += find_feasible(source, hedgefile, nodefile, outfile, reference_files)
        with open(outfile,'a') as o:
            o.write("\n\n")
    with open(outfile,'a') as o:
        o.write(summary_str)

iterate_find_feasible(wnt_ls,hedge_file,node_file, out, ref_ls)
