#this script prints a human readable summary of the output from the ICHA

#GLOBALS
#keeps track of the index files for lookups
#modifies which file we're investigating (outprefix)

outprefix = "wnt3a_wnt5a_bcat_nucleus"
hyperedge_data = "signaling-hypergraph-hyperedges.txt"

node_file = "signaling-hypergraph-hypernodes.txt"
hedge_file = "signaling-hypergraph-hyperedges.txt"
elements_file = "ncipid-elements.txt"
reactions_file = "ncipid-reactions.txt"
subpathways_file = "ncipid-subpathways.txt"
controls_file = "ncipid-controls.txt"
complexes_file = "ncipid-complexes.txt"

ref_ls = [elements_file, complexes_file, subpathways_file, controls_file, reactions_file]



#accession functions
#get the filename for various ilp output files

def order_file(outprefix,solnum=1):
    return outprefix+'-'+str(solnum)+'.order-variables'

def edge_file(outprefix,solnum=1):
    return outprefix+'-'+str(solnum)+'.edge-variables'

def node_file(outprefix,solnum=1):
    return outprefix+'-'+str(solnum)+'.node-variables'

def results_file(outprefix):
    return outprefix+'_icILP_results.txt'



#helper functions

def iterate_replace(old, new, ls):
    c = ls[:]
    for i in range(len(c)):
        c[i] = c[i].replace(old, new)
    return c

def find_entry(id, reference_files):
    #given a pid_number, searches through ncipid-elements.txt to find the entry associated with the id
    #returns that entry
    for filename in reference_files:
        with open(filename,'r') as f:
            for line in f:
                if line.split('\t')[0] == id:
                    return line.strip() + '\t' + filename + '\n'
    s = id + " not found.\n"
    return s




#ic_ILP output interpreting functions

def scroll_through_order(outprefix, reference_files=None):
    #goes through the order variables and puts them in order
    #this lets you see when each node in the path was visited
    with open(order_file(outprefix), 'r') as order:
        save = []
        order.readline()
        order.readline()
        for line in order:
            ls = line.split('\t')
            if float(ls[1]) != 0:
                save.append([ls[0],float(ls[1])])
        save.sort(key=lambda x:x[1])
        for e in save:
            print(e)
            if reference_files:
                print(find_entry(e[0][2:],reference_files))


def st_edgevars(outprefix):
    #goes through the edge variables file and finds each hyperedge
    #that was used
    #stores each hedge as a list: [tail, head] with delimiter ;
    #returns a list of all the hedges
    with open(edge_file(outprefix),'r') as edge:
        save = []
        edge.readline()
        edge.readline()
        for line in edge:
            ls = line.split('\t')
            if int(ls[3]) == 1:
                save.append(iterate_replace('|',';',ls[1:-1]))

        return save


def lookup_hyperedges(outprefix, hyperedge_file, reference_files = None):
    #goes through the hyperedge tail and head pids to find the associated
    #reaction pids. looks up the reaction pids.
    #prints a summary for each hyperedge in the hyperpath
    in_hyperpath = st_edgevars(outprefix)
    with open(hyperedge_file,'r') as hf:
        hf.readline()
        for e in in_hyperpath:
            for line in hf:
                ls = line.strip().split('\t') #tab delineated list of lookup line components, ls[0] is tail and ls[1] is head
                if e[0] == ls[0] and e[1] == ls[1]:
                    print('####################### NEW HYPEREDGE ##################')
                    print(ls)
                    print(e)
                    print(line)
                    if reference_files:
                        print("REACTION LOOKUP: "+find_entry(line.strip().split('\t')[-1],reference_files)) #looks up the reaction pid
                        print
                        print("TAIL NODE LOOKUP: ")
                        for node in e[0].split(';'):
                            print('\t'+find_entry(node,reference_files))
                        print("---------------------")

                        print("HEAD NODE LOOKUP: ")
                        for node in e[1].split(';'):
                            print('\t'+find_entry(node,reference_files))
                        print("---------------------")
                elif e[0] in ls[0] and e[1] == ls[1]:
                    print("####################### NEW HYPEREDGE #################")
                    print("THIS EDGE IS A CHEAT")
                    print(e)
                    print(line)
                    if reference_files:
                        print("REACTION LOOKUP: "+find_entry(line.strip().split('\t')[-1],reference_files))
                        print
                        print("TAIL NODE LOOKUP: ")
                        print('\t'+find_entry(e[0],reference_files))
                        print("CHEATED NODES LOOKUP: ")
                        for node in ls[0].split(';'):
                            if node != e[0]:
                                print('\t'+find_entry(node,reference_files))
                        print('--------------------')

                        print("HEAD NODE LOOKUP: ")
                        for node in e[1].split(';'):
                            print('\t'+find_entry(node,reference_files))
                        print('--------------------')
            hf.seek(1)



#run the following functions to get a summary

print('ORDER VARIABLES:')
scroll_through_order(outprefix, ref_ls)
print
print
print

print('HEDGE LOOKUP:')
lookup_hyperedges(outprefix, hyperedge_data, ref_ls)
