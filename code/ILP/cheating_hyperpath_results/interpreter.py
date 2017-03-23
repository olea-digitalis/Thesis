#this file is for automated lookup of the cheating hyperpath ILP output


outprefix = "wnt3a_MYF5"
hyperedge_data = "signaling-hypergraph-hyperedges.txt"

node_file = "signaling-hypergraph-hypernodes.txt"
hedge_file = "signaling-hypergraph-hyperedges.txt"
elements_file = "ncipid-elements.txt"
reactions_file = "ncipid-reactions.txt"
subpathways_file = "ncipid-subpathways.txt"
controls_file = "ncipid-controls.txt"
complexes_file = "ncipid-complexes.txt"

ref_ls = [elements_file, complexes_file, subpathways_file, controls_file, reactions_file]

def order_file(outprefix,solnum=1):
    return outprefix+'-'+str(solnum)+'.order-variables'

def edge_file(outprefix,solnum=1):
    return outprefix+'-'+str(solnum)+'.edge-variables'

def node_file(outprefix,solnum=1):
    return outprefix+'-'+str(solnum)+'.node-variables'

def results_file(outprefix):
    return outprefix+'_icILP_results.txt'

def iterate_replace(old, new, ls):
    c = ls[:]
    for i in range(len(c)):
        c[i] = c[i].replace(old, new)
    return c


def scroll_through_order(outprefix):
    with open(order_file(outprefix), 'r') as order:
        save = []
        order.readline()
        order.readline()
        for line in order:
            ls = line.split('\t')
            if float(ls[1]) != 0:
                save.append([ls[0],float(ls[1])])
        save.sort(key=lambda x:x[1])
        print(save)

def st_edgevars(outprefix):
    with open(edge_file(outprefix),'r') as edge:
        save = []
        edge.readline()
        edge.readline()
        for line in edge:
            ls = line.split('\t')
            if int(ls[3]) == 1:
                save.append(iterate_replace('|',';',ls[1:-1]))
        return save

def lookup_hyperedge_pids(outprefix, hyperedge_file):
    in_hyperpath = st_edgevars(outprefix)
    with open(hyperedge_file,'r') as hf:
        hf.readline()
        for e in in_hyperpath:
            for line in hf:
                ls = line.strip().split('\t')
                if e[0] == ls[0] and e[1] == ls[1]:
                    print(ls)
                    print(e)
                    print(line)
                elif e[0] in ls[0] and e[1] == ls[1]:
                    print("THIS EDGE IS A CHEAT")
                    print(e)
                    print(line)
            hf.seek(1)



def original_hedges(edge_ls, cheat_ls):
    new = []
    for e in edge_ls:
        if e not in cheat_ls:
            new.append(e)
    return new


def hedge_lookup(hedge,outprefix, hyperedge_file, cheat=False):
    #given an edge ID from the ILP, this finds the corresponding entry in the
    #hyperedge data file and returns it.
    #since multiple hyperedges in the data file may correspond to a single
    #ILP hyperedge (due to differences in regulators, multiple candidates for
    #cheats, etc) returns a list of lines which have the same head and tail
    #(if the hedge is not a cheat) or same head and similar element in the tail
    #(if the hedge is a cheat)
    with open(edge_file(outprefix),'r') as ef:
        found = None
        ef.readline()
        ef.readline()
        for line in ef:
            if hedge in line:
                found = line.strip().split('\t')[:-1]
        if found == None:
            print("Error. Hyperedge " + hedge + " not found in ILP logfile.")
            return None

    ilp_id = found[0]
    head = found[1].replace('|',';')
    tail = found[2].replace('|',';')
    head_set = set(head.split(';'))
    tail_set = set(tail.split(';'))
    
    with open(hyperedge_file, 'r') as hf:
        hyperedges = []
        hf.readline()
        for line in hf:
            ls = line.strip().split('\t')
            c_head = set(ls[0].split(';'))
            c_tail = set(ls[1].split(';'))
            if not cheat and head_set == c_head and tail_set == c_tail:
                hyperedges.append([ilp_id,line])
            elif cheat and head in c_head and tail_set == c_tail:
                hyperedges.append([ilp_id,line])
        if len(hyperedges) == 0:
            print("Error. Hyperedge " + hedge + " was found in ILP logfile but not hyperedge data file.")
            print(head)
            print(tail)
            return None
    return hyperedges

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


def reaction_lookup(hyperedges, reference_files):
    #use this on a single list of candidates
    #returns the reference entry for a given line in the hyperedge datafile
    save = []
    for h in hyperedges:
        ls = h[1].strip().split('\t')
        pid = ls[-1]
        entry = find_entry(pid, reference_files)
        save.append([h[0],entry])
    return save

def hyperpath_lookup(all_hedges, cheat_ls, outprefix, hyperedge_file, reference_files, outfile='out.txt'):
    hedge_ls = original_hedges(all_hedges, cheat_ls)
    with open(outfile,'w') as o:
        o.write("########################\n")
        o.write("Non-cheat hedges\n")
        o.write("########################\n")
        for h in hedge_ls:
            a = hedge_lookup(h,outprefix,hyperedge_file)
            candidate_reactions = reaction_lookup(a, reference_files)
            o.write(candidate_reactions[0][0][2:] + ':\n') #ILP id for the hedge
            o.write( str(h))
            o.write(str(a))
            o.write('\n')
            for c in candidate_reactions:
                o.write(c[1])
            o.write('\n')

        o.write("\n")
        o.write("#######################\n")
        o.write("Cheat hedges\n")
        o.write("#######################\n")
        for h in cheat_ls:
            a = hedge_lookup(h,outprefix,hyperedge_file,cheat=True)
            candidate_reactions = reaction_lookup(a, reference_files)
            o.write(candidate_reactions[0][0][2:] + ':\n')
            o.write(str(a))
            o.write('\n')
            for c in candidate_reactions:
                o.write(c[1])
            o.write('\n')
    return
            





edge_ls = ['e10599', 'e7923', 'e2867', 'e2139', 'e10938','e9555']
cheat_ls = ['e10599', 'e7923', 'e10938', 'e9555']

hyperpath_lookup(edge_ls, cheat_ls, outprefix, hyperedge_data, ref_ls)

print(reaction_lookup(hedge_lookup('e2867',outprefix, hyperedge_data),ref_ls))
print
print(reaction_lookup(hedge_lookup('e2139',outprefix,hyperedge_data),ref_ls))


h = hedge_lookup('e10599',outprefix,hyperedge_data,cheat=True)
print(reaction_lookup(h,ref_ls))
