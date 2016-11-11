import random

#random hypergraph generator

def erdos_rennie_like(num_nodes, num_hedges, max_hedge_width, self_loops=False):
    """
    generates a hypergraph using a model similar to the erdos rennie model for standard graphs.
    initializes all the nodes, randomly chooses the width for the head and tail (up to max_hedge_width), randomly picks an appropriate number of nodes to go in the hedge
    self_loops is a boolean that determines whether hedges can have the same node in the head and tail
    """
    if max_hedge_width > num_nodes:
        print("Error: the requested construction is flawed! There are not enoough nodes to support the requested max_hedge_width!")
        return -1
    
    if ((max_hedge_width * 2) > num_nodes) and (self_loops == False):
        print("Error: the requested construction is flawed! There are not enough nodes to support the requested max_hedge_width without self loops!")
        return -1
    
    node_ls = [] 
    hedge_ls = [] #hedges will be encoded as ordered lists, where entry 0 is a list containing nodes in the tail and entry 1 is a list containing nodes in the head
    hedge_set = set() #to avoid repeated edges I use a set to keep track of the edges. The set will be converted back to a list when it is returned.
    
    for n in range(num_nodes): #initializes nodes
        node_ls.append(str(n))
    
    hedges_full = False
    while hedges_full == False:
        tail_width = random.randint(1,max_hedge_width)
        head_width = random.randint(1,max_hedge_width)
        new_hedge = hedge_maker(tail_width, head_width, node_ls, self_loops)
        hedge_set.add(new_hedge)
        if len(hedge_set) >= num_hedges:
            hedges_full = True
    
    for hedge in hedge_set:    #formats the hedges into lists instead of tuples containing frozensets
        tail_ls = []
        head_ls = []
        for node in hedge[0]:
            tail_ls.append(node)
        
        for node in hedge[1]:
            head_ls.append(node)
        
        formatted_hedge = [tail_ls, head_ls]
        hedge_ls.append(formatted_hedge)
        
    return node_ls, hedge_ls


def hedge_maker(tail_width, head_width, node_ls, self_loops):
    """
    helper function for erdos_rennie_like() that selects the nodes that belong in the tail and head of the currently generating hedge
    """
    tail_set = set()
    head_set = set()

    tail_full = False
    head_full = False
    
    while tail_full == False: 
        selection = random.choice(node_ls)
        tail_set.add(selection)
        if len(tail_set) == tail_width:
            tail_full = True



    
    while head_full == False:
        selection = random.choice(node_ls)

        if self_loops == True:              #if we don't want self loops, this checks to make sure the selection isn't one of the tail nodes
            head_set.add(selection)
        else: 
            if selection not in tail_set:
                head_set.add(selection)

        if len(head_set) == head_width:
            head_full = True

    frozen_head = frozenset(head_set)
    frozen_tail = frozenset(tail_set)
    
    return (frozen_tail, frozen_head) #returns as a tuple of frozensets so that it will only be added to the edge_set if it's a unique edge


    
    
    
def set_to_list(s):
    ls = []
    for item in s:
        ls.append(item)

def export(fileprefix, hedges):
    """
    takes a list of nodes and a list of hedges and exports them to a .txt file with the given prefix
    """
    with open(fileprefix + '.txt', 'w') as f:
        for h in hedges:
            s = ""
            for node in h[0]: #each node in the tail
                s += str(node) + "|"
            s = s[:-1]
            s += '\t'
            for node in h[1]: #each node in the head
                s += str(node) + "|"
            s = s[:-1]
            s += '\t'
            s += '1' + '\n'   #assigns weight for the hedge, currently always set to 1
            f.write(s)



def main():
    """
    nodes, hdiag1 = erdos_rennie_like(600, 36000, 3, False)
    print(hdiag1)
    export('diag1',hdiag1)
    nodes, hdiag2 = erdos_rennie_like(600,180000,3)
    nodes,hdiag4 = erdos_rennie_like(600,300000,2)
    nodes,hdiag5 = erdos_rennie_like(100,8333,3)
    nodes,hdiag6 = erdos_rennie_like(100,8333,4)
    nodes,hdiag7 = erdos_rennie_like(300,75000,3)
    nodes,hdiag8 = erdos_rennie_like(450,168750,3)
    export('diag2',hdiag2)
    export('diag4',hdiag4)
    export('diag5',hdiag5)
    export('diag6',hdiag6)
    export('diag7',hdiag7)
    export('diag8',hdiag8)



    
    nodes, hdiag3 = erdos_rennie_like(600,300000,3)
    export('diag3',hdiag3)
    """

    """
    nodes, hd3 = erdos_rennie_like(100,8333,5)
    export('d3',hd3)

    nodes, hd5 = erdos_rennie_like(100,8333,6)
    export('d5',hd5)

    nodes, hd6 = erdos_rennie_like(100,8333,7)
    export('d6',hd6)
    """

    """
    nodes, sparse1 = erdos_rennie_like(600, 1200, 3)
    export('sparse_diag1', sparse1)

    nodes, sparse2 = erdos_rennie_like(600, 2400, 3)
    export('sparse_diag2',sparse2)

    nodes, sparse3 = erdos_rennie_like(600, 5800, 3)
    export('sparse_diag3',sparse3)

    nodes, sparse4 = erdos_rennie_like(600,11600, 3)
    export('sparse_diag4',sparse4)

    nodes, sparse5 = erdos_rennie_like(600,23200, 3)
    export('sparse_diag5',sparse5)
    """

    nodes, size1 = erdos_rennie_like(100, 500, 3)
    nodes, size2 = erdos_rennie_like(200,1000,3)
    nodes,size3 = erdos_rennie_like(300,1500,3)
    nodes,size4 = erdos_rennie_like(400,2000,3)
    nodes,size5 = erdos_rennie_like(500,2500,3)

    export('size_diag1',size1)
    export('size_diag2',size2)
    export('size_diag3',size3)
    export('size_diag4',size4)
    export('size_diag5',size5)

if __name__ == '__main__':
    main()
