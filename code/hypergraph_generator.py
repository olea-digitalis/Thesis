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
    hedge_ls = [] #edges will be encoded as ordered lists, where entry 0 is a list containing nodes in the tail and entry 1 is a list containing nodes in the head
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

        
        
def main():
    print(erdos_rennie_like(50, 10, 3, False))
    
    
    

if __name__ == '__main__':
    main()