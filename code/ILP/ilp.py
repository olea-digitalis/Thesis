from xml.dom import minidom
import sys
import cplex

EPSILON=0.0001
CONSTANT=1000000

'''
def make_bshrub_ilp_2(H,sources,targets,outfile):

	out = open(outfile,'w')

	## write objective
	out = writeObjective(H,out)

	## write constraints
	out.write('Subject To\n')
	c = 0
	out,c = writeConstraint_IfHedgeThenHnode(H,out,c)

	for t in targets:
		out.write('c%d_sources_out:' % (c))
		for s in sources:
			for hedge in H.get_forward_star(s):
				out.write(' + %s' % (pi(t,hedge)))
		out.write(' = 1\n')
		c+=1

		out.write('c%d_targets_in:' % (c))
		for hedge in H.get_backward_star(t):
			out.write(' + %s' % (pi(t,hedge)))
		out.write(' = 1\n')
		c+=1

		for hedge in H.hyperedge_id_iterator():
			if hedge in H.get_backward_star(t):
				continue
			## q_e \leq \sum_{heads of e} \sum_{forward star of heads} q_e'
			## q_e - \sum_{heads of e} \sum_{forward star of heads} q_e' \leq 0
			out.write('c%d_forward_flow: %s' % (c,pi(t,hedge)))
			for v in H.get_hyperedge_head(hedge):
				for hedge2 in H.get_forward_star(v):
					out.write(' - %s' % (pi(t,hedge2)))
			out.write(' <= 0\n')
			c+=1

		for hedge in H.hyperedge_id_iterator():
			skip = False
			for s in sources:
				if hedge in H.get_forward_star(s):
					skip = True
			if skip:
				continue

			## q_e \leq \sum_{tails of e} \sum_{backward star of tails} q_e'
			## q_e - \sum_{tails of e} \sum_{backward star of heads} q_e' \leq 0
			out.write('c%d_backward_flow: %s' % (c,pi(t,hedge)))
			for v in H.get_hyperedge_tail(hedge):
				for hedge2 in H.get_backward_star(v):
					out.write(' - %s' % (pi(t,hedge2)))
			out.write(' <= 0\n')
			c+=1

		## if q_e = 1 then alpha_e = 1
		## q_e <= a_e
		for hedge in H.hyperedge_id_iterator():
			out.write('c%d_if_pi_then_alpha: %s - %s <= 0\n' % (c,pi(t,hedge),a(hedge)))
			c+=1

	## if a_v = 1 then some incident a_e must be 1.
	## a_v < \sum_{incident e} a_e
	for hnode in H.node_iterator():
		incident = H.get_forward_star(hnode).union(H.get_backward_star(hnode))
		out.write('c%d_incident_nodes: %s' % (c,a(hnode)))
		for hedge in incident:
			out.write(' - %s' % (a(hedge)))
		out.write(' <= 0\n')
		c+=1



	#out,c = writeConstraint_fixTails(targets,out,c)

	## write bounds
	out = writeBounds(H,out,targets=targets)

	out.write('End\n')
	out.close()
	print 'Wrote to %s' % (outfile)
	return

def make_bshrub_ilp_predecessors(H,sources,targets,outfile):

	out = open(outfile,'w')

	## write objective
	out = writeObjective(H,out)

	## write constraints
	out.write('Subject To\n')
	c = 0
	out,c = writeConstraint_IfHedgeThenHnode(H,out,c)
	out,c = writeConstraint_IfHedgeThenPredecessor(H,out,c)
	out,c = writeConstraint_IfPredecessorThenHedge(H,out,c)
	out,c = writeConstraint_IfHnodeThenPredecessor(H,sources,out,c)
	out,c = writeConstraint_fixTails(targets,out,c)

	## write bounds
	out = writeBounds(H,out)

	out.write('End\n')
	out.close()
	print 'Wrote to %s' % (outfile)
	return
	'''

def make_shortesthyperpath_ilp(H,source,target,outfile):
	out = open(outfile,'w')

	## write objective
	out = writeObjective(H,out)

	## write constraints
	out.write('Subject To\n')
	c = 0
	out,c = writeConstraint_IfHedgeThenIncidents(H,out,c)
	out,c = writeConstraint_IfHnodeThenBackwardStar(H,source,out,c)
	out,c = writeConstraint_FixValues(H,{target:1},out,c)
	out,c = writeConstraint_OrderVariables(H,out,c)
	out = writeBinaryBounds(H,out)
	out.write('End\n')
	out.close()

	print 'Wrote to %s' % (outfile)
	return

def writeObjective(H,out,minimize=True):
	if minimize:
		out.write('Minimize\n')
	else:
		out.write('Maximize\n')

	for hedge in H.hyperedge_id_iterator():
		out.write(' + ' + a(hedge))
	out.write('\n')
	return out

def writeConstraint_IfHedgeThenIncidents(H,out,c):
	'''
	Writes the constraint that if hyperedge e is in the 
	solution, then all hypernodes incident to e must also
	be in the solution.  That is, for $I(e) = H(e) \cup T(e)$,
	\sum_{u I(e)} a_u >= |I(e)| a_e 
	\sum_{u I(e)} a_u - |I(e)| a_e >= 0

	'''
	for hedge in H.hyperedge_id_iterator(): # for all e \in E
		incident_set = H.get_hyperedge_head(hedge).union(H.get_hyperedge_tail(hedge))
		out.write('c%d_if_hedge_then_incidents: ' % (c))
		for hnode in incident_set: # for all u \in T(e) \cup H(e)
			out.write('+ %s ' % a(hnode))
		out.write(' - %d %s >= 0\n' % (len(incident_set),a(hedge)))
		c+=1

	return out,c

def writeConstraint_IfHnodeThenBackwardStar(H,s,out,c):
	'''
	Writes the constraint that if hypernode u is in the solution,
	then there must be at least one hyperedge in the backwards star
	that is also in the solution. This holds for all hypernodes
	except for s.
	\sum_{e \in BS(v)} a_e >= a_v
	\sum_{e \in BS(v)} a_e - a_v >= 0
	'''
	for hnode in H.node_iterator():
		if hnode == s:
			continue
		out.write('c%d_if_hnode_then_backwardstar: ' % (c))
		for hedge in H.get_backward_star(hnode):
			out.write('+ %s ' % (a(hedge)))
		out.write('- %s >= 0\n' % (a(hnode)))
		c+=1

	return out,c

def writeConstraint_FixValues(H,tofix,out,c):
	'''
	Writes the constraint that fixes variables in the tofix dictionary.
	The variable tofix is a dictionary of {hypernode: <0 or 1>}.
	'''
	for t in tofix:
		out.write('c%d_fixed: %s = %d\n' % (c,a(t),tofix[t]))  
		c+=1 #increment constraint counter
	return out,c

def writeConstraint_OrderVariables(H,out,c):
	''' 
	Writes the constraint that the order variables
	must be smaller in the tail than the order
	variables in the head.

	Order Bounds: a_v >= o_v >= 0
	o_v >= 0 and a_v - o_v >= 0
	'''
	for hnode in H.node_iterator():
		out.write('c%d_order_bounds: %s >= 0\n' % (c,o(hnode)))
		c+=1 

		out.write('c%d_order_bounds: %s - %s >= 0\n' % (c,a(hnode),o(hnode)))
		c+=1 #increment constraint counter

	'''
	Order Constraints:
	o_u <= o_v - \epsilon + C (1-a_e) forall u,v in (T(e),H(e)) forall e \in E 
	o_u <= o_v - \epsilon + C (1-a_e) 
	o_u - o_v + C a_e <= C - \epsilon
	'''
	for hedge in H.hyperedge_id_iterator(): # for all e \in E
		for u in H.get_hyperedge_tail(hedge):
			for v in H.get_hyperedge_head(hedge):
				out.write('c%d_order_constraint: %s - %s + %d %s <= %f\n' % \
							  (c,o(u),o(v),CONSTANT,a(hedge),CONSTANT-EPSILON))
				c+=1

	return out,c

def writeBinaryBounds(H,out):
	'''
	Specify all the alpha variables as binary.
	'''
	out.write('Binary\n')
	for hnode in H.node_iterator():
		out.write(' %s\n' % (a(hnode)))
	for hedge in H.hyperedge_id_iterator(): # for all e \in E
		out.write(' %s\n' % (a(hedge)))
	return out

################################

def solveILP(H,nodeset,lpfile,outprefix,numsols,subopt=False,verbose=False):
	print '\nSolving ILP...'

	print '\n' + '-'*20 + 'Cplex Output Start' + '-'*20
	ilp = cplex.Cplex()
	ilp.read(lpfile)

	numsolsfound = 1
	numoptobjective = 0
	maxobj = None
	allvars = []
        times = []
	while(numsolsfound < numsols+1):
		## Solve ILP
		print '-'*10 + 'Looking for Solution %d' % (numsolsfound) + '-'*10 + '\n'
                start = ilp.get_dettime()
		ilp.solve()
                finish = ilp.get_dettime()
                times.append(finish-start)
		print '-'*10 + 'Done Looking for Solution %d' % (numsolsfound) + '-'*10 + '\n'
		
		if ilp.solution.pool.get_num()>0:
			print 'Solution Found.'
			
			objective = ilp.solution.pool.get_objective_value(0)
			if numsolsfound == 1:
				maxobj = objective
				numoptobjective+=1
				print 'Max Objective of %d' % (objective)
			elif objective != maxobj and not subopt:
				print 'Solution (obj=%d) does not have max objective of %d: quitting.' % (objective,maxobj)
				break
			elif objective != maxobj and subopt:
				print 'Solution (obj=%d) does not have max objective of %d.' % (objective,maxobj)
			else:
				print 'ANOTHER OPTIMAL OBJECTIVE=%d' % (objective)
				numoptobjective+=1

			# Incumbent solution is 0 in the pool:
			ilp.solution.pool.write('%s-%d.sol' % (outprefix,numsolsfound),0)

			# Get variables 
			variables = getILPSolution(H,nodeset,outprefix,numsolsfound,objective,verbose)
			allvars.append(variables)

			# Add constraint for this solution.
			# See http://www-01.ibm.com/support/docview.wss?uid=swg21399929
			#ones = [var for var in variables if variables[var] == 1]
			#zeros = [var for var in variables if variables[var] == 0]
			#ind = ones + zeros
			#val = [-1.0]*len(ones)+[1.0]*len(zeros)
			#eq = cplex.SparsePair(ind,val)
			#ilp.linear_constraints.add(lin_expr = [eq], senses = ['G'], rhs = [1-len(ones)], names = ['solution%d' % (numsolsfound)])

			# all we need is the ones.
			ones = [var for var in variables if variables[var] == 1 and var[0]=='a' and var.split('_')[1] not in nodeset]
			val = [1.0]*len(ones)
			eq = cplex.SparsePair(ones,val)
			ilp.linear_constraints.add(lin_expr = [eq], senses = ['L'], rhs = [len(ones)-1], names = ['solution%d' % (numsolsfound)])
		else:
			print 'Infeasible Solution. quitting.'
			break
		numsolsfound+=1

	print '-'*20 + 'Cplex Output End' + '-'*20 + '\n'
	print '%d solutions found' % (numsolsfound-1)
	return numsolsfound-1,numoptobjective,allvars, times



################################
def getILPSolution(H,nodeset,outprefix,num,objective,verbose):
	print '\nGetting ILP Solution for Solution # %d in Pool' % (num)

	# parse xml
	xml = minidom.parse('%s-%d.sol' % (outprefix,num))
	cplexsol = xml.getElementsByTagName('CPLEXSolution')[0]
   
	# get variables
	elements = cplexsol.getElementsByTagName('variables')[0].getElementsByTagName('variable')
	variables = {}
	for v in elements:
		variables[str(v.getAttribute('name'))] = float(v.getAttribute('value'))

	eout = open('%s-%d.edge-variables' % (outprefix,num),'w')
	eout.write('# Objective = %s\n' % (objective))
	eout.write('#name\ttail\thead\tvalue\n')
	nout = open('%s-%d.node-variables' % (outprefix,num),'w')
	nout.write('# Objective = %s\n' % (objective))
	nout.write('#name\tnode\tvalue\n')
	oout = open('%s-%d.order-variables' % (outprefix,num),'w')
	oout.write('# Objective = %s\n'% (objective))
	oout.write('#name\tvalue\n')
	if verbose:
		print 'VARIABLES:'
	
	numrounded = 0
	for v in variables:
		if v[0] == 'o': # order variables.
			oout.write('%s\t%s\n' % (v,variables[v]))
			continue

		# round variable if necessary.
		if variables[v] != 0 and variables[v] != 1:
			numrounded+=1
			variables[v] = int(variables[v]+0.5)
		else:
			variables[v] = int(variables[v])

                row = [v[0],v[2:]]
		if v[0] == 'a': # node
			if verbose and variables[v] != 0:
				print v,'-->',row[1],'=',variables[v]
				
			if row[1] in nodeset: # node
				nout.write('%s\t%s\t%d\n' % (v,row[1],variables[v]))
			else: # edge
				eout.write('%s\t%s\t%s\t%d\n' % (v,'|'.join(H.get_hyperedge_tail(row[1])),\
					'|'.join(H.get_hyperedge_head(row[1])),variables[v]))		  
	eout.close()
	nout.close()
	oout.close()
	
	print ' %d variables had to be rounded.' % (numrounded)
	print ' wrote edge variables to file %s' % ('%s-%d.edge-variables' % (outprefix,num))
	print ' wrote node variables to file %s' % ('%s-%d.node-variables' % (outprefix,num))
	print ' wrote order variables to file %s' % ('%s-%d.order-variables' % (outprefix,num))

	return variables

####################
def writeConstraint_IfHedgeThenPredecessor(H,out,c):
	'''
	Writes the constraint that if hyperege e is in the solution,
	then it must act as at least one predecessor for a hypernode
	in the head of e. That is,
	a_e \leq \sum_{u \in H(e)} \pi_{e,u} 		for all e \in E
	a_e - \sum_{u \in H(e)} \pi_{e,u} \leq 0 	for all e \in E
	'''
	for hedge in H.hyperedge_id_iterator(): # for all e \in E
		out.write('c%d_if_hedge_then_predecessor: %s' % (c,a(hedge)))
		for hnode in H.get_hyperedge_head(hedge):
			out.write(' - %s' % (pi(hnode,hedge)))
		out.write(' <= 0\n')
		c+=1

	return out,c

def writeConstraint_IfPredecessorThenHedge(H,out,c):
	'''
	Writes the constraint that if a predecessor variable is in
	the solution, then the hyperege e is in the solution.  That is,
	\pi_{e,u} \leq a_e 			for all e \in E, u \in H(e)
	\pi_{e,u} - a_e \leq 0		for all e \in E, u \in H(e)
	'''
	for hedge in H.hyperedge_id_iterator(): # for all e \in E
		for hnode in H.get_hyperedge_head(hedge): #\forall u \in H(e)
			out.write('c%d_if_predecessor_then_hedge: %s - %s <= 0\n' % (c,pi(hnode,hedge),a(hedge)))
			c+=1

	return out,c

def writeConstraint_IfHnodeThenPredecessor(H,sources,out,c):
	'''
	Writes the constraint that if hypernode u is in the solution,
	then there must be exactly one predecessor variable specified.
	This only holds true for hypernodes that are not in S.
	a_u = \sum_{e:u \in H(e)} \pi_{e,u} 		for all u \in U \setminus S
	a_u - \sum_{e:u \in H(e)} \pi_{e,u} = 0		for all u \in U \setminus S
	'''
	for hnode in H.node_iterator():
		if hnode in sources:
			continue
		out.write('c%d_if_hnode_then_predecessor: %s' % (c,a(hnode)))
		for hedge in H.get_backward_star(hnode):
			out.write(' - %s' % (pi(hnode,hedge)))
		out.write(' = 0\n')
		c+=1

	return out,c

def writeConstraint_fixTails(targets,out,c):
	'''
	Sets alpha values of targets to 1.
	a_u = 1 		for all u \in T
	'''
	for hnode in targets:
		out.write('c%d_fix_tails: %s = 1\n' % (c,a(hnode)))
		c+=1
	return out,c

def writeBounds(H,out,targets=None):
	out.write('Binary\n')
	for hnode in H.node_iterator():
		out.write(a(hnode)+'\n')
	for hedge in H.hyperedge_id_iterator():
		out.write(a(hedge)+'\n')
		if targets:
			for t in targets:
				out.write(pi(t,hedge)+'\n')
		else:
			for hnode in H.get_hyperedge_head(hedge):
				out.write(pi(hnode,hedge)+'\n')
	return out

def read_solution(solfile):

	alphas = {}
	pis = {}

	xmldoc = minidom.parse(solfile)
	varlist = xmldoc.getElementsByTagName('variable')
	for variable in varlist:
		name = variable.getAttribute('name').split('_')
		value = int(float(variable.getAttribute('value'))+0.5)

		if name[0] == 'a':
			alphas[name[1]] = value
		elif name[0] == 'pi':
			pis[name[1]] = value
		else:
			sys.exit('not tracking %s' % (name[0]))
	
	return alphas,pis

def a(name):
	return 'a_%s' % (name)

def o(name):
	return 'o_%s' % (name)

def pi(hnode,hedge):
	return 'pi_%s_%s' % (hnode,hedge)
