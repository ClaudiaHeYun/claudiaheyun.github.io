# Title: The S_n-equivariant rational homology of the tropical moduli spaces Delta_{2,n}
# Author: Claudia He Yun
# Date: 8/7/2020
# Sage version: 9.0
# Note: for n = 4, 5, 6, 7, 8. Described in Section 3.2 of [arXiv link]

# Usage: top_degree_homology_rep(n)

import copy
import time
from sympy.combinatorics import Permutation as sympyPermutation
from sympy.combinatorics.generators import symmetric

#===========================================construct theta graphs===============================================
class ThetaGraph:
    def __init__(self, markings, edges):
        # markings is a list of lists: [[first branch],[second branch],[third branch],[left,right]]
        # edges is a list of lists: [[first branch],[second branch],[third branch]]
        self.markings = list(markings)
        self.edges = list(edges)
        self.sign = 1

    def see(self):
        print(str(self.sign)+"\n"+"markings: "+str(self.markings)+"\n"+"edges: "+str(self.edges))
    
    def changeSign(self, new):
        self.sign = int(new)

    # put in standard order:
    # number of markings in the first branch >= number of markings in the second >= number of markings in the third
    # when one of left and right exists, put it on the RHS
    # when both left and right exist, put the smaller # on the LHS
    # when neither of them exists, first number of the first branch < last number of the first branch
    # first # of the first branch < first # of the second branch < first # of the third branch
    
    def isNormal(self):
        # boolean that determines if a graph is in normal form        
        boo = True
        
        # if it's not true that #markings on top >= #middle >= #bottom, false
        if len(self.markings[0]) < len(self.markings[1]) or len(self.markings[1]) < len(self.markings[2]):
            boo = False
            
        # if the first condition is true, now we consider reflection & switching
        elif self.markings[3] == [-1,-1]:
            
            if self.markings[0][0] > self.markings[0][-1]:
                    boo = False
            
            elif len(self.markings[0]) == len(self.markings[1]) == len(self.markings[2]):
                if self.markings[0][0] != 0:
                    boo = False
                if self.markings[1][0] > self.markings[2][0]:
                    boo = False
            
            elif len(self.markings[0])== len(self.markings[1]):
                if self.markings[0][0] > self.markings[1][0] and self.markings[0][-1] > self.markings[1][0]:
                    boo = False
                elif self.markings[0][0] > self.markings[1][-1] and self.markings[0][-1] > self.markings[1][-1]:
                    boo = False
            
            elif len(self.markings[1]) == len(self.markings[2]):
                if self.markings[1][0] > self.markings[2][0]:
                    boo = False
        
        elif self.markings[3][0] > -1 and self.markings[3][1]==-1:
            boo = False
            
        elif self.markings[3][1] > -1 and self.markings[3][0]==-1:
            if len(self.markings[0]) == len(self.markings[1]) == len(self.markings[2]):
                if self.markings[0][0] > self.markings[1][0] or self.markings[1][0] > self.markings[2][0]:
                    boo = False
            
            elif len(self.markings[0])== len(self.markings[1]):
                if self.markings[0][0] > self.markings[1][0]:
                    boo = False
            
            elif len(self.markings[1]) == len(self.markings[2]):
                if self.markings[1][0] > self.markings[2][0]:
                    boo = False
        else:
            if self.markings[3][0] > self.markings[3][1]:
                boo = False
            elif len(self.markings[0]) == len(self.markings[1]) == len(self.markings[2]):
                if self.markings[0][0] > self.markings[1][0] or self.markings[1][0] > self.markings[2][0]:
                    boo = False
            
            elif len(self.markings[0])== len(self.markings[1]):
                if self.markings[0][0] > self.markings[1][0]:
                    boo = False
            
            elif len(self.markings[1]) == len(self.markings[2]):
                if self.markings[1][0] > self.markings[2][0]:
                    boo = False
        
        return boo
    
    def __eq__(self, other):
        if self.markings == other.markings:
            return True
        else:
            return False

#===================================useful functions=======================================================
def deepcopy(graph):
    newMarkings = copy.deepcopy(graph.markings)
    newEdges = copy.deepcopy(graph.edges)
    newGraph = ThetaGraph(newMarkings,newEdges)
    newGraph.changeSign(graph.sign)
    return newGraph

def flatten(listOfLists):
    newList = []
    for babyList in listOfLists:
        newList = newList+babyList
    return newList

def reflect(graph):
    newGraph = deepcopy(graph)
    for marking in newGraph.markings:
        marking.reverse()
    for edge in newGraph.edges:
        edge.reverse()
        
    return newGraph

def switch(graph):
    configurations = []
    # ()
    configurations.append(graph)
    
    # (12)
    graph12 = deepcopy(graph)
    graph12.markings[0], graph12.markings[1] = graph12.markings[1], graph12.markings[0]
    graph12.edges[0], graph12.edges[1] = graph12.edges[1], graph12.edges[0]
    configurations.append(graph12)
    
    # (13)
    graph13 = deepcopy(graph)
    graph13.markings[0], graph13.markings[2] = graph13.markings[2], graph13.markings[0]
    graph13.edges[0], graph13.edges[2] = graph13.edges[2], graph13.edges[0]
    configurations.append(graph13)
    
    # (23)
    graph23 = deepcopy(graph)
    graph23.markings[1], graph23.markings[2] = graph23.markings[2], graph23.markings[1]
    graph23.edges[1], graph23.edges[2] = graph23.edges[2], graph23.edges[1]
    configurations.append(graph23)
    
    # (123)
    graph123 = deepcopy(graph)
    graph123.markings[0], graph123.markings[1], graph123.markings[2] = graph123.markings[2], graph123.markings[0], graph123.markings[1]
    graph123.edges[0], graph123.edges[1], graph123.edges[2] = graph123.edges[2], graph123.edges[0],  graph123.edges[1]
    configurations.append(graph123)
    
    # (132)
    graph132 = deepcopy(graph)
    graph132.markings[0], graph132.markings[1], graph132.markings[2] = graph132.markings[1], graph132.markings[2], graph132.markings[0]
    graph132.edges[0], graph132.edges[1], graph132.edges[2] = graph132.edges[1], graph132.edges[2],  graph132.edges[0]
    configurations.append(graph132)

    return configurations

def normalize(graph):
    if not graph.isNormal():
        symmetries = switch(graph)+switch(reflect(graph))
        for sym in symmetries:
            if sym.isNormal():
                return sym
    else:
        return graph

def contract(graph):
        # the first and last markings of each branch become the left and right
        contraction = []
        for i in range(3):
            if len(graph.markings[i])>1:
                if graph.markings[3][0] == -1:
                    # make a deep copy, contract the left edge, append to the list
                    newMarkings = copy.deepcopy(graph.markings)
                    movingMark = newMarkings[i].pop(0)
                    newMarkings[3][0] = movingMark
                    newEdges = copy.deepcopy(graph.edges)
                    removedEdge = newEdges[i].pop(0)
                    
                    # sign of the new graph, (-1)^removedEdge
                    sign = 1
                    if removedEdge % 2 == 1:
                        sign = -1

                    # adjust the numbering of the edges
                    for k in range(len(newEdges)):
                        for j in range(len(newEdges[k])):
                            if newEdges[k][j]>removedEdge:
                                newEdges[k][j] = newEdges[k][j]-1

                    newGraph = ThetaGraph(newMarkings, newEdges)
                    newGraph.changeSign(sign)
                    contraction.append(newGraph)
                
                if graph.markings[3][1] == -1:
                    # make a deep copy, contract the right edge, append to the list
                    newMarkings = copy.deepcopy(graph.markings)
                    movingMark = newMarkings[i].pop()
                    newMarkings[3][1] = movingMark
                    newEdges = copy.deepcopy(graph.edges)
                    removedEdge = newEdges[i].pop()
                    
                    # sign of the new graph, (-1)^removedEdge
                    sign = 1
                    if removedEdge % 2 == 1:
                        sign = -1

                    # adjust the numbering of the edges
                    for k in range(len(newEdges)):
                        for j in range(len(newEdges[k])):
                            if newEdges[k][j]>removedEdge:
                                newEdges[k][j] = newEdges[k][j]-1

                    newGraph = ThetaGraph(newMarkings, newEdges)
                    newGraph.changeSign(sign)
                    contraction.append(newGraph)
                    
        return contraction

# a and b are two edge sets (permutations of numbers, lists of lists).
def parity(a,b):
    pi = sympyPermutation(flatten(a))
    sigma = sympyPermutation(flatten(b))
    numberOfTranspositions = len(pi.transpositions())+len(sigma.transpositions())
    parity = 1
    if numberOfTranspositions % 2 == 1:
        parity = -1
    
    return parity

# outputs 0 if graphs a and b are not the same, +/-1 if they are, sign determined by their edge permutation
def compare(a,b):
    output = 0
    normalA = normalize(a)
    normalB = normalize(b)
    if normalA.markings == normalB.markings:
        output = parity(normalA.edges,normalB.edges)
    
#    aMarkings = a.markings
#    aEdges = a.edges
#    symmetriesOfB = switch(b)+switch(reflect(b))
#    for graph in symmetriesOfB:
#        if aMarkings == graph.markings:
#            # then we need to consider the permutation of the edges
#            output = parity(aEdges,graph.edges)    

    return output

# permute the markings of the graph with pi 
def permute(pi,graph):
    newGraph = deepcopy(graph)
    for i in range(len(newGraph.markings)):
        for j in range(len(newGraph.markings[i])):
            if newGraph.markings[i][j] != -1:
                newGraph.markings[i][j] = pi(newGraph.markings[i][j])
    
    return newGraph

# now we need to generate all the distinct theta graphs.
# a,b,c are positive integers. a = #markings on the top branch, b middle, c bottom
# boo1, boo2 are booleans. boo1 true means there is a marking on the left, boo2 true means marking on the right
def generateThetaGraph(a,b,c,boo1,boo2,pi):
    perm = pi.list()
    markings = []
    
    firstBranch = []
    for i in range(a):
        firstBranch.append(perm.pop(0))
    markings.append(firstBranch)
    
    secondBranch = []
    for i in range(b):
        secondBranch.append(perm.pop(0))
    markings.append(secondBranch)
    
    thirdBranch = []
    for i in range(c):
        thirdBranch.append(perm.pop(0))
    markings.append(thirdBranch)
    
    leftAndRight = []
    if boo1:
        leftAndRight.append(perm.pop(0))
    else:
        leftAndRight.append(-1)
    if boo2:
        leftAndRight.append(perm.pop(0))
    else:
        leftAndRight.append(-1)
    markings.append(leftAndRight)
    
    edges = []
    firstRow = []
    for i in range(a+1):
        firstRow.append(i)
    edges.append(firstRow)
    
    secondRow = []
    for i in range(b+1):
        secondRow.append(a+1+i)
    edges.append(secondRow)
    
    thirdRow = []
    for i in range(c+1):
        thirdRow.append(a+1+b+1+i)
    edges.append(thirdRow)
    
    graph = ThetaGraph(markings,edges)
    
    return graph

#=================================generate basis for chain groups=============================================
def top_deg_theta_shape(n):
    if n == 4:
        return [[2,1,1]]
    if n == 5:
        return [[3,1,1],[2,2,1]]
    if n == 6:
        return [[4,1,1],[3,2,1],[2,2,2]]
    if n == 7:
        return [[5,1,1],[4,2,1],[3,3,1],[3,2,2]]
    if n == 8:
        return [[6,1,1],[5,2,1],[4,3,1],[4,2,2],[3,3,2]]

def lower_deg_theta_shape(n):
    if n == 4:
        return [[1,1,1]]
    else:
        return top_deg_theta_shape(n-1)

def generate_top_deg_basis(n):
    shapes = top_deg_theta_shape(n)
    basis = {}
    counter = 0

    for shape in shapes:
        for pi in symmetric(int(n)):
            graph = generateThetaGraph(shape[0],shape[1],shape[2],False,False,pi)
            if graph.isNormal():
                basis.update({repr(graph.markings):(counter,graph)})
                basis.update({counter:graph})
                counter = counter + 1
    return basis

def generate_lower_deg_basis(n):
    shapes = lower_deg_theta_shape(n)
    basis = {}
    counter = 0

    for shape in shapes:
        for pi in symmetric(int(n)):
            graph = generateThetaGraph(shape[0],shape[1],shape[2],False,True,pi)
            if graph.isNormal():
                basis.update({repr(graph.markings):(counter,graph)})
                basis.update({counter:graph})
                counter = counter + 1
    return basis

#=====================================generate top degree boundary map===========================================
# takes in bases for the top and lower degree chains groups, ouput top deg boundary map
def top_boundary(chainTop, chainLower):
    map = matrix(QQ,int(len(chainLower)/2),int(len(chainTop)/2),sparse=True)
    for j in range(int(len(chainTop)/2)):
        contractions = contract(chainTop[j])
        for graph in contractions:
            nGraph = normalize(graph)
            index = chainLower[repr(nGraph.markings)][0]
            entry = parity(nGraph.edges,graph.edges)
            map[index,j] = entry*graph.sign
    return map

#============================================more helper functions=================================================


# convert a sage permutation to a sympy permutation
def convert_convention_1to0(permutation):
    one_line = list(permutation)
    new_one_line = [one_line[i]-1 for i in range(len(one_line))]
    return sympyPermutation(new_one_line)

# convert a sympy permutation to a sage permutation
def convert_convention_0to1(s_permutation):
    one_line = list(s_permutation)
    new_one_line = [one_line[i]+1 for i in range(len(one_line))]
    return Permutation(new_one_line)

# generate *num* random vectors of size *size* with density *density*
def random_vec(size, num, density):
    rand_mat = random_matrix(QQ,size,num, density=density, sparse=True)
    output = [rand_mat.column(i) for i in range(num)]
    return output

# rank of a set of vectors in Q^dim(vec)
def rank(vectors):
    dim = len(vectors[0])
    total_space = VectorSpace(QQ, dim)
    span = total_space.subspace(vectors)
    return span.dimension()

# given vector *vector* in C_top with basis *dimTop* and permutation *perm*, produce *perm**vec* \in C_top
def permute_vec(perm, vec, dimTop):
    dim = len(vec)
    new_vec = vector(QQ, dim, sparse=True)
    for i in range(dim):
        if vec[i] != 0:
            permuted = permute(perm, dimTop[i])
            normalized = normalize(permuted)
            index = dimTop[repr(normalized.markings)][0]
            sign = parity(normalized.edges,dimTop[repr(normalized.markings)][1].edges)
            new_vec[index] = sign*vec[i]
    return new_vec

# extract the 11 entry of the matrix rep of the Specht module
def entry_dictionary(partition):
    spc = SymmetricGroupRepresentation(partition, 'specht')
    n = partition.size()
    mydict = {}
    for perm in Permutations(n):
        entry_11 = spc.representation_matrix(perm)[0,0]
        mydict.update({str(perm):entry_11})
    return mydict

# image of *vec* in C_top with basis *dimTop* under the map p_11 with respect to the reference representation 
# indexed by *partition*
def proj_11(dimTop, partition, vec, entry_11):
    n = partition.size()
    dim = len(vec)
    im = vector(QQ, dim, sparse=True)
    for perm in Permutations(n):
        entry = entry_11[str(perm)]
        if entry != 0:
            im += entry*permute_vec(convert_convention_1to0(perm), vec, dimTop)
    return im

# C_top has basis *dimTop* and *coefficient* many copies of irrep indexed
# this function computes how many of those copies are in the kernel of the top boundary map *dTop* restricted
# to that isotypic component
def kernel_size(dimTop, dTop, coefficient, partition, entry_11, density):
    dim = int(len(dimTop)/2)
    random_vec_set = random_vec(dim, coefficient, density)

    image_start = time.time()
    im_set = []
    loop_start = time.time()
    for i in range(coefficient):
        im_set.append(proj_11(dimTop, partition, random_vec_set[i], entry_11))
        loop_fin = time.time()
        print(f'{i+1} vector images out of {coefficient} calculated', loop_fin-loop_start)
        print('======================================================')
        loop_start = loop_fin
    image_fin = time.time()
    print('calculated all images', image_fin-image_start)
    
    ld_start = time.time()
    linear_dependency = rank(im_set)
    ld_fin = time.time()
    print('check that images form a basis', linear_dependency==coefficient)
    print(ld_fin-ld_start)
    
    rk_start = time.time()
    diff_set = [dTop*im_set[i] for i in range(coefficient)]
    rk = rank(diff_set)
    rk_fin = time.time()
    print('calculated rank of images under boundary map', rk_fin-rk_start)
    
    print('size of kernel of dTop on this component is', coefficient-rk)
    return coefficient-rk

# gives the decomposition into irreducibles for the top chain group. Start with the coefficient for [1^n].
def top_chain_group_rep(n):
    if n == 4:
        return [0,1,0,1,0]
    if n == 5:
        return [0,1,2,3,3,3,1]
    if n == 6:
        return [1,4,6,6,8,13,4,7,10,5,2]
    if n == 7:
        return [0,7,14,17,21,41,30,23,24,46,19,16,22,9,3]
    if n == 8:
        return [2,12,32,50,20,32,108,118,98,76,56,162,90,128,22,66,116,50,44,32,12,0]

def density_random_vec(n):
    if n == 4:
        return 1
    if n == 5:
        return 0.5
    if n == 6:
        return 0.1
    if n == 7:
        return 0.01
    if n == 8:
        return 0.001

def top_degree_homology_rep(n):
    # generate bases for chain groups
    basis_start = time.time()
    chainTop = generate_top_deg_basis(n)
    chainLower = generate_lower_deg_basis(n)
    basis_fin = time.time()
    print("generated bases for chain groups", basis_fin-basis_start)

    # generate top degree boundary map
    dTop_start = time.time()
    dTop = top_boundary(chainTop, chainLower)
    dTop_fin = time.time()
    print("generated the top boundary map", dTop_fin-dTop_start)
    
    partitions = Partitions(n).list()
    partitions.reverse()

    coefficients = top_chain_group_rep(n)

    density = density_random_vec(n)
    rep = []
    for i in range(len(partitions)):
        print(partitions[i])
        start = time.time()
        entry_11 = entry_dictionary(partitions[i])
        if coefficients[i] != 0:
            k = kernel_size(chainTop, dTop, coefficients[i], partitions[i], entry_11, density)
        else:
            k = 0
        rep.append(k)
        end = time.time()
        print(f'{partitions[i]} component completed', end-start)
        print('++++++++++++++++++++++++++++++++++++++++++++')
    print('representation afforded by the top degree homology is')
    print(rep)
    return rep
    