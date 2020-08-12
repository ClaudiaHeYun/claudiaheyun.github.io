import copy
from sympy.combinatorics import Permutation
from sympy.combinatorics.generators import symmetric

class ThetaGraph:
    def __init__(self, markings, edges):
        # markings is a list of lists: [[first branch],[second branch],[third branch],[left,right]]
        # edges is a list of lists: [[first branch],[second branch],[third branch]]
        self.markings = markings
        self.edges = edges
        self.sign = 1
    
    @classmethod
    def fromMarkings(cls, markings):
        # edges are ordered in increasing order
        first_br = range(0, len(markings[0])+1)
        second_br = range(len(markings[0])+1, len(markings[0])+len(markings[1])+2)
        third_br = range(len(markings[0])+len(markings[1])+2, len(markings[0])+len(markings[1])+
                              len(markings[2])+3)
        edges = [first_br, second_br, third_br]
        return cls(markings, edges)

    def see(self):
        print(str(self.sign)+"\n"+"markings: "+str(self.markings)+"\n"+"edges: "+str(self.edges))
    
    
    def changeSign(self, new):
        self.sign = int(new)

    # put in standard order:
    # number of markings in the first branch >= number of markings in the second >= number of markings in the third
    # when one of left and right exists, put it on the RHS, because -1 < every possible marking
    # when both left and right exist, put the smaller # on the LHS
    # when neither of them exists, first number of the first branch < last number of the first branch
    # first # of the first branch < first # of the second branch < first # of the third branch

    # this seems pretty complicated, not sure if it's better than the current solution
    
    # boolean that determines if a graph is in normal form
    def isNormal(self):
        boo = True
        
        # if there isn't marking on the left or right, we already manually put
        # the branches in a way such that the top has the most markings, then
        # the middle, and the bottom one has the fewest. So we only need to
        # worry reflections and switching branches when they have the same 
        # number of markings.
        
        # ha, we still need to consider the markings on branches, because
        # although I specified the numbers when I produced theta graphs, it is
        # not nec. true that after contractions they still have the desired
        # number on each branch
        
        # if it's not true that #markings on top >= #middle >= #bottom, false
        if len(self.markings[0]) < len(self.markings[1]) or len(self.markings[1]) < len(self.markings[2]):
            boo = False
            
        # #top >= #middle >= #bottom; now we consider reflection & switching
        
        # there are no markings on the side vertices
        # in this case, we really want the upper left marking to be the smallest, and that determines the reflection
        # and then we decide switching branches
        elif self.markings[3] == [-1,-1]:
            
            # if the first marking from the left on the first branch is larger than the last marking, false
            if self.markings[0][0] > self.markings[0][-1]:
                    boo = False
            
            elif len(self.markings[0]) == len(self.markings[1]) == len(self.markings[2]):
                if (self.markings[0][0] > self.markings[1][0] or self.markings[0][0] > self.markings[1][-1]
                   or self.markings[0][0] > self.markings[2][0] or self.markings[0][0] > self.markings[2][-1]):
                    boo = False
                elif self.markings[1][0] > self.markings[2][0]:
                    boo = False
            
            elif len(self.markings[0])== len(self.markings[1]):
                if self.markings[0][0] > self.markings[1][0] or self.markings[0][0] > self.markings[1][-1]:
                    boo = False
            
            elif len(self.markings[1]) == len(self.markings[2]):
                if self.markings[1][0] > self.markings[2][0]:
                    boo = False
        
        # now there is at least one side marking
        # always want the left marking to be smaller than right marking
        # this includes the case where there is nothing on the left, because then it has marking -1
        elif self.markings[3][0] > self.markings[3][1]:
            boo = False
        
        # if there is only one side marking
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
            if len(self.markings[0]) == len(self.markings[1]) == len(self.markings[2]):
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
    #newGraph = graph.deepcopy()
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
    #graph12 = graph.deepcopy()
    graph12 = deepcopy(graph)
    graph12.markings[0], graph12.markings[1] = graph12.markings[1], graph12.markings[0]
    graph12.edges[0], graph12.edges[1] = graph12.edges[1], graph12.edges[0]
    configurations.append(graph12)
    
    # (13)
    #graph13 = graph.deepcopy()
    graph13 = deepcopy(graph)
    graph13.markings[0], graph13.markings[2] = graph13.markings[2], graph13.markings[0]
    graph13.edges[0], graph13.edges[2] = graph13.edges[2], graph13.edges[0]
    configurations.append(graph13)
    
    # (23)
    #graph23 = graph.deepcopy()
    graph23 = deepcopy(graph)
    graph23.markings[1], graph23.markings[2] = graph23.markings[2], graph23.markings[1]
    graph23.edges[1], graph23.edges[2] = graph23.edges[2], graph23.edges[1]
    configurations.append(graph23)
    
    # (123)
    #graph123 = graph.deepcopy()
    graph123 = deepcopy(graph)
    graph123.markings[0], graph123.markings[1], graph123.markings[2] = graph123.markings[2], graph123.markings[0], graph123.markings[1]
    graph123.edges[0], graph123.edges[1], graph123.edges[2] = graph123.edges[2], graph123.edges[0],  graph123.edges[1]
    configurations.append(graph123)
    
    # (132)
    #graph132 = graph.deepcopy()
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

# a and b are two edge sets (permutations of numbers, lists of lists). This determines the sign of the permutation that takes one to the other.
# example: a = [0,1,2,3,4,5,6,7,8], b = [1,0,2,3,4,5,6,7,8], then (0,1) takes one to the other, so the sign is -1.
def parity(a,b):
    pi = Permutation(flatten(a))
    sigma = Permutation(flatten(b))
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
    if normalA == normalB:
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
    #newGraph = graph.deepcopy()
    newGraph = deepcopy(graph)
    for i in range(len(newGraph.markings)):
        for j in range(len(newGraph.markings[i])):
            if newGraph.markings[i][j] != -1:
                newGraph.markings[i][j] = pi(newGraph.markings[i][j])
    
    return newGraph

# Given a permutation in S_n, this produces a theta graph with markings according to the given structure of the theta graph and the permutation.
# a,b,c are positive integers. a markings on the top branch, b on the middle, c on the bottom
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
    
    graph = ThetaGraph.fromMarkings(markings)
    
    return graph