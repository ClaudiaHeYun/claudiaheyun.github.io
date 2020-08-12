import os
os.system('sage --preparse ThetaGraph.sage')
os.system('mv ThetaGraph.sage.py ThetaGraph.py')
import time
from ThetaGraph import *

#-------------Part I: calculating S_n-equivariant homology----------

#generating theta graphs; producing boundary maps

t1 = time.time()
dimTop = {}
counter = 0

for pi in symmetric(5):
    graph = generateThetaGraph(3,1,1,False,False,pi)
    if graph.isNormal():
        dimTop.update({repr(graph.markings):(counter,graph)})
        dimTop.update({counter:graph})
        counter = counter + 1

for pi in symmetric(5):
    graph = generateThetaGraph(2,2,1,False,False,pi)
    if graph.isNormal():
        dimTop.update({repr(graph.markings):(counter,graph)})
        dimTop.update({counter:graph})
        counter = counter + 1

dimTopM1 = {}
counter = 0

for pi in symmetric(5):
    graph = generateThetaGraph(2,1,1,False,True,pi)
    if graph.isNormal():
        dimTopM1.update({repr(graph.markings):(counter,graph)})
        dimTopM1.update({counter:graph})
        counter = counter + 1

dimTopM2 = {}
counter = 0
for pi in symmetric(5):
    graph = generateThetaGraph(1,1,1,True,True,pi)
    if graph.isNormal():
        dimTopM2.update({repr(graph.markings):(counter,graph)})
        dimTopM2.update({counter:graph})
        counter = counter + 1

t2 = time.time()
print 'generating lists of theta graphs'
print t2-t1

# generating boundary maps
dTop = matrix(QQ,int(len(dimTopM1)/2),int(len(dimTop)/2))
for j in range(int(len(dimTop)/2)):
    contractions = contract(dimTop[j])
    for graph in contractions:
        nGraph = normalize(graph)
        index = dimTopM1[repr(nGraph.markings)][0]
        entry = parity(nGraph.edges,graph.edges)
        dTop[index,j] = entry*graph.sign

dLower = matrix(QQ,int(len(dimTopM2)/2),int(len(dimTopM1)/2))
for j in range(int(len(dimTopM1)/2)):
    contractions = contract(dimTopM1[j])
    for graph in contractions:
        nGraph = normalize(graph)
        index = dimTopM2[repr(nGraph.markings)][0]
        entry = parity(nGraph.edges,graph.edges)
        dLower[index,j] = entry*graph.sign

t3 = time.time()
print 'generating boundary maps as matrices'
print t3-t2

# computing kernels
kerTop = dTop.right_kernel()
kerLower = dLower.right_kernel()
t4 = time.time()
print 'computing kernels of differentials'
print t4-t3

#calculating the character of the chain group and the top degree homology group

permList = [[0,1,2,3,4],[1,0,2,3,4],[1,0,3,2,4],[1,2,0,3,4],[1,2,0,4,3],[1,2,3,0,4],[1,2,3,4,0]]

def trace(mat):
    tr = 0
    for i in range(mat.nrows()):
        tr = tr+mat[i,i]
    return tr

t5 = time.time()

C_top = []    # chain group C_top as S_n rep
H_top = []    # top degree homology as S_n rep
dim = int(len(dimTop)/2)
for permutation in permList:
    perm = Permutation(permutation)
    perm_matrix = matrix(QQ,dim,dim)
    for i in range(dim):
        permuted = permute(perm,dimTop[i])
        normalized = normalize(permuted)
        row_index = dimTop[repr(normalized.markings)][0]
        entry = parity(normalized.edges,permuted.edges)
        
        perm_matrix[row_index,i] = entry
    C_top.append(trace(perm_matrix))
    
    # This calculates how each permutation permute basis elements
    tr=0
    for i in range(kerTop.dimension()):
        permuted_ker = perm_matrix*kerTop.basis()[i]
        linear_combo = kerTop.coordinate_vector(permuted_ker)
        tr=tr+linear_combo[i]
    H_top.append(tr)
    print(str(perm)+' done')
t6 = time.time()

print 'character of top degree chain group',C_top
print 'character of top degree homology',H_top
print t6-t5

# do the same operation for kerLower
kerLower_char = []
dimLower = int(len(dimTopM1)/2)
for permutation in permList:
    perm = Permutation(permutation)
    perm_matrix = matrix(QQ,dimLower,dimLower,sparse=True)
    for i in range(dimLower):
        permuted = permute(perm,dimTopM1[i])
        normalized = normalize(permuted)
        row_index = dimTopM1[repr(normalized.markings)][0]
        entry = parity(normalized.edges,permuted.edges)
        
        perm_matrix[row_index,i] = entry
        
    tr=0
    for i in range(kerLower.dimension()):
        permuted_ker = perm_matrix*kerLower.basis()[i]
        linear_combo = kerLower.coordinate_vector(permuted_ker)
        tr=tr+linear_combo[i]
    kerLower_char.append(tr)
    print(str(perm)+' done')
t7 = time.time()

print 'character of the kernel of lower degree differential',kerLower
print t7-t6

#------------------Part II; Calculating irreducible subrepresentations-----------------

#characters of irreducible representations

# character function of V_41
def chi_41(pi):
    c11111 = Permutation([0,1,2,3,4]).cycle_structure
    c2111 = Permutation([1,0,2,3,4]).cycle_structure
    c221 = Permutation([1,0,3,2,4]).cycle_structure
    c311 = Permutation([1,2,0,3,4]).cycle_structure
    c32 = Permutation([1,2,0,4,3]).cycle_structure
    c41 = Permutation([1,2,3,0,4]).cycle_structure
    c5 = Permutation([1,2,3,4,0]).cycle_structure
    
    if pi.cycle_structure == c11111:
        return 4
    elif pi.cycle_structure == c2111:
        return 2
    elif pi.cycle_structure == c221:
        return 0
    elif pi.cycle_structure == c311:
        return 1
    elif pi.cycle_structure == c32:
        return -1
    elif pi.cycle_structure == c41:
        return 0
    elif pi.cycle_structure == c5:
        return -1

def chi_32(pi):
    c11111 = Permutation([0,1,2,3,4]).cycle_structure
    c2111 = Permutation([1,0,2,3,4]).cycle_structure
    c221 = Permutation([1,0,3,2,4]).cycle_structure
    c311 = Permutation([1,2,0,3,4]).cycle_structure
    c32 = Permutation([1,2,0,4,3]).cycle_structure
    c41 = Permutation([1,2,3,0,4]).cycle_structure
    c5 = Permutation([1,2,3,4,0]).cycle_structure
    
    if pi.cycle_structure == c11111:
        return 5
    elif pi.cycle_structure == c2111:
        return 1
    elif pi.cycle_structure == c221:
        return 1
    elif pi.cycle_structure == c311:
        return -1
    elif pi.cycle_structure == c32:
        return 1
    elif pi.cycle_structure == c41:
        return -1
    elif pi.cycle_structure == c5:
        return 0

def chi_311(pi):
    c11111 = Permutation([0,1,2,3,4]).cycle_structure
    c2111 = Permutation([1,0,2,3,4]).cycle_structure
    c221 = Permutation([1,0,3,2,4]).cycle_structure
    c311 = Permutation([1,2,0,3,4]).cycle_structure
    c32 = Permutation([1,2,0,4,3]).cycle_structure
    c41 = Permutation([1,2,3,0,4]).cycle_structure
    c5 = Permutation([1,2,3,4,0]).cycle_structure
    
    if pi.cycle_structure == c11111:
        return 6
    elif pi.cycle_structure == c2111:
        return 0
    elif pi.cycle_structure == c221:
        return -2
    elif pi.cycle_structure == c311:
        return 0
    elif pi.cycle_structure == c32:
        return 0
    elif pi.cycle_structure == c41:
        return 0
    elif pi.cycle_structure == c5:
        return 1

# producing matrices for all permutations in S_5 in terms of theta graphs

start=time.time()
all_perm_mat = {}
for perm in symmetric(5):
    perm_matrix = matrix(QQ,dim,dim)
    for i in range(dim):
        permuted = permute(perm,dimTop[i])
        normalized = normalize(permuted)
        row_index = dimTop[repr(normalized.markings)][0]
        entry = parity(normalized.edges,permuted.edges)
        
        perm_matrix[row_index,i] = entry
    all_perm_mat[str(perm)] = perm_matrix
fin=time.time()
print 'producing matrices for element in S_5 wrt theta graphs took', fin-start

#calculating the projection matrix for each irreducible component using projection formula from Serre

dim = int(len(dimTop)/2)
projection_41 = matrix(QQ,kerTop.dimension(),kerTop.dimension())
projection_32 = matrix(QQ,kerTop.dimension(),kerTop.dimension())
projection_311 = matrix(QQ,kerTop.dimension(),kerTop.dimension())
start = time.time()
for pi in symmetric(5):
    # how does pi act on the kernel of the top differential?
    perm_matrix = all_perm_mat[str(pi)]
    action = []
    for i in range(kerTop.dimension()):
        permuted_ker = perm_matrix*kerTop.basis()[i]
        linear_combo = kerTop.coordinate_vector(permuted_ker)
        action.append(linear_combo)
    mat = matrix(QQ,action)
    projection_41 += chi_41(pi)*mat
    projection_32 += chi_32(pi)*mat
    projection_311 += chi_311(pi)*mat
end = time.time()
print('calculating projection matrix used', end-start)

#producing the irreducible representations

subrep_41 = projection_41*kerTop
subrep_32 = projection_32*kerTop
subrep_311 = projection_311*kerTop

#---------------Part III: constructing isomorphism from V_311 to S^311-------------------

#producing a more intuitive basis; call this the standard basis of V_311

# this permutation corresponds to (01234). cf. Section 4.3
perm = Permutation([0,3,4,2,1])
e0 = subrep_311.basis()[0]
new_e0 = all_perm_mat[str(perm)]*e0

std_basis=[new_e0]

# these permutations form S_{0,1,2}
std_perm=[Permutation([1,0,2,3,4]),Permutation([2,1,0,3,4]),Permutation([0,2,1,3,4]),Permutation([1,2,0,3,4]),Permutation([2,0,1,3,4])]

for perm in std_perm:
    permuted_e = all_perm_mat[str(perm)]*new_e0
    std_basis.append(permuted_e)

# perform base change

base_change = matrix(QQ,6,6)
for i in range(6):
    linear_combo = subrep_311.coordinate_vector(std_basis[i])
    base_change[i]=linear_combo
base_change = base_change.transpose().inverse()

# produce a matrix in terms of the std_basis for each element of the group
v311 = {}
v311_inverse = {}
start = time.time()
for perm in symmetric(5):
    perm_mat = matrix(QQ,6,6)
    for i in range(6):
        e = std_basis[i]
        permuted_e = all_perm_mat[str(perm)]*e
        linear_combo = subrep_311.coordinate_vector(permuted_e)
        c_basis_vec = matrix(QQ,linear_combo).transpose()
        std_basis_vec = base_change*c_basis_vec
        perm_mat[i] = std_basis_vec.transpose()
    perm_mat = perm_mat.transpose()
    v311[str(perm)] = perm_mat
    v311_inverse[str(perm)] = perm_mat.inverse()
fin = time.time()
print 'producing matrices wrt our standard basis took', fin-start

# now we need to do the same but using basis of the Specht module in terms of polytabloids.
# manually calculate the matrices of transpositions
sp311_inverse = {str(Permutation([0,1,2,3,4])):matrix.identity(QQ,6)}
sp311 = {str(Permutation([0,1,2,3,4])):matrix.identity(QQ,6)}
m01 = matrix(QQ,[[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[-1,-1,0,-1,0,0],[1,0,-1,0,-1,0],[0,1,1,0,0,-1]])
m12 = matrix(QQ,[[1,0,0,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,0,0,-1]])
m23 = matrix(QQ,[[0,1,0,0,0,0],[1,0,0,0,0,0],[0,0,-1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,0,1],[0,0,0,0,1,0]])
m34 = matrix(QQ,[[-1,0,0,0,0,0],[0,0,1,0,0,0],[0,1,0,0,0,0],[0,0,0,0,1,0],[0,0,0,1,0,0],[0,0,0,0,0,1]])
m02 = matrix(QQ,[[1,0,0,0,0,0],[-1,-1,0,-1,0,0],[1,0,-1,0,-1,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,-1,-1,-1]])
m03 = matrix(QQ,[[-1,-1,0,-1,0,0],[0,1,0,0,0,0],[0,-1,-1,0,0,1],[0,0,0,1,0,0],[0,0,0,-1,-1,-1],[0,0,0,0,0,1]])
m04 = matrix(QQ,[[-1,0,1,0,1,0],[0,-1,-1,0,0,1],[0,0,1,0,0,0],[0,0,0,-1,-1,-1],[0,0,0,0,1,0],[0,0,0,0,0,1]])
m13 = matrix(QQ,[[0,0,0,1,0,0],[0,1,0,0,0,0],[0,0,0,0,0,-1],[1,0,0,0,0,0],[0,0,0,0,-1,0],[0,0,-1,0,0,0]])
m14 = matrix(QQ,[[0,0,0,0,-1,0],[0,0,0,0,0,-1],[0,0,1,0,0,0],[0,0,0,-1,0,0],[-1,0,0,0,0,0],[0,-1,0,0,0,0]])
m24 = matrix(QQ,[[0,0,-1,0,0,0],[0,-1,0,0,0,0],[-1,0,0,0,0,0],[0,0,0,0,0,1],[0,0,0,0,1,0],[0,0,0,1,0,0]])
transposition_mat = {(0,1):m01, (0,2):m02, (0,3):m03, (0,4):m04, (1,2):m12, (1,3):m13, (1,4):m14, (2,3):m23, (2,4):m24, (3,4):m34}

start = time.time()
for perm in symmetric(5):
    if perm != Permutation([0,1,2,3,4]):
        decomp = perm.transpositions()
        mat = transposition_mat[decomp[-1]]
        for i in range(len(decomp)-2,-1,-1):
            tup = decomp[i]
            mat = transposition_mat[tup]*mat
        sp311[str(perm)] = mat
        sp311_inverse[str(perm)] = mat.inverse()
fin = time.time()
print 'producing matrices wrt to polytabloids took', fin-start

# now use the formula in Serre's book
h = matrix(QQ,[[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
h_0 = matrix(QQ,6,6)
for perm in symmetric(5):
    h_0 += sp311_inverse[str(perm)]*h*v311[str(perm)]

print 'Isomorphism from V_311 to S^311 is given by'
print h_0