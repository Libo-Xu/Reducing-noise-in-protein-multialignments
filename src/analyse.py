#! /usr/bin/python2.7
import sys
import os
import tempfile
import dendropy
from dendropy import Tree
from dendropy.calculate import treecompare
import pdb

tns = dendropy.TaxonNamespace()

def evaluate(ref,file_name):

    # To store the data during the process, we create two temporary files.
    tmp1 = tempfile.mkstemp()
    tmp2 = tempfile.mkstemp()

    # Use the commands of fastprot and fnj.
    # The output of the FastPhylo programs is in file 'tmp2'.
    os.system("fastprot -m -o " + tmp1[1] + " " + file_name)
    os.system("fnj -O newick -m FNJ -o " + tmp2[1] + " " + tmp1[1])
    
    #Use Dendropy to compare the trees.
    in_tree = Tree.get_from_stream(os.fdopen(tmp2[0]), schema='newick', taxon_namespace=tns)
    ref_tree = Tree.get_from_path(ref, schema='newick', taxon_namespace=tns)
    sym_diff = treecompare.symmetric_difference(ref_tree, in_tree)

    return sym_diff


####### Main function #######

# The input should be the name of sub-directorys in the 'result' directory, like 'asymmetric_0.5'.
sub_name = sys.argv[1]

path1 = '/home/xlbbbbbb/project/appbio11' + '/' + sub_name
path2 = '/home/xlbbbbbb/project/result' + '/' + sub_name

# Create 3 empty lists to store the symmetric difference between the original alignment tree and reference tree,
# the symmetric difference between noise-reduced tree and reference tree,
# and the difference between these 2 values.
original_sd = []
reduced_sd = []
diff = []

# Find the reference tree.
for filename in os.listdir(path1):
    if filename.endswith('.tree') and filename[0] != '.':
        ref_path = path1 + '/' + filename
#pdb.set_trace()
# Compare each files in the sub-directory and store the results in the corresponding list.
for filename in os.listdir(path1):
    if filename.endswith('.msl') and filename[0] != '.':
        detail_path1 = path1 + '/' + filename
        original = evaluate(ref_path,detail_path1)
        original_sd.append(original)
        for _filename in os.listdir(path2):
            if _filename == 'unnoisy_' + str(filename):
                detail_path2 = path2 + '/' + _filename
                reduced = evaluate(ref_path,detail_path2)  
                reduced_sd.append(reduced) 
                diff.append(original - reduced)




####### Begin to analyse the data #######

# The recovery frequency
f1 = int(0)
f2 = int(0)
for sd in original_sd:
    if sd == 0:
        f1 = f1 + 1
for sd in reduced_sd:
    
    if sd == 0:
        f2 = f2 + 1
print(sub_name + ' :  original:' + str(f1) + ' reduced:' + str(f2))

# The average of the difference between the original_sd and reduced_sd of each sub-directory.
# And the proportion of the case where the symmetric difference of the tree after noise removal is smaller than that of the original tree.
i = int(0)
_sum = int(0)
for d in diff:
    _sum = _sum + d
    if d > 0: 
        i = i + 1
avg = float(_sum)/300
print('Average of difference: %.2f'%avg)
print('Proportion: %.2f%%'%(float(i)/300*100))





        
