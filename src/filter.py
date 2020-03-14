#! /usr/bin/python2.7
import os
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import pdb

def read_file(file_name):
    try:
        seqs = []
        alignment = AlignIO.read(file_name, 'fasta') # Control: If the format of the content is wrong, the program would jump to 'except' here.
        for record in alignment:
            if len(record.seq) == 0:
                raise Exception  # Control: Check whether the sequence is empty.
            seqs_empty = SeqRecord(Seq(''), record.id, record.name, record.description)   
            seqs.append(seqs_empty)  # Create a list without sequence to store the unnoisy column.
    
    except:
        print('Error occurs, Please check the content of the input file.')
        sys.exit()  

    return seqs, alignment


def find_noise(column):
    # Count the amount of every amino acid and store them in the dictionary.
    aa_dict = {'A':0,'R':0,'N':0,'D':0,'C':0,'Q':0,'E':0,'G':0,'H':0,'I':0,'L':0,'K':0,'M':0,'F':0,'P':0,'S':0,'T':0,'W':0,'Y':0,'V':0, '-':0}
    for l in column:
        aa_dict[l] = aa_dict[l] + 1    

    # Start to find the noisy column.
    noisy = False
    aa_total = int(0)
    indels = int(0) 
    aa_unique = int(0)
    aa_morethan2 = int(0)
    for aa in aa_dict:
        if aa != '-': 
            aa_total = aa_total + aa_dict[aa] # Count the total amount of amino acids.
            if aa_dict[aa] == 1:
                aa_unique = aa_unique + 1 # Count the amount of unique amino acids.
            if aa_dict[aa] > 2:
                aa_morethan2 = aa_morethan2 + 1 # Count the amount of amino acids that appear more than twice.     
        else:
            indels = aa_dict['-'] # Count the amount of indels.    
    
    # Check the column according to the 3 conditions.
    if float(indels)/len(column) > 0.5: # Condition 1
        noisy = True
    if float(aa_unique)/aa_total >= 0.5: # Condition 2
        noisy = True
    if aa_morethan2 == 0: #Condition 3
        noisy = True

    return noisy

def noise_remove(new_path, result_path, file_name):
    seqs,alignment = read_file(new_path + '/' + file_name)
    ALL = True # Check whether all the columns are noisy and removed. 
    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i] # Get every column in the alignment.
        #If the function find_noise returns 'False', add the column to the empty list.
        if not find_noise(column): 
            ALL = False
            for j in range(len(seqs)):
                seqs[j] = seqs[j] + column[j]
    if ALL:
        print('Error occurs, all columns are noisy')    
        sys.exit()
    else:         
        outfile_name = result_path + '/unnoisy_' + str(file_name)
        output = MultipleSeqAlignment(seqs)
        AlignIO.write(output, outfile_name, 'fasta')

####### Main function #######

# The input is 'appbio11', the name of the directory of test data.
dir_name = sys.argv[1] 
path = '/home/xlbbbbbb/project/' + dir_name
s = []

# Store the names of sub-directories in a empty list.
for sub_name in os.walk(path):
    sub = sub_name[0][len(path)+1:]
    if sub != '':
        s.append(sub)

# Create the 'result' directory.
os.makedirs('/home/xlbbbbbb/project/result')
#pdb.set_trace()
for sub in s:
    new_path = path + '/' + sub
    result_path = '/home/xlbbbbbb/project/result/' + sub 
    os.makedirs(result_path) # Create sub-directories in the 'result' directory.
    for file_name in os.listdir(new_path):
        if file_name.endswith('.msl') and file_name[0] != '.':
            # Filter the alignment and write the output to a new file. 
            noise_remove(new_path, result_path, file_name)   
    















