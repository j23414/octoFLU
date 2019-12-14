# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 21:25:20 2019

@author: Mazeller
"""

import sys
import getopt
import dendropy
#from pathlib import Path
from collections import OrderedDict
from operator import itemgetter   

# Print help text, including parameters
def helpText():
    print("""Usage: treedist.py -i inputtree.tre

Depends on dendrop.py, assumes newick format.
    """)
        
    
# treedist function
def treedist(treeFile = False, columnAnnotated = [5], outfile=None):
    """Determines clade based on nearest neighbor classifier
    
    Given a tree where references have been annotated with phylogenetic clade
    information and delimited by pipes, print out the non-piped strains and their
    nearest neighbor reference strain clade.
    
      :param treeFile The name of a nexus file containing query and reference strains
      :param columnAnnotated List of positions in the pipe-delimited reference strain names
      :param outfile The name of an output file. This is optional, if None then will print to stdout
    
    """
    if outfile is not None:
        ofile = open(outfile, 'a+')
    else:
        ofile=None
    
    #Error handeling for required args
    if not treeFile:
        print("Tree file is required.\n")
        helpText()
        sys.exit(2)
    
    # attempt at dealing with windows vs linux style paths
    #treeFile = Path(treeFile)

    #Attempt to load tree file
    tree = dendropy.Tree.get(path = treeFile, schema="newick")
    
    #Grab distances
    pdm = tree.phylogenetic_distance_matrix()
    pdma = pdm.as_data_table()

    #Iterate through col names. If missing the "tag", find the nearest neighbor that has a tag.
    for idx1, taxon1 in enumerate(tree.taxon_namespace):
        
        #Check if clade designation is present, skip if so
        defSplit = str(taxon1).split("|")
        if(len(defSplit) >= max(columnAnnotated)):    #Assumption! One delimiter (JC: can we make this more robust?) 
            continue
        #print(str(taxon1))
        
        #Find the nearest neighbor with a clade label otherwise. Starting at index 1 incase 100% identity label match
        dist = pdma._data[str(taxon1)[1:-1]]
        orderedDist = OrderedDict(sorted(dist.items(), key = itemgetter(1)))
        for distance in orderedDist:
            compSplit = distance.split("|")
            if(len(compSplit) >= max(columnAnnotated)):    
                outString = str(taxon1)[1:-1].replace(" ","_") + "\t"
                for i in columnAnnotated:
                    outString += compSplit[i - 1] + "\t"     #Index shift, minus 1 (probably could use map instead of for)
                print(outString, file=ofile)
                break
            
    if outfile is not None:
        ofile.close()


#Grab R from commandline arguments, or print usage
#Main execution
def main():
    argv = sys.argv[1:]

    try:
        opts, args = getopt.getopt(argv,"i:c:h",["input=","column=","help"])
    except getopt.GetoptError:
        helpText()
        sys.exit(2)
    
    #Assign local variables    
    treeFile = False
    columnAnnotated = [5]
    
    #Process command line args    
    for opt, arg in opts:
        if opt in ('-h',"help"):
            helpText()
            sys.exit()
        elif opt in ("-i", "-input"):
            treeFile = arg
        elif opt in ("-c", "-column"):
            columnAnnotated = list(map(int, arg.split(",")))
            
    #Call treedist function
    treedist(treeFile=treeFile, columnAnnotated=columnAnnotated)

    
            
if __name__ == "__main__": main()
