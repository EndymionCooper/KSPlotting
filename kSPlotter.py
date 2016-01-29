#!/path/to/python

##########################################################################################
#                                                                                        #
#   Created by Endymion D. Cooper (endymion.dante.cooper@gmail.com)                      #
#   This is version 1.0 created on 15th Sept 2015                                        #
#                                                                                        #
#   kSPlotter.py takes a table of blast hits and plots distributions of pairwise         #
#   kS values for duplicated genes within a genome.                                      #
#                                                                                        #
#   Dependancies: written for python 2.7, builtins and standard libraries. Requires      #
#   biopython. The following must also be in your path:                                  #
#      - muscle                                                                          #
#      - codeml (from paml)                                                              #
#                                                                                        #
#   Citation information (and versions I use) for dependencies.                          #
#   MUSCLE v3.8.31                                                                       #
#   Edgar, R.C. Nucleic Acids Res 32(5), 1792-97.                                        #
#   PAML v4.8                                                                            #
#   Yang, Z. (2007) Molecular Biology and Evolution 24: 1586-1591                        #
#                                                                                        #
#   Make sure the path to your python installation is correctly specified above.         #
#                                                                                        #
#   To get usage information: kSPlotter.py --help                                        #
#                                                                                        #
##########################################################################################

import sys
import os
from Bio import SeqIO, Seq
from Bio.Alphabet import generic_dna
import subprocess
import argparse

##########################################################################################
#   Parse command line arguments.
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, 
   description='''

    ------------------------------------------------------------------------------------- 
   |                                  calculate_KS.py                                    |
    ------------------------------------------------------------------------------------- 

   Created by Endymion D. Cooper (endymion.dante.cooper@gmail.com)                      
   This is version 1.0 created on 16th Sept 2015                                        
                                                                                        
   kSPlotter.py takes a table of blast hits and generates kS values for inferred gene 
   duplications
                                                                                        
   Dependancies: written for python 2.7, builtins and standard libraries. Requires      
   biopython to process fasta files.
   
   First gene families are inferred by naive clustering on blast hits.
   Blast output should be formatted using < -outfmt "7 std qlen slen" >
   
   Next all pairwise kS values are computed within each gene family using:
       1. MUSCLE to generate a sequence alignment for each pair using amino acid sequences
       2. No longer required: [PAL2NAL to generate nucleotide alignments using the amino acid alignments]
       3. PAML (codeml) to calculate maximum likelihood estimates of the kS scores.

   Redundancy in kS scores within a gene family is reduced to approximate gene 
   duplications using a UPGMA like clustering method (following Maere et al. PNAS 2005
   DOI:10.1073/pnas.0501102102).                                                        
                                                                                        
   Maere et al. (2005) Supporting methods:                                              
                                                                                        
       " Correction for redundant KS values.                                            
         A gene family of n members originates from n-1 retained single gene            
         duplications, whereas the number of possible pairwise comparisons (KS          
         measurements) within a gene family is n(n-1)/2. To correct for the             
         redundancy of KS values when building the age distribution for duplicated      
         genes, we constructed tentative phylogenetic trees for each gene family        
         with an average linkage clustering algorithm using KS as a distance            
         measure, similar to the approach adopted by Blanc and Wolfe (1). Starting      
         from each gene as a separate cluster, the two clusters with the lowest         
         mean inter-cluster KS value (i.e. the mean of all observed KS values           
         (edges) between two clusters) were iteratively merged. The splits in           
         the resulting average linkage tree represent the n-1 retained duplication      
         events. For each split, the m KS measurements between the two merged           
         gene clusters were added to the KS distribution with a weight 1/m. In          
         other words, all KS estimates for a particular duplication event were          
         added to the KS distribution, while the total weight of a single               
         duplication event sums up to one.  "                                           
                                                                                        
   Two methods of calculating mean KS between sub-clusters are available.               
   Specify M1 (mode 1) or M2 (mode 2).                                                  
   The differences are explained below with an example.                                 
                                                                                        
   Example data:                                                                        
   <Cluster ID>  <Gene1 ID>  <Gene2 ID>  <KS>                                           
     cluster1        G1          G2      0.003                                          
     cluster1        G1          G3      0.326                                          
     cluster1        G1          G4      0.563                                          
     cluster1        G2          G3      0.245                                          
     cluster1        G2          G4      0.637                                          
     cluster1        G3          G4      0.476                                          
                                                                                        
   Denote:  - the KS score between two sequences/subclusters as e.g. G1:G2.             
            - a subcluster of sequences as e.g. G1/2                                    
                                                                                        
   First iteration:                                                                     
   The smallest KS value is between G1 & G2 therefore join as G1/2. Tree:               
        G1 __                                                                           
        G2 __|  <-- Duplication event 1. KS score G1:G2 = 0.003                         
        G3                                                                              
        G4                                                                              
   And recalculate the KS scores by averaging (note that mode 1 = mode 2 here):         
    -------------------------------------------------------------------------           
   | New Pairs |            Mode 1            |            Mode 2            |          
   |-----------|------------------------------|------------------------------|          
   | G1/2  G3  | ((G1:G3)+(G2:G3))/2 = 0.2855 | ((G1:G3)+(G2:G3))/2 = 0.2855 |          
   | G1/2  G4  | ((G1:G4)+(G2:G4))/2 = 0.6    | ((G1:G4)+(G2:G4))/2 = 0.6    |          
   | G3    G4  | (raw value)         = 0.476  | (raw value)         = 0.476  |          
    -------------------------------------------------------------------------           
   Second iteration:                                                                    
   The smallest KS value is between G1/2 and G3 therefore join as G1/2/3. Tree:         
        G1 __                                                                           
        G2 __|__                                                                        
        G3 _____|  <-- Duplication event 2. KS score G1/2:G3 = 0.2855                   
        G4                                                                              
   And recalculate the KS scores by averaging (mode 1 =/= mode 2):                      
    ----------------------------------------------------------------------------------  
   |  New Pairs |                Mode 1               |            Mode 2             | 
   |------------|-------------------------------------|-------------------------------| 
   | G1/2/3  G4 | ((G1:G4)+(G2:G4)+(G3:G4))/3 = 0.559 | ((G1/2:G4)+(G3:G4))/2 = 0.538 | 
    ----------------------------------------------------------------------------------  
                                                                                        
   Third iteration 3:                                                                   
   The smallest KS value is between G1/2/3 and G4 therefore join as G1/2/3/4. Tree:     
        G1 __                                                                           
        G2 __|__                                                                        
        G3 _____|__                                                                     
        G4 ________|   <-- Duplication event 3. KS score G1/2/3:G4 = 0.559 (mode 1)     
                                                                   = 0.538 (mode 2)     
    ------------------------------------------------------------------------------------- 
''')
parser.add_argument('-b', nargs=1, required=True,
               help='specify the input blast table',metavar='<filename>')
parser.add_argument('-aa', nargs=1, default=['translate'],
               help='specify the input amino acid sequence file',
               metavar='<filename>')
parser.add_argument('-nt', nargs=1, required=True,
               help='specify the input nucleotide sequence file',
               metavar='<filename>')
parser.add_argument('-o', nargs=1, required=True,
               help='output file prefix',metavar='<filename>')
parser.add_argument('-I','--identity', nargs=1, default=['50'],
               help='filter blast table using sequence identity [default = 50]',
               metavar='<int>')
parser.add_argument('-E','--e-value', nargs=1, default=['1e-10'],
               help='filter blast table using the e-value [default = 1e-10]',
               metavar='<e-value>')
parser.add_argument('-L','--length_ratio', nargs=1, default=['0.8'],
               help='''filter blast table using ratio of query and subject lengths, 
               where the ratio is shortest divided by longest [default = 0.8]''',
               metavar='<float>')
parser.add_argument('-A','--alignment', nargs=1,default=['0.8'],
               help='''filter blast table using a minimum alignment length, where the
               alignment length is expressed as a fraction of the shortest of the 
               query and subject sequences [default = 0.8]''',
               metavar='<int>')
parser.add_argument('-R','--reduce', nargs=1,default=['M1'],
               help='mode for reducing reduntant kS values from gene families',
               metavar='<M1/M2>')
parser.add_argument('-ER','--error_tolerance', nargs=1,default=['10'],
               help='number of errors to allow before terminating',
               metavar='<int>')
args = parser.parse_args()
input_blast_table = args.b[0]
input_amino_acids = args.aa[0]
input_CDS = args.nt[0]
output_prefix = args.o[0]
identity_filter = float(args.identity[0])
e_value_filter = float(args.e_value[0])
length_ratio_filter = float(args.length_ratio[0])
alignment_length_filter = float(args.alignment[0])
redundant_kS_mode = str(args.reduce[0])
error_tolerance = int(args.error_tolerance[0])

##########################################################################################
#   Read inputs:   1. Blast table.
#                  2. Amino acid sequence file
#                  3. Corresponding CDS

# 1. Blast table.
file = open(input_blast_table, "r")
file_lines = file.readlines()
file.close()

# 2. Amino acid sequence file
if input_amino_acids != 'translate':
   aa_file = open(input_amino_acids,"r")
   aa_dict = SeqIO.to_dict(SeqIO.parse(aa_file, "fasta"))
   aa_file.close()
else:
   print "\nNo amino acid file provided. Translating nucleotides..."
   nt_file = open(input_CDS,"r")
   aa_file = open(output_prefix+"translated_amino_acids.fasta","w")
   for record in SeqIO.parse(nt_file, "fasta"):
      print >>aa_file,">"+record.id
      print >>aa_file,Seq.Seq(str(record.seq), generic_dna).translate(table=1)
   nt_file.close()
   aa_file.close()
   aa_file = open(output_prefix+"translated_amino_acids.fasta","r")
   aa_dict = SeqIO.to_dict(SeqIO.parse(aa_file, "fasta"))
   aa_file.close()
   os.remove(output_prefix+"translated_amino_acids.fasta")


# 3. Corresponding CDS file
nt_file = open(input_CDS,"r")
nt_dict = SeqIO.to_dict(SeqIO.parse(nt_file, "fasta"))
nt_file.close()

# log file to hold stdout & stderr
log_file = open(output_prefix+"log.txt","w")

# error counter
ERRORS = 0

##########################################################################################
#   Extract pairs from blast table.
def parse_blast(lines):
   pairs = []
   for line in lines:
      # Ignore comment lines in blast table.
      if line.startswith("#"):
         continue
      else:
         ls = line.split()
         # Alignment length
         a_length = float(ls[3])
         # Query length
         q_length = float(ls[12])
         # Subject length
         s_length = float(ls[13])
         e_value = float(ls[10])
         identity = float(ls[2])
         if ((min(s_length,q_length)/max(s_length,q_length) >= length_ratio_filter) and 
            (a_length/min(s_length,q_length) >= alignment_length_filter) and
            (identity >= identity_filter) and
            (e_value <= e_value_filter)):
            if set([ls[0],ls[1]]) not in pairs:
               if len(list(set([ls[0],ls[1]]))) > 1:
                  pairs.append(set([ls[0],ls[1]]))
               else:
                  continue
         else:
            continue
   return pairs

##########################################################################################
#   Build gene family clusters from blast pairs.
def build_clusters(temp):
   lt1 = len(temp)
   for x in temp:
      temp_list = []
      for y in temp:
         if x != y:
            if len(x & y) > 0:
               for item in list(x | y):
                  temp_list.append(item)
               temp.remove(y)
      if len(temp_list) > 0:
         temp.append(set(temp_list))
         temp.remove(x)
   if lt1 > len(temp):
      return build_clusters(temp)
   else:
      return temp

##########################################################################################
#   Write a codeml.ctl file.
def write_codeml_ctl(input, output):
   CODEML = open("codeml.ctl","w")
   print >> CODEML, "seqfile = ",input
   print >> CODEML, "treefile = stewart.trees"
   print >> CODEML, "outfile = ",output
   print >> CODEML, "noisy = 0"
   print >> CODEML, "verbose = 0"
   print >> CODEML, "runmode = -2"
   print >> CODEML, "seqtype = 1"
   print >> CODEML, "CodonFreq = 2"
   print >> CODEML, "clock = 0"
   print >> CODEML, "aaDist = 0"
   print >> CODEML, "model = 0"
   print >> CODEML, "NSsites = 0  "
   print >> CODEML, "icode = 0"
   print >> CODEML, "Mgene = 0"
   print >> CODEML, "fix_kappa = 0"
   print >> CODEML, "kappa = 2"
   print >> CODEML, "fix_omega = 0"
   print >> CODEML, "omega = .4"
   print >> CODEML, "fix_alpha = 1"
   print >> CODEML, "alpha = 0"
   print >> CODEML, "Malpha = 0"
   print >> CODEML, "ncatG = 8"
   print >> CODEML, "getSE = 0"
   print >> CODEML, "RateAncestor = 1"
   print >> CODEML, "Small_Diff = .5e-6"
   print >> CODEML, "cleandata = 1"
   print >> CODEML, "method = 0"
   CODEML.close()

##########################################################################################
#   Extract a pair of amino acid sequence records from input, and
#   extract corresponding pair of nucleotide sequence records from input.
def get_aa_seq_pair(seq_1,seq_2):
   sequence_pair = []
   sequence_pair.append(aa_dict[seq_1])
   sequence_pair.append(aa_dict[seq_2])
   return sequence_pair
def get_nt_seq_pair(seq_1,seq_2):
   sequence_pair = []
   sequence_pair.append(nt_dict[seq_1])
   sequence_pair.append(nt_dict[seq_2])
   return sequence_pair


##########################################################################################
#   get a nucleotide alignment from an amino acid alignment
def aa2nt_aln(aa_aln,nt_fasta,outfile):
	# store the nucleotide sequences in a dictionary
	nt_file = open(nt_fasta,"r")
	nt_dict = SeqIO.to_dict(SeqIO.parse(nt_file, "fasta"))
	nt_file.close()
	# store the nucleotide sequences in a dictionary
	aa_file = open(aa_aln,"r")
	aa_dict = SeqIO.to_dict(SeqIO.parse(aa_file, "fasta"))
	aa_file.close()

	# read through an aa sequence one site at a time
	# if the site is not a gap insert the corresponding codon into a new nt sequence
	# if it is a gap, insert three gap characters
	seq_name=1
	for seq in aa_dict:
		new_seq=""
		counter=0
		for character in aa_dict[seq]:
			if character != '-':
				if character != '*':
					new_seq = new_seq+nt_dict[seq].seq[counter:counter+3]
					counter = counter+3
				else:
					new_seq = new_seq+"---"
					counter = counter+3					
			else:
				new_seq = new_seq+"---"
		print >>outfile, ">seq"+str(seq_name)
		print >>outfile, new_seq
		seq_name=seq_name+1
	outfile.close()

##########################################################################################
#   Calculate the kS scores
def get_KS(clusters):
   global ERRORS
   counter = 0
   output_list = []
   for cluster in clusters:
      if ERRORS > error_tolerance:
         print >> sys.stderr, "\n\nEXECUTITION HALTED. TOO MANY ERRORS.\n\n"
         sys.exit()
      counter = counter + 1
      cluster_ID = "CL"+str(100000+counter)
      print "Processed cluster ",counter," of ", len(clusters), "clusters"
      temp_dict ={}
      for seq_1 in list(cluster):
         for seq_2 in list(cluster):
            if seq_1 != seq_2 and frozenset([seq_1,seq_2]) not in temp_dict:
               try:
                  sequence_pair_aa = get_aa_seq_pair(seq_1,seq_2)
                  output_aa = open(str(seq_1)+"_"+str(seq_2)+"_aa.fasta", "w")
                  SeqIO.write(sequence_pair_aa, output_aa, "fasta")
                  output_aa.close()
                  sequence_pair_nt = get_nt_seq_pair(seq_1,seq_2)
                  output_nt = open(str(seq_1)+"_"+str(seq_2)+"_nt.fasta", "w")
                  SeqIO.write(sequence_pair_nt, output_nt, "fasta")
                  output_nt.close()
                  try:
                     # Generate an amino acid sequence alignment using muscle.
                     subprocess.call(["muscle",
                                 "-in", str(seq_1)+"_"+str(seq_2)+"_aa.fasta",
                                 "-out", str(seq_1)+"_"+str(seq_2)+"_aa.aln"],
                                 stdout=log_file,stderr=log_file)
                     # Generate corresponding nucleotide alignment using PAL2NAL
                     #pal2nal_out = open(str(seq_1)+"_"+str(seq_2)+"_nt.aln", "w")
                     aa2nt_aln_out = open(str(seq_1)+"_"+str(seq_2)+"_nt.aln", "w")
                     aa2nt_aln(str(seq_1)+"_"+str(seq_2)+"_aa.aln",str(seq_1)+"_"+str(seq_2)+"_nt.fasta",aa2nt_aln_out)
                     #subprocess.call(["pal2nal.pl",
                     #            str(seq_1)+"_"+str(seq_2)+"_aa.aln", 
                     #            str(seq_1)+"_"+str(seq_2)+"_nt.fasta",
                     #            "-output", 
                     #            "paml", str(seq_1)+"_"+str(seq_2)+"_nt.aln"],
                     #            stdout=pal2nal_out,stderr=log_file)
                     #pal2nal_out.close()
                     if os.stat(str(seq_1)+"_"+str(seq_2)+"_nt.aln").st_size != 0:
                        try:
                           # Calculate kS values using codeml
                           write_codeml_ctl(str(seq_1)+"_"+str(seq_2)+"_nt.aln",
                                       "codeml.tmp")
                           subprocess.call(["codeml"],
                                       stdout=log_file,stderr=log_file)
                           KS_file = open("2ML.dS", "r")
                           KS_file_lines = KS_file.readlines()
                           KS_file.close()
                           KS = KS_file_lines[2].split()[1]
                           temp_dict[frozenset([seq_1,seq_2])]=KS
                        except:
                           print >> sys.stderr,"*****************************************"
                           print >> sys.stderr,"****** ERROR CALCULATING kS SCORES ******"
                           print >> sys.stderr," offending pair", seq_1, seq_2
                           print >> sys.stderr,'    -->  Using large value for kS [=999]'
                           print >> sys.stderr,"*****************************************"
                           ERRORS=ERRORS+1
                           temp_dict[frozenset([seq_1,seq_2])]=str(999)
                     else:
                        print >> sys.stderr,"*****************************************"
                        print >> sys.stderr,"********* ERROR RUNNING PAL2NAL *********"
                        print >> sys.stderr," offending pair", seq_1, seq_2
                        print >> sys.stderr,'    -->  Using large value for kS [=999]'
                        print >> sys.stderr,'Often this error is caused by conflict'
                        print >> sys.stderr,'between nucleotide and amino acid sequences'
                        print >> sys.stderr,"*****************************************"
                        ERRORS=ERRORS+1
                        temp_dict[frozenset([seq_1,seq_2])]=str(999)
                  except:
                     print >> sys.stderr,"*****************************************"
                     print >> sys.stderr,"********* ERROR RUNNING MUSCLE **********"
                     print >> sys.stderr," offending pair", seq_1, seq_2
                     print >> sys.stderr,'    -->  Using large value for kS [=999]'
                     print >> sys.stderr,"*****************************************"
                     ERRORS=ERRORS+1
                     temp_dict[frozenset([seq_1,seq_2])]=str(999)
               except:
                  print >> sys.stderr,"*****************************************"
                  print >> sys.stderr,"****** ERROR EXTRACTING SEQUENCES *******"
                  print >> sys.stderr," offending pair", seq_1, seq_2
                  print >> sys.stderr,'    -->  Using large value for kS [=999]'
                  print >> sys.stderr,'Check that your fasta file is correctly' 
                  print >> sys.stderr,'formatted and that the nucleotide and'
                  print >> sys.stderr,'amino acid headers match'
                  print >> sys.stderr,"*****************************************"
                  ERRORS=ERRORS+1
                  temp_dict[frozenset([seq_1,seq_2])]=str(999)
               # Clean up temporary files.
               try:
                  os.remove(str(seq_1)+"_"+str(seq_2)+"_aa.fasta")
               except:
                  continue
               try:
                  os.remove(str(seq_1)+"_"+str(seq_2)+"_nt.fasta")
               except:
                  continue
               try:
                  os.remove(str(seq_1)+"_"+str(seq_2)+"_aa.aln")
               except:
                  continue
               try:
                  os.remove(str(seq_1)+"_"+str(seq_2)+"_nt.aln")
               except:
                  continue
               try:
                  os.remove("2ML.dN")
                  os.remove("2ML.dS")
                  os.remove("2ML.t")
                  os.remove("2NG.dN")
                  os.remove("2NG.dS")
                  os.remove("2NG.t")
                  os.remove("codeml.ctl")
                  os.remove("codeml.tmp")
                  os.remove("rst")
                  os.remove("rst1")
                  os.remove("rub")
               except:
                  continue

      # Prepare per cluster output.
      for key in temp_dict:
         output_list.append([cluster_ID,list(key)[0],list(key)[1],temp_dict[key]])
   return output_list   


##########################################################################################
#   Reduce redundancy in the kS scores

def reduce_redundant_KS(output_file_lines):
   data_dict = {}
   for item in output_file_lines:
      if item[0] not in data_dict:
         data_dict[item[0]]=[item]
      else:
         data_dict[item[0]].append(item)
   return data_dict

def remove_redundant_KS(lines):
    # list to hold current clusters
    clusters = []
    # list to hold sequence names
    sequences = []
    # dictionary to hold cluster definitions
    cluster_defs = {}
    # dictionary to hold all KS scores
    all_KS = {}
    # dictionary to hold KS scores for current clusters
    current_KS = {}
    # list for final KS scores
    KS = []
    # a counter to append to subcluster IDs
    counter = 0

    # populate the lists and dictionaries with initial values.
    for line in lines:
        L2 = line
        L3 = frozenset([L2[1],L2[2]])
        # populate the KS value dictionary
        all_KS[L3] = L2[3]
        # populate the cluster and sequence lists
        if L2[1] not in clusters:
            clusters.append(L2[1])
            sequences.append(L2[1])
        if L2[2] not in clusters:
            clusters.append(L2[2])
            sequences.append(L2[2])
    # populate the cluster definitions
    for i in xrange(0,len(clusters)):
        cluster_defs[clusters[i]]=[clusters[i]]

    # define a function to get the KS scores for the current clusters
    # KS scores held in a dictionary and sent to get_min_KS to extract smallest KS value
    def get_current_KS():
        dict = {}
        temp = []
        for x in xrange(0,len(clusters)):
            for y in xrange(0,len(clusters)):
                if x != y:
                    # This is the simplest case.
                    pair = frozenset([clusters[x],clusters[y]])
                    if pair in all_KS:
                        if pair not in dict:
                            dict[pair] = all_KS[pair]
                        else:
                            continue
                    # Here for merged clusters. Put pairs in temporary list. 
                    else:
                        if pair not in temp:
                            temp.append(pair)
                        else:
                            continue
                else:
                    continue
        # For merged clusters, calculate the mean of all KS between clusters.
        for key in temp:
            pairs = []
            values = []
            # Values in first cluster of pair
            for x in xrange(0,len(cluster_defs[list(key)[0]])):
                # Values in second cluster of pair
                for y in xrange(0,len(cluster_defs[list(key)[1]])):
                    combi = frozenset([list(cluster_defs[list(key)[0]])[x],
                                   list(cluster_defs[list(key)[1]])[y]])
                    if combi not in pairs:
                        pairs.append(combi) 
                        values.append(float(all_KS[combi]))
            dict[key]=str(sum(values)/len(values))
            # If mode is M2 we update all_KS to treat subclusters as if they were terminal
            # sequences for next iteration.
            if redundant_kS_mode == 'M2':
               all_KS[key]=str(sum(values)/len(values))
        # If mode is M2 we update cluster definitions to treat subclusters as if they were
        # terminal sequences for next iteration.
        if redundant_kS_mode == 'M2':
           for i in xrange(0,len(clusters)):
              cluster_defs[clusters[i]]=[clusters[i]]
        return dict

    # get the min KS in current clustering
    def get_min_KS(ks_dict):
        min_KS = {}
        pair = min(ks_dict, key = ks_dict.get)
        min_KS[pair] = ks_dict[pair]
        # Append the KS score to the final list of KS scores.
        KS.append(ks_dict[pair])
        return min_KS

    # update the current clusters
    # having identified the smallest KS, those sequences/sub-clusters are merged and the 
    # cluster list modified.
    def update_clusters(cluster_list,min_KS):
        temp = []
        new_clusters = []
        for x in xrange(0,len(clusters)):
            for key in min_KS:
                # If KS value for a pair of sequences.
                if list(key)[0] in sequences and list(key)[1] in sequences:
                    cluster_defs["sub_clust_"+str(counter)]=list(key)
                    if clusters[x] in list(key):
                        temp.append(clusters[x])
                    else:
                        continue
                # If KS value for a sequence and a cluster of sequences.
                elif list(key)[0] in sequences or list(key)[1] in sequences:
                    seq_list = []
                    if list(key)[1] in sequences:
                        seq_list.append(str(list(key)[1]))
                        if len(cluster_defs[list(key)[0]]) > 1:
                            for y in xrange(0,len(cluster_defs[list(key)[0]])):
                                if cluster_defs[list(key)[0]][y] not in seq_list:
                                    seq_list.append(cluster_defs[list(key)[0]][y])
                    if list(key)[0] in sequences:
                        seq_list.append(str(list(key)[0]))
                        if len(cluster_defs[list(key)[1]]) > 1:
                            for y in xrange(0,len(cluster_defs[list(key)[1]])):
                                if cluster_defs[list(key)[1]][y] not in seq_list:
                                    seq_list.append(cluster_defs[list(key)[1]][y])
                    cluster_defs["sub_clust_"+str(counter)]=seq_list
                    if clusters[x] in list(key):
                        temp.append(clusters[x])
                    else:
                        continue
                # If KS value for two clusters of sequences.
                else:
                    seq_list = []
                    if len(cluster_defs[list(key)[0]]) > 1:
                        for y in xrange(0,len(cluster_defs[list(key)[0]])):
                            if cluster_defs[list(key)[0]][y] not in seq_list:
                                seq_list.append(cluster_defs[list(key)[0]][y])
                    if len(cluster_defs[list(key)[1]]) > 1:
                        for y in xrange(0,len(cluster_defs[list(key)[1]])):
                            if cluster_defs[list(key)[1]][y] not in seq_list:
                                seq_list.append(cluster_defs[list(key)[1]][y])
                    cluster_defs["sub_clust_"+str(counter)]=seq_list
                    if clusters[x] in list(key):
                        temp.append(clusters[x])
                    else:
                        continue
        for y in xrange(0,len(clusters)):
            if clusters[y] not in temp:
                new_clusters.append(clusters[y])
            else:
                if "sub_clust_"+str(counter) not in new_clusters:
                    new_clusters.append("sub_clust_"+str(counter))
        return new_clusters

    # An iterator to sequentially merge sequences/sub-clusters with shortest KS.
    while len(clusters) > 1:
        # if mode is M2 we update the sequence list so sub-clusters are treated as 
        # terminal sequences in next iteration.
        if redundant_kS_mode == 'M2':
           sequences = clusters
        counter = counter + 1
        current_KS = get_current_KS()
        min_KS = get_min_KS(current_KS)
        clusters = update_clusters(clusters,min_KS)
    return KS

##########################################################################################
#   Write the output file.
def write_output_file(output_list,data_dict):
   OUTFILE = open(output_prefix+"_ALL_KS.txt","w")   
   for xyz in output_list:
      print >>OUTFILE, xyz[0]+"\t"+xyz[1]+"\t"+xyz[2]+"\t"+xyz[3]
   OUTFILE.close()
   OUT1 = open(output_prefix+"_KS_by_cluster.txt", "w")
   OUT2 = open(output_prefix+"_REDUCED_KS.txt", "w")
   print "\nReducing redundant KS values, this will take a litte while...\n"
   for key in data_dict:
       KS = remove_redundant_KS(data_dict[key])
       print >> OUT1,key,KS
       for i in xrange(0,len(KS)):
          print >> OUT2,KS[i]
   OUT1.close()
   OUT2.close()

##########################################################################################
#   Write an R script to plot ks values.
def draw_plot():
   kS_scores_to_plot = str(output_prefix+"_REDUCED_KS.txt")
   RPLOT = open("rplot.R","w")
   working_directory = str(os.getcwd())
   print >> RPLOT,'data.dir <- "'+working_directory+'"'
   print >> RPLOT,'setwd(data.dir)'
   print >> RPLOT,'library(ggplot2)'
   print >> RPLOT,'M1 <- read.table("'+kS_scores_to_plot+'", header=FALSE)'
   print >> RPLOT,'ggplot() +'
   print >> RPLOT,'  ggtitle("kS plot") +'
   print >> RPLOT,'  xlab("KS") +'
   print >> RPLOT,'  coord_cartesian() +'
   print >> RPLOT,'  layer('
   print >> RPLOT,'    mapping=aes(x=M1$V1[M1$V1 <= 5]),'
   print >> RPLOT,'    stat="bin", stat_params=list(binwidth=0.05),'
   print >> RPLOT,'    geom="histogram", geom_params=list(alpha=0.5,fill="blue"))'
   print >> RPLOT,'ggsave(file="'+output_prefix+'kS_plot.pdf", width=4, height=4)'
   RPLOT.close()
   subprocess.call(["R", "CMD","BATCH","rplot.R"])   
   os.remove("rplot.R")
   os.remove("Rplots.pdf")
   os.remove("rplot.Rout")

##########################################################################################
#   Main function calls.
print "\nParsing blast table, this will take a little while..."
try:
   pairs = parse_blast(file_lines)
except:
   print >> sys.stderr, "\n\nERROR PARSING BLAST TABLE\n"
   BE = str('Blast output should be formatted using < -outfmt "7 std qlen slen" >\n\n')
   print >> sys.stderr,BE
   sys.exit()
print "\nBuilding blast clusters, this will take a little while..."
clusters = build_clusters(pairs)
print "\nGetting KS values, this takes quite a while...\n"
output_file_lines = get_KS(clusters)
data_dict = reduce_redundant_KS(output_file_lines)
write_output_file(output_file_lines,data_dict)
log_file.close()
try:
   draw_plot()
except:
   print >> sys.stderr, "\nERROR. COULD NOT DRAW PLOT\n"

if ERRORS >0:
   print "\n\nYour run completed with",ERRORS,"errors"
   print "Check log.txt for details\n\n"
else:
   print "\n\nYour run completed\n\n"
