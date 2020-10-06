
# David L Patton


import re #regex import
import sys #import system tools
import numpy 
import csv
import pandas
from Bio.Seq import Seq #import biopython

#Read In Fasta File and Convert the string to a Char list to operate on
my_file = open("Meso.fasta")
file_contents = my_file.read().replace('\n', '')
fasta_array = list(file_contents)
my_file.close()


#Read in the WT BaseCall File
#++++++++++++++++++++++++++++++
reader = csv.reader(open("Mesoplasma_Florum_WT.csv", "rb"), delimiter=",")
x = list(reader)
BaseCall = numpy.array(x)
#print(BaseCall)
#++++++++++++++++++++++++++++++

seqlength = len(fasta_array)
column_names = ['T', 'G', 'C', 'A']
row_names    = ['TxT', 'TxG', 'TxC', 'TxA', 'GxT', 'GxG', 'GxC', 'GxA', 'CxT', 'CxG', 'CxC', 'CxA', 'AxT', 'AxG', 'AxC', 'AxA']
ori = 0
term = 400000
#print BaseCall[0][:] #Access individual row in Basecall File

Neighbor_Array_Chromosome = numpy.zeros((16,4))
Neighbor_Array_left = numpy.zeros((16,4))
Neighbor_Array_right = numpy.zeros((16,4))




#Loop to take the BaseCall file and retrieve the left and right neighbors of each mutation
#output_file = open("Meso_output_chromosome.txt", "w")
#output_file_left = open("Meso_output_left.txt", "w")
#output_file_right = open("Meso_output_right.txt", "w")
output_file_diag_pos = open("Meso_output_Diag.txt", "w")
for row in range(BaseCall.shape[0]):
	pos = BaseCall[row][0]	#Grabs position in each row
	refpos = int(BaseCall[row][0])	#coerces position to int for retrieval purposes
	ConcNuc = BaseCall[row][1]	#concensus nucleotide
	MutNuc = BaseCall[row][2]	#mutant nucleotide
	concfasta = fasta_array[refpos]
	#refseqnuc = fasta_array[refpos]	#reference nucleotide for diagnostics
	left = fasta_array[refpos-1]	#left codon nucleotide
	right = fasta_array[refpos+1]	#right codon nucleotide
	left_test = fasta_array[refpos-2] #correct
	right_test = fasta_array[refpos]  #correct
	concfasta_test = fasta_array[refpos-1] #correct

 	output_file_diag_pos.write('Position in Basecall file: ' + str(pos) + ' Position in fasta(index is n-1): ' + str(refpos)+ '\n') 
	output_file_diag_pos.write('Consensus(file): ' + ConcNuc+ ' '+ 'Consensus(Fasta): '+ concfasta_test + '\n')
	output_file_diag_pos.write('Default Codon: ' +str(left)+str(concfasta)+str(right)+'   New codon method: ' +str(left_test)+str(concfasta_test)+str(right_test)+ '\n\n')

