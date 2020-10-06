# Test
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
output_file = open("Meso_output_chromosome.txt", "w")
output_file_left = open("Meso_output_left.txt", "w")
output_file_right = open("Meso_output_right.txt", "w")
output_file_diag_pos = open("Meso_output_Diag.txt", "w")
for row in range(BaseCall.shape[0]):
	pos = BaseCall[row][0]	#Grabs position in each row
	refpos = int(BaseCall[row][0])	#coerces position to int for retrieval purposes
	ConcNuc = BaseCall[row][1]	#concensus nucleotide
	MutNuc = BaseCall[row][2]	#mutant nucleotide

	left = fasta_array[refpos-2] #correct
	right = fasta_array[refpos]  #correct
	concfasta = fasta_array[refpos-1] #correct

	left_test = fasta_array[refpos-1]
	right_test = fasta_array[refpos+1]
	concfasta_test = fasta_array[refpos]
	#left_test = fasta_array[refpos-2]
	#right_test = fasta_array[refpos]
	#concfasta_test = fasta_array[refpos-1]

	#left_test = fasta_array[refpos-1]
	#right_test = fasta_array[refpos+1]
	#concfasta_test = fasta_array[refpos]
	if(concfasta != ConcNuc):
		print "DISCREPANCY DETECTED"
		output_file_diag_pos.write("DISCREPANCY DETECTED AT: " + refpos + '\n')

 	output_file_diag_pos.write('Position in Basecall file: ' + str(pos) + ' Position in fasta(index is n-1): ' + str(refpos)+ '\n') 
	output_file_diag_pos.write('Consensus(file): ' + ConcNuc+ ' '+ 'Consensus(Fasta): '+ concfasta + '\n')
	output_file_diag_pos.write('Current Codon Method (grabs correct consensus, incorrect neighbors): ' +str(left)+str(concfasta)+str(right)+'   New Codon Method: ' +str(left_test)+str(concfasta_test)+str(right_test)+ '\n\n')
	#output_file_diag_pos.write('Default Codon: ' +str(left)+str(concfasta)+str(right)+'   New codon method: ' +str(left_test)+str(concfasta_test)+str(right_test)+ '\n\n')
	#output_file_diag_pos.write('Default Codon: ' + left+concfasta+right+'   New codon method: ' +left_test+concfasta_test+right_test)

	#write the output chromosome information
	#output_file.write('Position: ' + pos + ' Left: '+ left + ' Consensus: ' + ConcNuc +  ' Mutant: ' + MutNuc + ' Right: ' + right + '\n')

	rowswitch = left+right
	colswitch = MutNuc
	#print 'row: ' + rowswitch
	#print 'col: ' + colswitch
	rowiter = ''
	coliter = ''
	
	#Switch to determine Column
	if colswitch == 'T':
		coliter = 0
	elif colswitch == 'G':
		coliter = 1
	elif colswitch == 'C':
		coliter = 2
	elif colswitch == 'A':
		coliter = 3
	
	#Switch to determine the row
	if rowswitch == 'TT':
		rowiter = 0
	elif rowswitch == 'TG':
		rowiter = 1
	elif rowswitch == 'TC':
		rowiter = 2
	elif rowswitch == 'TA':
		rowiter = 3
	elif rowswitch == 'GT':
		rowiter = 4
	elif rowswitch == 'GG':
		rowiter = 5
	elif rowswitch == 'GC':
		rowiter = 6
	elif rowswitch == 'GA':
		rowiter = 7
	elif rowswitch == 'CT':
		rowiter = 8
	elif rowswitch == 'CG':
		rowiter = 9
	elif rowswitch == 'CC':
		rowiter = 10
	elif rowswitch == 'CA':
		rowiter = 11
	elif rowswitch == 'AT':
		rowiter = 12
	elif rowswitch == 'AG':
		rowiter = 13
	elif rowswitch == 'AC':
		rowiter = 14
	elif rowswitch == 'AA':
		rowiter = 15


	#print 'rowiter: ' + str(rowiter)
	#print 'coliter: ' + str(coliter)
	#coliter = int(coliter)
	#rowiter = int(rowiter)
	increment = Neighbor_Array_Chromosome[rowiter][coliter]
	increment = increment + 1
	Neighbor_Array_Chromosome[rowiter][coliter] = increment

	if refpos >= ori and refpos <= term:
		increment = Neighbor_Array_right[rowiter][coliter]
		increment = increment + 1
		Neighbor_Array_right[rowiter][coliter] = increment
		output_file_right.write('Position: ' + pos + ' Left: '+ left + ' Consensus: ' + ConcNuc +  ' Mutant: ' + MutNuc + ' Right: ' + right + '\n')

	if refpos >= term or refpos <= ori:
		increment = Neighbor_Array_left[rowiter][coliter]
		increment = increment + 1
		Neighbor_Array_left[rowiter][coliter] = increment
		output_file_left.write('Position: ' + pos + ' Left: '+ left + ' Consensus: ' + ConcNuc +  ' Mutant: ' + MutNuc + ' Right: ' + right + '\n')


Array_Output = pandas.DataFrame(Neighbor_Array_Chromosome, columns=column_names, index=row_names)
Array_Output_left = pandas.DataFrame(Neighbor_Array_left, columns=column_names, index=row_names)
Array_Output_right = pandas.DataFrame(Neighbor_Array_right, columns=column_names, index=row_names)
#Array_file = open("Mesoplasma_Florum_WT_Chromosome.csv", "w")
Array_Output.to_csv("Mesoplasma_Florum_WT_Chromosome.csv", sep = ",")
Array_Output_left.to_csv("Mesoplasma_Florum_WT_left.csv", sep = ",")
Array_Output_right.to_csv("Mesoplasma_Florum_WT_right.csv", sep = ",")

output_file.close()
output_file_left.close()
output_file_right.close()
#sum_elements = sum(Neighbor_Array_Chromosome)
#print Array_Output

#print 'Sum of all recorded mutants in matrix: '
#print sum_elements





#==================================
#Scratch code:

#Check if Nucleotides are equal at given position
#if refseqnuc == ConcNuc:
#	print 'You have a match!'
#	print 'reference: '+ refseqnuc + ' ' + 'Basecall: ' + ConcNuc 
#else:
#	print 'error'
#	print  'Position: ' + refpos + 'reference: '+ refseqnuc + ' ' + 'Basecall: ' + ConcNuc