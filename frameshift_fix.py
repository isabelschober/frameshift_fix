#!/usr/bin/env python

import sys
import os
import argparse

'''
usage: frameshifts.py -t <table.tbl> -f <frameshift_folder>


Correct frameshifts in table.tbl according to 
sqn.frameshifts.problems.short.txt

To obtain frameshift_folder upload sqn file to 
https://www.ncbi.nlm.nih.gov/genomes/frameshifts/frameshifts.cgi, download the results, save them and unzip them. 


by Isabel
'''


def frameshifts_list(tbl,fs_long):
	#store all identified frameshifts
	
	desc_fs=open(fs_long)	#open submission check result file
	del_list=[] #list of locus_tags of entries that will later be deleted
	
	#make a dictionary of all frameshifts in the submission check result file 
	fs_dic={}
	
	for line in desc_fs:
		
		if line.startswith("gnl"):	#locus_tags and strand of the two genes that are one pseudo gene
			locs=[]
			for item in line.strip().split("\t"):
				locs.append(item.strip().split("|")[-1].split("(")[0])
				strand=item[-2]
		
		elif line.startswith("s"):	#annotation of the pseudo gene
			annotation=line.split("|")[-1].strip()
		
		elif line.startswith("["):	#coordinates of the two genes that are one pseudo gene
			coo=[]
			for item in line.strip().split("\t\t"):
				for coordinate in item.split("..."):
					coo.append(int(coordinate.replace("[","").replace("]","")))
			
			#save the frameshift
			loc_keep=min(locs) #the locus_tag of the pseudogene
			
			#check if there are additional entries between the two frameshiften genes
			nr1=int(locs[0].split("_")[1])	#the locus_tags as ints
			nr2=int(locs[1].split("_")[1])
			ident=locs[0].split("_")[0] #the letter part of the locus_tags
			nr_len=len(locs[0].split("_")[1]) #the number of digits in the locus_tags
			if abs(nr1-nr2)!=1: #if there are additional entries, add them to the list of locus_tags to be deleted later
				locs=[]
				for i in range(min(nr1,nr2),max(nr1,nr2)+1):
					locs.append(ident+"_"+str(i).zfill(nr_len))
						
			locs.remove(loc_keep)
			fs_dic[loc_keep]=[locs,coo,annotation, strand] #store the locus to keep, the loci to delete, the coordinates and the annotation and the strand
			for loc in locs:
				if loc not in del_list:
					del_list.append(loc) #store the loci to delete additionally in the delete list
	
	
	#change dic entries for when more than two genes are one pseudogene
	for key in sorted(fs_dic.keys(),reverse=True):
		if key in del_list: #if a locus of a pseudogene is additionally in the delete list -> merge two pseudogenes
			for key1 in fs_dic.keys():
				if key in fs_dic[key1][0]: #find the second pseudogene locus that has the other locus in its delete list
					for coordinate in fs_dic[key][1]:
						if coordinate not in fs_dic[key1][1]:
							fs_dic[key1][1].append(coordinate) #add the coordinates of the other locus to the pseudogene locus that has the other locus in its delete list
					for loc in fs_dic[key][0]:
						if loc not in fs_dic[key1][0]:
							fs_dic[key1][0].append(loc) #add additonal loci to the delete list
			del fs_dic[key] #delete the now redundant second pseudogene

	return fs_dic, del_list
	
def frameshifts_fix(entries, fs_dic, del_list):
	#delete the entries of frameshifted genes and add pseudogene entries
	
	for i in range(len(entries)-1,-1,-1): #go through all feature entries of the original table file !backwards! because we will delete entries
		loc=""
		feature=entries[i][0].split("\t")[-1].replace("\n","") #type of feature
		for item in entries[i]:
			if "locus_tag" in item:
				loc=item.split("\t")[-1].strip() #get locus_tag of entry
		if loc in del_list:
			del entries[i] #delete the entry if it is in the list of entries to be deleted
		elif loc in fs_dic.keys() and feature=="CDS":
			del entries[i] #also delete CDS entries of the pseudogenes, for which we will only keep the gene entries
		elif loc in fs_dic.keys(): #instead of the original gene entry, construct a new pseudogene entry
			c1=min(fs_dic[loc][1])+1 #coordinates are the lowest and highest coordinates from the list of genes that are being merged to a pseudogene
			c2=max(fs_dic[loc][1])+1
			if fs_dic[loc][3]=="+": #write pseudogene entry on forward strand
				entries[i]=[str(c1)+"\t"+str(c2)+"\tgene\n\t\t\tpseudo\n\t\t\tlocus_tag\t"+loc+"\n\t\t\tnote\t"+"frameshift; "+fs_dic[loc][2]+"\n"]	
			elif fs_dic[loc][3]=="-": #write pseudogene entry on reverse strand
				entries[i]=[str(c2)+"\t"+str(c1)+"\tgene\n\t\t\tpseudo\n\t\t\tlocus_tag\t"+loc+"\n\t\t\tnote\t"+"frameshift; "+fs_dic[loc][2]+"\n"]	
			print("Corrected "+loc, end="") #print out the locus_tags of entries that changes were made to
			for del_loc in fs_dic[loc][0]:
				print(" "+del_loc, end="")
			print()
		else:
			pass

	return entries #corrected list of entries
			
def main(argv):

	parser = argparse.ArgumentParser(description="Correct frameshifts in table.tbl according to sqn.frameshifts.problems.short.txt\n\nTo obtain frameshift_folder upload sqn file to https://www.ncbi.nlm.nih.gov/genomes/frameshifts/frameshifts.cgi, download the results, save them and unzip them. \n\nby Isabel", formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-t','--table', required=True)
	parser.add_argument('-f','--frameshifts', required=True)
	args = parser.parse_args()
	
	print("\nframeshift_fix v.2\n")

	tbl=args.table	#the input tbl file
	if args.frameshifts.endswith("/"):
		args.frameshifts=args.frameshifts[:-1]
	fs_long=args.frameshifts+"/"+args.frameshifts.split("/")[-1]+".sqn.frameshifts.problems.short.txt" #the frameshifts file

	fs_dic, del_list=frameshifts_list(tbl,fs_long) #store all identified frameshifts

	table=open(tbl) #open input tbl file
	table_new=open(tbl.replace(".tbl","_frameshifts_fixed.tbl"),"w")	#open new output table file
	
	#make list of entries in the old table file
	entries=[]

	for line in table:
		linesplit=line.split("\t")
		
		if line=="\n":
			print("Attention! Empty line in table file!")
			
		elif line.startswith(">"): #save header as separate entry
			if "entry" in locals(): #if this is not the first header in the file, save the previous entry
				entries.append(entry)
			entry=[line]	
			
		elif linesplit[0]!="":	#start of entry in the old table file
			entries.append(entry)	#save previous entry to entries list
			entry=[line]	#start new entry with first line
		
		else:
			entry.append(line)
			
	entries.append(entry) #save last entry
	
	entries=frameshifts_fix(entries, fs_dic, del_list) #correct frameshift in list of entries
	
	#print new table file
	for entry in entries:
		for line in entry:
			table_new.write(line)
	table_new.close()
	
	print("\ndone!\n")
	
if __name__ == "__main__":
   main(sys.argv)
