#!/usr/bin/env python

import sys
import os
import argparse

'''
usage: frameshifts.py -t <table.tbl> -f <frameshift_folder>

correct frameshifts in table.tbl according to 
sqn.frameshifts.problems.short.txt

to obtain frameshift_folder upload sqn file to 
https://www.ncbi.nlm.nih.gov/genomes/frameshifts/frameshifts.cgi, download the results, save them and unzip them. 


by Isa
'''


def frameshifts(tbl,fs_long):
	
	desc_fs=open(fs_long)	#open submission check result file
	
	#make a dictionary of all frameshifts in the submission check result file 
	fs_dic={}
	
	for line in desc_fs:
		if line.startswith("gnl"):	#locus_tags of the two genes that are one pseudo gene
			locs=[]
			for item in line.split("\t"):
				locs.append(item.strip().split("|")[-1].split("(")[0])
		elif line.startswith("s"):	#annotation of the pseudo gene
			annotation=line.split("|")[-1]
		elif line.startswith("["):	#coordinates of the two genes that are one pseudo gene
			loc_2=max(locs) #save in dic with the higher locus_tag as key
			fs_dic[loc_2]=[locs,line.strip().split("\t\t"),annotation[1:].strip()]
	
	#change dic entries for when more than two genes are one pseudo gene
	key_last=0
	for key in sorted(fs_dic.keys()):
		if key_last!=0:
			if fs_dic[key][0][0]==fs_dic[key_last][0][1]:	#first locus_tag in key is second locus_tag in key_last
				fs_dic[key][1][0]=fs_dic[key_last][1][0]	#change first coordinates in key to first coordinates in key_last
			elif fs_dic[key][0][1]==fs_dic[key_last][0][0]:	#second locus_tag in key is first locus_tag in key_last
				fs_dic[key][1][1]=fs_dic[key_last][1][1]	#change second coordinates in key to second coordinates in key_last
			key_last=key
		else:
			key_last=key


	table=open(tbl) #open table file
	table_new=open(tbl.replace(".tbl","_frameshifts_fixed.tbl"),"w")	#open new table file
	
	#make list of entries to be printed to the new table file
	CDSs=[]
	fs="no"	#fs -> is current entry frameshift or not
	for line in table:
		#print line
		linesplit=line.split("\t")
		if line=="\n":
			print "Attention! Empty line in table file!"
		elif ">" in linesplit[0]: #save header as entry
			if "CDS" in locals():
				CDSs.append(CDS)
			CDS=[line]
		elif linesplit[0]!="" and linesplit[2]!="gene\n" and linesplit[2]!="sig_peptide\n":	#start of entry in the old table file
			CDSs.append(CDS)	#save previous entry to entries list
			CDS=[line]	#start new entry with first line
			fs="no"	#entry is not (yet known to be) a frameshift
			
		elif len(linesplit)>3 and linesplit[3]=="locus_tag" and linesplit[4].strip() in fs_dic.keys(): #entry is a frameshift
			fs="yes"	
			
			#since frameshifts were saved with the second locus_tag as key, now delete previous entry from entries list,
			#after getting it's locus_tag for the frameshift
			for item in CDSs[-1]:	
				for itemsplit in item.split("\n"):
					if "locus_tag" in itemsplit:
						loc=itemsplit.split("\t")[-1].strip()
			del CDSs[-1]
			
			key=linesplit[4].strip()	#key is locus_tag of current entry

			loc1=fs_dic[key][0][0]		#locus_tags of frameshift from dic
			loc2=fs_dic[key][0][1]
				
			nr1=int(loc1.split("_")[-1])	#the locus_tags as ints
			nr2=int(loc2.split("_")[-1])
			print nr1, nr2, loc1, loc2
			if nr1<nr2:	#coordinates for forward strand frameshift
				c1=str(int(fs_dic[key][1][0].split("...")[0].replace("[",""))+1)
				c2=str(int(fs_dic[key][1][1].split("...")[1].replace("]",""))+1)
			else:	#coordinates for reverse strand frameshift
				c2=str(int(fs_dic[key][1][1].split("...")[0].replace("[",""))+1)
				c1=str(int(fs_dic[key][1][0].split("...")[1].replace("]",""))+1)
			
			#make frameshift entry for new table file
			CDS=[c1+"\t"+c2+"\tgene\n\t\t\tpseudo\n\t\t\tlocus_tag\t"+loc+"\n\t\t\tgene_desc\t"+fs_dic[key][2]+"\n\t\t\tnote\tframeshift\n"]		
			
		elif fs=="no":
			CDS.append(line)
		else:
			pass
	
	# print new table file
	CDSs.append(CDS)
	for CDS in CDSs:
		for line in CDS:
			table_new.write(line)
	
	table_new.close()
			
def main(argv):

	parser = argparse.ArgumentParser()
	parser.add_argument('-t','--table', required=True)
	parser.add_argument('-f','--frameshifts', required=True)
	args = parser.parse_args()
	
	tbl=args.table	#the tbl file
	if args.frameshifts.endswith("/"):
		args.frameshifts=args.frameshifts[:-1]
	fs_long=args.frameshifts+"/"+args.frameshifts.split("/")[-1]+".sqn.frameshifts.problems.short.txt" #the frameshifts file

	frameshifts(tbl,fs_long)
	
if __name__ == "__main__":
   main(sys.argv)
