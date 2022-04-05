# frameshift_fix
Fixes pseudogenes/frameshifts detected by NCBI's 'Microbial Genome Submission Check' in feature table files of assembled and annotated genomes.


Frameshifted pseudogenes occur when a SNP introduces a frameshift, which leads to an early stop codon. While we then assume these to be non-functional pseudogenes, in automatic gene prediction and annotation, this often leads to two shorter gene(-fragment)s to be called and annotated with the same function.

The NCBI developed an online tool called **'Microbial Genome Submission Check'**, which finds consecutive genes that should actually be combined as pseudogenes. It can be found at https://www.ncbi.nlm.nih.gov/genomes/frameshifts/frameshifts.cgi.

## Preparation

This tool currently only takes single contigs as input, so we usually only use it for the chromosome. If your genome contains plasmids, you have to separate the chromosome from them before uploading your genome to the website. To do that extract the entries corresponding to your genome's chromosome in the fasta and feature table files and save them seprately, e.g. as **<your_genome_chr.fsa>** and **<your_genome_chr.tbl>**. 

To create the file needed for submission to the 'Microbial Genome Submission Check', use NCBI's [**tbl2asn**](https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/) tool on your chromosome-only files. 
```bash
tbl2asn -V b -a s -Z <your_genome_chr.err> -i <your_genome_chr.fsa>
```
This outputs (among others), a file named your_genome_chr.sqn, which can now be uploaded at https://www.ncbi.nlm.nih.gov/genomes/frameshifts/frameshifts.cgi. Keep the tick for 'Run slow checks' and click 'Check submission' and stay on the page or copy the Job ID/ bookmark the page. When it is done, download the zipped folder containing the results into the folder that contains your genome files and extract it (<submission_check_results_folder>).

## Usage

**Attention:** Use the feature table file of the complete genome (your_genome.tbl) in this command, not the chromosome-only file. (The chromosome-only files are not needed anymore)

```bash
frameshift_fix.py -t <your_genome.tbl> -f <submission_check_results_folder>
```

## Output

The script outputs a file named **<your_genome_frameshifts_fixed.tbl>**. Delete the previous your_genome.tbl file or rename it to keep it. Then delete the 'frameshifts_fixed' part from the name of your_genome_frameshifts_fixed.tbl and run tbl2asn to transfer the changes to the GenBank and Sequin files. 
