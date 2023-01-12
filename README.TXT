#######Trips DB shelve#######
This shelve file will act as a central database for trips with one top level key for "file_ids",one top level key for "study_ids" one for "user_ids",and one key for "organism_ids". Should be stored
in the same place as "__init.py__" under the name trips_db.shelve. 
Each top  level key in "user_ids" will be an unique integer its values will be: "username" and "password", where password is a hash of the password plus a salt
Each top level key in "file_ids" will be an integer its values will be a dictionary with the following keys :"file_name", "file_description", "organism_id","study_id","file_type"(riboseq/rnaseq).
Each top level key in "study_ids" will be an unique integer with the following keys:"organism_id","GSE","Title","Year","Paper","Srp","Description","Adapters". 
Each top level key in "organism_id" will be an unique integer with the following keys: "organism_name","transcriptome_list", "gwips_databasename", "gwips_clade", "gwips_organism", "gwips_database", "default_transcript"

As shelve files are binary a seperate script is made to allow for easy manual modification of this file: modify_trips_db.py



#######0 or 1 based#######

Trips reads should be 0-based as currently using pysam to extract reads which always converts to 0-based regardless of input file (sam files are 1 based, bam files are 0 based) so should the extraction
tool change in future new reads should be 0 based. Co-ordinates in the "transcriptome shelve" file for things like exon junctions and cds_start/cds_stop should also be 0-based. 
Cds_start should correspond to the A of the AUG, cds_stop should correspond to the A/G of the stop codon.
E.g if the 5' leader seq and start codon was AGCTCGGCTAGCAUG the cds_start should be 12


This should also be the case on the frontend so position 0 on the plots should correspond to the first nucleotide in the sequence and the cds_start/stop markers should point to the first and last positions
of the start and stop codons. Exon junctions should point to the first nucleotide of the exon.

#######ADDING A NEW ORGANISM FOR TRIPS PROCESSING#########

1. Download a GTF file and fasta file.
    If you want pseudo utr's it's best to get the fasta file from UCSC (use microbes.ucsc.edu for bacteria or gwips) table browser by selecting CDS exons and add the desired length upstream/downstream
    Note: (When using table browser and selecting 300 bases upstream/downstream i had issues with some genes not showing up. Instead what I did was to use the option: One FASTA record per region
    (exon, intron, etc.) with extra bases upstream (5') and extra downstream (3'). This was ok as it was a bacterial genome that only had one entry per gene but would not work with some genomes)    

2. "Clean" the fasta file.
    This is done by running the script clean_fasta_headers.py in the ribopipe/scripts folder. This does a number of things to the headers such as removing anything after space or "-" characters
    replacing "(" and ")" with underscores,convert to uppercase, and optionally replaces a user defined term in the header (set as a variable at the top of the script). This is in an attempt to ensure that all the 
    headers will match those in the gtf file. However some manual intervention may still be required. The output of this will be a fasta file with the ".cleaned" extension appended to the filename. 
    This is what should be used from now on. 

3. Create the bowtie indices using the cleaned fasta file above.
    Using the bowtie-build command.

4. Create a "transcriptome shelve" file.
    This is just a dictionary of transcripts with info such as sequence, cds_start, cds stop, length, exon junctions. It is stored as a shelve object to make it easy to access later on. To create this run the script
    setup_transcriptome_shelve.py in the ribopipe/scripts directory. Pass it the following arguments in this order <organism> <gtf_file> <cleaned_fasta_file> <pseudo_utr_length>. Organism should have
    underscores in place of spaces. Make sure the pseudo utr length matches the pseudo utr length in step 1, if no pseudo utr's just enter 0. Change permissions to 771 before moving.

5. Create a "transcriptome to genome shelve".
    Another shelve object, this time used to convert transcriptomic to genomic co-ordinates (used to tell if a read is unambigouos). To create just run the script create_transcriptomic_to_genomic_dictionary.py
    in ribopipe/scripts passing the gtf file as the first argument, and passing the pseudo utr length as the second.  Change permissions to 771 before moving.

6. Add an entry to scripts/user_defined_info.py 
    In the trips_organism_files_dict add your organism as a key, the value is a list with 4 items in the following order: path to fasta file, path to gtf file, path to bowtie indices path to transcriptome shelve

7. In the same file as above add an entry to tran_to_genome_shelves
    Add the organism as a key and the path to the "transcriptome to genome shelve" as the value.



