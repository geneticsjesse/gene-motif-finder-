# BINF_Assignment3.py
# Jesse Wolf 0830233

import sys
import re
import os
from Bio import SeqIO
import random 

# Setting the text file of our gene IDs to our first argument, the text file of our promoters to the second argument, and the folder containing the gff3 and fasta file to the third argument that we receive from the command line. In this script, we assume that only the gff3 and fasta files from chromosome 1 will be used.
geneIDs_file=sys.argv[1]
promoters_file=sys.argv[2]
gff3_fasta_folder=sys.argv[3]

# Assigning the path of this directory to the base_path variable and assigning the path to each file in this directory to their respective variables.
base_path = os.path.dirname(os.path.abspath(__file__))
geneID_file_path = os.path.join(base_path, geneIDs_file)
promoters_file_path = os.path.join(base_path, promoters_file)

# Assigning the path to the folder containing our gff3 and fasta file to the gff3fasta_folder_path variable.
gff3fasta_folder_path = os.path.join(base_path, gff3_fasta_folder)

 # Creating a list of our file paths to use in the subsequent function.
path_list = [geneID_file_path, promoters_file_path, gff3_fasta_folder]

# Creating a function to check our list of file paths.
def list_files_check(path_list):
    """ This function checks the availability of each file and asks for new user input if a file is absent."""
    # Create an empty list of paths to append to.
    available_paths = []
    # Create a for loop for each file in our existing path list.
    for file in path_list:
        # Create an if statement if a given file does not exist.
        if not os.path.exists(file):
            print("The file '" + file + "' is not available.")
            # Ask for new file(s).
            new_path_list = input("Please enter the file paths for the genes .txt file and promoters .txt file, as well as the path to the directory that contains the fasta file and gff3 file, separated by a comma and a space, in the following format: /path/to/genes_file, /path/to/promoters_file, /path/to/fasta_and_gff3_file: ")
            # Split the new path list by a comma and a space.
            path_list = new_path_list.split(", ")
            # Create a new variable equal to this function called on our initial path list.
            available_paths = list_files_check(path_list)
        else:
            # If all files are available, append them to the available paths list.
            available_paths.append(file)
    return available_paths

available_paths = list_files_check(path_list)        

# Using os.scandir to ensure we have two files, one of fasta and one of gff3 format and if we do, create a path for both.
with os.scandir(available_paths[2]) as entries:
    for entry in entries:
        if entry.is_file():
            if entry.name.endswith(".fa"):
                fasta_file_path = os.path.join(available_paths[2], entry.name)
            elif entry.name.endswith(".gff3"):
                gff3_file_path = os.path.join(available_paths[2], entry.name)

# Printing the genes text file, promoters text file, and path associated with the fasta and .gff3 file that are being used, separated by a line break.
print ("\n" + "Gene names were retrieved from the file " + available_paths[0] + "." + "\n")
print ("Promoter sequences were retrieved from the file " + available_paths[1] + "."+ "\n") 
print ("Fasta and .gff3 files were retrieved from the folder: " + available_paths[2] + "."+ "\n")

# Opening all relevant files.
genes = open (available_paths[0])
promoters = open (available_paths[1])
gff3 = open (gff3_file_path)
# Using BioPython to read in our fasta file.
fasta = SeqIO.read(fasta_file_path, "fasta")

# Reading each line into a list for the above opened files.
gff3_datalines = gff3.readlines()
promoters_datalines = promoters.readlines()
genes_datalines = genes.readlines()

# Creating an empty list to append our genes and promoters to.
promoters_list=[] 
genes_list=[]

# Creating a for loop for each promoter in our promoters file.
for i in promoters_datalines:
    # Removing line breaks from our promoters file and coercing all characters to uppercase, and appending to a list.
    promoters_list.append(i.rstrip("\n").upper())

# Creating a for loop for each gene in our genes file.
for i in genes_datalines:
    # Removing line breaks from our genes file and appending to a list.
    genes_list.append(i.rstrip("\n")) 

# Creating an empty dictionary for 1) All of the genes found in the provided gff3 file and 2) All of the genes within our genes_list that are also within the provided gff3 file.
gene_dict={}
gene_dict_subset={}

# Create a function that will create a dictionary containing gene names and their respective promoter sequences using our gff3 file. We use start and end sites for positive and negative strands, respectively, to obtain 400 nucleotides upstream/downstream as the promoter sequence. We also search each promoter sequence to locate any sequence of 100 N's, which represent an unknown length of missing DNA and if any such patterns are located, all nucleotides upstream are removed.

def dict_creator():
    """This function creates a dictionary from all genes found in our .gff3 file and a dictionary that matches gene names found in our fasta file with start sites of five prime UTRs in our gff3 file."""
# Creating a for loop for each line in our gff3 file.
    for line in gff3_datalines:
        # Ensuring we ignore lines that start with a #.
        if not line.startswith("#"):
            # Splitting lines using tab delimiters.
            gene_list=line.split("\t")
            # Creating an if statement to check if our third column is a five prime UTR.
            if gene_list[2]=="gene":
                # Using multiple splits within the ninth column to generate a gene name to match against our list of genes.
                gene_name = gene_list[8].split(":")[1].split(";")[0]
                # Assign the start site to the key of our dictionary.
                gene_dict[gene_name]=gene_list[3]
                # Creating an if statement to see if we are looking at the positive chromosome.
                if gene_list[6] == "+":
                    # If we are, we want the start of the sequence.
                    gene_pos_start = int(gene_list[3])
                    # Here, we set the key of a given gene to the 400 nucleotides before the start site.
                    gene_dict[gene_name] = fasta.seq[(gene_pos_start-400):gene_pos_start]
                    # Using re.split to remove any strings of 100 N's or more from our promoters.
                    N_matchs = re.split("(N{100,})", str(gene_dict[gene_name]))
                    # Replacing the value of a given gene name in our dictionary (our promoter sequences) with the last element of our re.split, which is only the nucleotides downstream of the N's.
                    gene_dict[gene_name]=N_matchs[-1]
                # Creating a for loop for each gene in list of genes.
                    for gene in genes_list:
                    # Creating an if statement to check if the gene name from our gff3 file matches the gene name in our list of genes.
                        if gene_name == gene:
                        # Creating an if statement to check if the gene name is already in our dictionary, if not, we add a key from our list of genes and the corresponding start site (the fourth column) of the 5 prime UTR.
                            if gene_name not in gene_dict_subset:
                                gene_dict_subset[gene_name]=fasta.seq[(gene_pos_start-400):gene_pos_start]
                                # Using re.split to remove any strings of 100 N's or more from our promoters.
                                N_matchs_subset = re.split("(N{100,})", str(gene_dict_subset[gene_name]))
                                # Replacing the value of a given gene name in our dictionary (our promoter sequences) with the last element of our re.split, which is only the nucleotides downstream of the N's.
                                gene_dict_subset[gene_name]=N_matchs_subset[-1]
                # Creating an if statement to see if we are looking at the negative chromosome.
                elif gene_list[6] =="-":
                    # If we are, we want the end of the sequence.
                    gene_neg_end = int(gene_list[4])
                    # Here, we set the key of a given gene to the reverse complement of the 400 nucleotides after the end site.
                    gene_dict[gene_name] = fasta.seq[(gene_neg_end):gene_neg_end+400].reverse_complement()
                    # Using re.split to remove any strings of 100 N's or more from our promoters.
                    N_matchs = re.split("(N{100,})", str(gene_dict[gene_name]))
                    # Replacing the value of a given gene name in our dictionary (our promoter sequences) with the last element of our re.split, which is only the nucleotides downstream of the N's.
                    gene_dict[gene_name]=N_matchs[-1]
                    # Creating a for loop for each gene in list of genes.
                    for gene in genes_list:
                    # Creating an if statement to check if the gene name from our gff3 file matches the gene name in our list of genes.
                        if gene_name == gene:
                        # Creating an if statement to check if the gene name is already in our dictionary, if not, we add a key from our list of genes and the corresponding start site.
                            if gene_name not in gene_dict_subset:
                                gene_dict_subset[gene_name]=fasta.seq[(gene_neg_end):gene_neg_end+400].reverse_complement()
                                # Using re.split to remove any strings of 100 N's or more from our promoters.
                                N_matchs_subset = re.split("(N{100,})", str(gene_dict_subset[gene_name]))
                                # Replacing the value of a given gene name in our dictionary (our promoter sequences) with the last element of our re.split, which is only the nucleotides downstream of the N's.
                                gene_dict_subset[gene_name]=N_matchs_subset[-1]    
    return [gene_dict, gene_dict_subset]
dict_creator()

# We can check how many genes we have in our gff3 file relative to the genes that are found in both the Zea mays .gff3 and genes file.
print ("There are " + str((len(gene_dict))) + " genes that are found in the Zea mays .gff3 file relative to " + str(len(gene_dict_subset))+ " genes that are found in both the Zea mays .gff3 and genes .txt file." + "\n")

# Using BLAST to check if our promoter sequences show up in Zea Mays. The below lines are commented out to avoid being run each time.
# We randomly subsample 5 of our promoters to use in the BLAST command to avoid a lengthy compute time/crashing.
# rand_sample = random.sample(list(gene_dict_subset.values()), 5)
# from Bio.Blast import NCBIWWW
# result_handle = NCBIWWW.qblast("blastn", "nt", rand_sample)
# with open('blast_results.xml', 'w') as save_file: 
#     blast_results = result_handle.read() 
#     save_file.write (blast_results)

# From the BLAST results, we see that a subset of the promoter sequences from our dictionary of known Zea mays genes are in fact, identified as Zea Mays.

# Now, we want to randomnly subsample an equivalent number (n=125) of genes in 5 samples. We can use the length of our subset dictionary to determine the sample size for each subsampling event.
# Creating a variable for the number of random samples we want. This can be changed in future iterations for increased reproducibility.
num_randoms = 5

# Creating a function that randomly samples gene names and adds them to a list.
def random_sampler():
    """ This function randomly samples gene names from our large dictionary and appends them to a master list that includes our subset dictionary """
    # Creating an empty list to append to.
    list_of_lists=[]
    # Appending our subset dictionary to the list.
    list_of_lists.append(gene_dict_subset)
    # Creating a for loop for each integer in the range of our num_randoms variable.
    for i in range(num_randoms):
        # Creating a variable that represents a random sample from our large dictionary, specifying the length as the length of our subset dictionary.
        random_samples=random.sample(list(gene_dict), len(gene_dict_subset))
        # Appending our 5 random sample list to our larger list.
        list_of_lists.append(random_samples)
    return list_of_lists
random_sampler()

# Create a function to count motifs.
def motif_finder (lists):
    """ This function counts the number of times a given motif is found within a promoter region """
    # Open a text file so that we can write to it.
    with open("motif_counts_output.txt","w") as x:
        # Creating a for loop for each promoter in our promoters list.
        for promoter in promoters_list:
            # Writing to our text file with the name of each promoter and using rstrip to remove the new line characters.
            x.write(promoter.rstrip("\n"))
            # Creating a for loop that uses each list in our larger list of lists.
            for list in lists:
                # Creating a variable for motif counts that we can add to within the subsequent for loop.
                motif_counts = 0
                # Creating a for loop that looks at each promoter in a given list.
                for i in list:
                    # Creating a variable that takes the value of the promoter from our dictionary.
                    promoter_sequence = gene_dict[i]
                    # Use .upper to coerce all motifs to capitals and rstrip to remove the new line characters.
                    promoter_motif = promoter.rstrip("\n").upper()
                    # Using re.findall to find all matches of the motifs provided, within each promoter sequence.
                    matchs = re.findall(promoter_motif, str(promoter_sequence))
                    # Creating a variable so that we can count the number of times a given motif was found within our promoter.
                    motif_counts += len(matchs)
                # Writing to our text file with the count of a given motif delimited by tabs.
                x.write("\t" + str(motif_counts))
            x.write ("\n")
# Calling our function using our random_sampler function as the argument.
motif_finder(random_sampler())

# Now that we have the average counts of our random samples, we can determine which motifs are over or under-expressed relative to our known Zea mays genes.
#The best way to determine if a specific motif is over or under-expressed is to compare the expression levels of the genes associated with the motif to the expression levels of the overall gene population. If the expression levels of genes associated with the motif are significantly higher or lower than the expression levels of the overall gene population, the motif can be considered to be over or under-expressed, respectively.

# Create a variable to add to within our for loop.
genes_subset=0
# Create a for loop to iterate over each gene found in the genes_list that was created from our gene file.
for gene in genes_list:
    # If the gene from our text file is also in our dictionary, add one to the value of our genes_subset variable.
    if gene in gene_dict:
        genes_subset+=1

# Create a function that looks at differential expression between our known Zea mays genes and our randomly sampled genes.
def expression_check(input):
    """ This function counts the number of times a given motif is found within a promoter region """
    # Open a file that we can write to.
    with open('differential_expression_output.txt', 'w') as outfile:
        # Take an input file and open it. In this case, our input file is a file of motif counts.
        with open (input, 'r') as infile:
            # Write column names for file accessibility.
            outfile.write ("Motif" + "\t" + "Zea_mays_count" + "\t" + "Random_count" + "\t" + "Relative_expression" + "\n")
            # Creating a for loop for each line of our input file.
            for line in infile:
                # Strip the line endings and split by tabs.
                line_list = line.rstrip("\n").split("\t")
                # Create a variable that ignores the first two columns of our input file.
                random_counts = line_list[2:]
                # Creating a variable that only takes the motif count from our known genes.
                gene_value = int(line_list[1])
                # Calculating the average count across all random samples.
                random_expression_avg = (sum(map(int, random_counts))/len(random_counts))
                # Creating an if statement that looks at all random same counts above 0.
                if random_expression_avg>0:
                    # Checking if the count of our known Zea mays gene is greater than 1.25x than our random gene count and if so, we denote this motif as over-expressed.
                    if gene_value > (1.25*random_expression_avg):
                        relative_expression = "This motif is over-expressed"
                    # Checking if the count of our known Zea mays gene is less than 0.75x than our random gene count and if so, we denote this motif as under-expressed.
                    elif gene_value < (0.75*random_expression_avg):
                        relative_expression = "This motif is under-expressed"
                    # If we don't satisfy either of the above if statements, we denote this motif as not significantly differentially expressed.
                    else:
                        relative_expression = "This motif is not significantly differentially expressed"
                    # Write to our outfile: the motif sequence, the count from our known Zea mays genes, the count from our randomly sampled genes, and whether the gene is significantly differentially expressed or not.
                    outfile.write(line_list[0] + "\t" + line_list[1] + "\t" + str(random_expression_avg) + '\t' + relative_expression + "\n")
                # Creating an if statement that looks at all random counts that are equal to 0.
                elif random_expression_avg==0:
                    # Create a threshold value from the number of a given gene within our subset dictionary.
                    threshold = (genes_subset*0.1)
                    # Checking to see if the count of our known Zea mays genes is above the threshold and if it is, the motif is over-expressed.
                    if gene_value > threshold:
                        relative_expression = "This motif is over-expressed"
                    else:
                        relative_expression = "This motif is not significantly differentially expressed"
                    # Write to our outfile: the motif sequence, the count from our known Zea mays genes, the count from our randomly sampled genes, and whether the gene is significantly differentially expressed or not.
                    outfile.write(line_list[0] + "\t" + line_list[1] + "\t" + str(random_expression_avg) + '\t' + relative_expression + "\n")
        # Return our output file.
        return (outfile)
# Running our function on the motif_counts.txt file. The below is commented out to avoid re-running every time the script is run.
expression_check("motif_counts_output.txt")

print("The analysis is now complete! You should have two new text files located in your working directory. If you ran the BLAST commands, you will have a third file (of .xml format).")
print("The first text file contains the motif counts for the genes that appear in " + geneIDs_file + " and 5 random samples.")
print("The second text file contains motifs, their counts, and whether or not they are differentially expressed between the selected and randomly sampled genes." + "\n")