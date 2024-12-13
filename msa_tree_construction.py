import os
import subprocess




def run_mafft(input_fasta, output_fasta, threads=3):
    """
    Runs MAFFT alignment on a given FASTA file  
    with a specified number of threads.

    Parameters:
        - input_fasta (str): Path to the input FASTA file.
        - output_fasta (str): Path where the aligned output will be saved.
        - threads (int): Number of threads to use (default is 3).
    """
    # MAFFT command with specified number of threads
    command = ["mafft", "--thread", str(threads), input_fasta]

    try:
        # Open the output file and run the MAFFT command
        with open(output_fasta, "w") as output_file:
            subprocess.run(command, stdout=output_file, check=True)
        print(f"MAFFT alignment completed successfully. Output saved to {output_fasta}")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during MAFFT alignment: {e}")




def create_phylogenetic_tree(alignment_file, outgroup_sequence, sample, threads,):
    """
    Create a phylogenetic tree using IQ-TREE by specifying the first sequence in the alignment as the outgroup.

    Parameters:
        - alignment_file (str): Path to the input alignment file.
        - outgroup_sequence (str): Identifier of the outgroup sequence in the alignment.
        - sample (str): Depicts the Sample number to which the alignment file belongs.
        - threads (int): Number of threads to use for running the IQTREE.
    """
    # Define the IQ-TREE command with the outgroup specified
    iqtree_cmd = [
        "iqtree2",
        "-s", alignment_file,
        "-o", outgroup_sequence,  # the outgroup sequence
        "-m", "GTR+F+I+G4",  # Substitution model
        "-nt", str(threads),  # Number of threads to use
        "-redo",  # Redo the analysis from scratch if it already exists
        "-fast",  # Use the fasttree algorithm for tree construction
        "-asr",  # Ancestral state reconstruction
        "-pre", f"../Samples/{sample}/{sample}" # Location to save the generated files
    ]

    # Execute the IQ-TREE command
    try:
        subprocess.run(iqtree_cmd, check=True)
        print("IQ-TREE analysis completed successfully.")
    except subprocess.CalledProcessError as e:
        print("Error running IQ-TREE:", e)



def construct_msa_and_phylogenetic_tree(folder_path, outgroup_sequence, sample_size, threads_to_be_used):
    """
    Construct the Multiple Sequence Alignment and the Phylogenetic tree for all the sampled genome sets

    Parameters:
        - folder_path (str): The path to the folder containing all the sampled genomesets sub-directories
        - outgroup_sequence (str): Identifier of the outgroup sequence.
        - sample_size (str): Number of Sampled genomes in each Genome set
        - threads_to_be_used (int): Number of threads to use for running the IQTREE.
    """
    
    
    # Iterate over all samples in the directory
    for sample in os.listdir(folder_path):
        print(sample)
        # Designate the full path for each multifasta file corresponding to each sampled genomeset
        input_fasta = os.path.join(folder_path, sample, f"SARS-CoV-2_{sample}_{sample_size}+1.fasta")
        # construct MSA output file path
        msa_output_fasta = os.path.join(folder_path, sample, f"SARS-CoV-2_{sample}_{sample_size}+1_msa.fasta")
        run_mafft(input_fasta, msa_output_fasta, threads=threads_to_be_used)
        # construct phylogenetic tree        
        create_phylogenetic_tree(msa_output_fasta, outgroup_sequence, sample, threads_to_be_used)


sample_size=4000 # No of Genomes per sample
outgroup_sequence = "NC_045512.2" # Outgroup Sequence ID
threads_to_be_used=30 # No of CPU threads to use       

# Run the Code to Generate MSA and Phylogenetic trees
construct_msa_and_phylogenetic_tree("../Samples/", outgroup_sequence,sample_size, threads_to_be_used)
