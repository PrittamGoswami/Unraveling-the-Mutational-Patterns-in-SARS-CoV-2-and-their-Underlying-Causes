import os
import subprocess





def run_mafft(input_fasta, output_fasta, threads=3):
    """
    Runs MAFFT alignment on a given FASTA file using a specified number of threads.

    Parameters:
    input_fasta (str): Path to the input FASTA file.
    output_fasta (str): Path where the aligned output will be saved.
    threads (int): Number of threads to use (default is 3).
    """
    # MAFFT command with the specified number of threads
    command = ["mafft", "--thread", str(threads), input_fasta]

    try:
        # Open the output file and run the MAFFT command
        with open(output_fasta, "w") as output_file:
            subprocess.run(command, stdout=output_file, check=True)
        print(f"MAFFT alignment completed successfully. Output saved to {output_fasta}")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during MAFFT alignment: {e}")






def create_phylogenetic_tree(alignment_file, outgroup_sequence, sample, threads=3):
    """
    Create a rooted IQ-TREE phylogenetic tree by specifying the first sequence in the alignment as the outgroup.

    Parameters:
        alignment_file (str): Path to the input alignment file.
        outgroup_sequence (str): Identifier of the outgroup sequence in the alignment.
        sample (str) : Denotes the sample number 
        threads (int): Number of threads to use (default is 3).

    Returns:
        str: Path to the output tree file containing the rooted phylogenetic tree.
    """
    # Define the IQ-TREE command with the outgroup specified
    iqtree_cmd = [
        "iqtree2",
        "-s", alignment_file,
        "-o", outgroup_sequence,  # Specify the outgroup sequence
        "-m", "GTR+F+I+G4",  # Substitution model 
        "-nt", str(threads),  # Number of threads 
        "-redo",  # Redo the analysis from scratch if the results already exists
        "-fast",  # use fasttree contrusction algorithm
        "-asr",  # ancestral state reconstruction
        "-pre", f"../Samples/{sample}/{sample}"
    ]

    # Execute the IQ-TREE command
    try:
        subprocess.run(iqtree_cmd, check=True)
        print("IQ-TREE analysis completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running IQ-TREE: {sample}", e)


def construct_msa_and_phylogenetic_tree(folder_path, outgroup_sequence, sample_size, threads_to_be_used):
    # Iterate over all samples in the directory
    for sample in os.listdir(folder_path):
        print(sample)
        # Create the full path for each item
        input_fasta = os.path.join(folder_path, sample, f"SARS-CoV-2_{sample}_{sample_size}+1.fasta")
        # construct MSA output file path
        msa_output_fasta = os.path.join(folder_path, sample, f"SARS-CoV-2_{sample}_{sample_size}+1_msa.fasta")
        run_mafft(input_fasta, msa_output_fasta, threads=threads_to_be_used)
        # construct phylogenetic tree        
        create_phylogenetic_tree(msa_output_fasta, outgroup_sequence, sample, threads_to_be_used)

# Run the Code
sample_size=3600
outgroup_sequence = "NC_045512.2"
threads_to_be_used=30        
construct_msa_and_phylogenetic_tree("../Samples", outgroup_sequence,sample_size, threads_to_be_used)
