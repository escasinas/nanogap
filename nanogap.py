# # # # # IMPORTS # # # # #

import argparse
import csv
import multiprocessing
import os
import subprocess

# # # # # FUNCTIONS # # # # # 

def make_output_dirs(main_out):
    """
    Make the relevant output directories.

        Params...
    (str) main_out: Name of the parent output directory e.g. ngap_output

        Returns...
    (dict) out_dict: Dictionary with key values as the path name of the directories made
    """
    
    # variables
    assemblies = f"{main_out}/assemblies"
    logs = f"{main_out}/logs"
    flye_logs = f"{logs}/flye"
    medaka_logs = f"{logs}/medaka"
    barrnap_logs = f"{logs}/barrnap"
    rrna_16s_logs = f"{logs}/16S"
    blast_logs = f"{logs}/blast"

    # make directories
    if not os.path.exists(main_out):
        os.mkdir(main_out) # ngap_output
        os.mkdir(assemblies) # ngap_output/assemblies
        os.mkdir(logs) # ngap_output/logs
        os.mkdir(flye_logs) # ngap_output/logs/flye
        os.mkdir(medaka_logs) # ngap_output/logs/medaka
        os.mkdir(barrnap_logs) # ngap_output/logs/barrnap
        os.mkdir(rrna_16s_logs) # ngap_output/logs/16S
        os.mkdir(blast_logs) # ngap_output/logs/blast

    out_dict = {
            "assemblies": assemblies,
            "logs": logs,
            "flye_logs": flye_logs,
            "medaka_logs": medaka_logs,
            "barrnap_logs": barrnap_logs,
            "rrna_16s_logs": rrna_16s_logs,
            "blast_logs": blast_logs
    }

    return out_dict


def run_flye(strain_name, fastq, flye_logs_dir, num_threads):
    """
    Run the Flye assembler.

        Params...
    (str) strain_name: Name of the strain
    (str) fastq: Name of the fastq file path
    (str) logs_dir: Name of the log directory e.g. ngap_output/logs
    (str) assemblies_dir: Name of the assemblies directory e.g. ngap_output/assemblies
    (int) num_threads: Number of threads Flye will use

        Returns...
    (str) draft_fasta: Path to the fasta file of the draft assembly  
    """
    print(f"\n# # # # # FLYE # # # # #\n")

    # variables
    log_strain = f"{flye_logs_dir}/{strain_name}" # ngap_output/logs/flye/<strain>
    
    # run flye on the command line
    flye_cmd = f"flye --nano-raw {fastq} -o {log_strain} --threads {num_threads}"
    subprocess.run(flye_cmd.split())

    # change the name from assembly.fasta to <strain>_draft.fasta 
    assembly_fasta = f"{log_strain}/assembly.fasta"
    draft_fasta = f"{log_strain}/{strain_name}_draft.fasta" # return
    
    # move command to change the name
    mv_cmd = f"mv {assembly_fasta} {draft_fasta}"
    subprocess.run(mv_cmd.split())

    # ngap_output/logs/flye/<strain>_draft.fasta
    return draft_fasta


def parse_flye_log(strain_name, flye_logs_dir):
    """
    Parse the assembly statistics from flye.log to get the relevant assembly information.

        Params...
    (str) strain_name: Name of thr strain
    (str) flye_logs_dir: Path to "ngap_output/logs/flye"

        Returns...
    (dict) assembly_info_dict: Dictionary containing assembly statistics
    """
    # ngap_output/logs/flye/<strain>/flye.out
    flye_log = f"{flye_logs_dir}/{strain_name}/flye.log"
    
    with open(flye_log, 'r') as inf:
        assembly_info = inf.readlines()[-8:-2]

    genome_size = int(assembly_info[0].strip().split()[-1])
    num_contigs = int(assembly_info[1].strip().split()[-1])
    n50 = int(assembly_info[2].strip().split()[-1])
    mean_coverage = int(assembly_info[-1].strip().split()[-1])

    assembly_info_dict = {
            "genome_size": genome_size,
            "num_contigs": num_contigs,
            "n50": n50,
            "mean_coverage": mean_coverage
    }

    return assembly_info_dict


def run_medaka(strain_name, fastq, fasta, model, medaka_logs_dir, assemblies_dir, num_threads):
    """
    Run Medaka for error correction.

        Params...
    (str) strain_name: Name of the strain
    (str) fastq: Path to the fastq file
    (str) fasta: Path to the draft fasta assembled by Flye
    (str) model: Medaka model
    (str) medaka_logs_dir: Path to "ngap_output/logs/medaka"
    (str) assemblies_dir: Path to "ngap_output/assemblies"
    (int) num_threads: Number of threads Medaka will use

        Returns...
    (str) assembly_fasta: Path to the error-corrected consensus fasta file
    """
    print("\n# # # # # MEDAKA # # # # #\n")

    # variables
    log_strain = f"{medaka_logs_dir}/{strain_name}" # ngap_output/logs/medaka/<strain>

    # run medaka_consensus on the command line
    medaka_cmd = f"medaka_consensus -i {fastq} -d {fasta} -o {log_strain} -t {num_threads} -m {model}"
    subprocess.run(medaka_cmd.split())

    # change the name from consensus.fasta to <strain>.fasta
    consensus_fasta = f"{log_strain}/consensus.fasta"
    strain_fasta = f"{log_strain}/{strain_name}.fasta" # ngap_output/logs/medaka/<strain>/<strain>.fasta

    # move command to change the name
    mv_cmd = f"mv {consensus_fasta} {strain_fasta}"
    subprocess.run(mv_cmd.split())

    # copy the consensus assembly to ngap_output/assemblies
    cp_cmd = f"mv {strain_fasta} {assemblies_dir}/"
    subprocess.run(cp_cmd.split())

    # ngap_output/assemblies/<strain>.fasta
    assembly_fasta = f"{assemblies_dir}/{strain_name}.fasta"
    return assembly_fasta


def run_barrnap(strain_name, fasta, barrnap_logs_dir, num_threads):
    """
    Run Barrnap to extract 16S rRNA sequence from the genome assembly.

        Params...
    (str) strain_name: Name of the strain
    (str) fasta: Path to the consensus fasta generated by Medaka
    (str) barrnap_logs_dir: Path to "ngap_output/logs/barrnap"
    (int) num_threads: Number of threads Barrnap will use

        Returns...
    (str) outseq: Path to the output fasta file containing rRNA sequences (16S, 18S, 23S). 

    """
    print("\n# # # # # BARRNAP # # # # #\n")

    # variables
    outseq = f"{barrnap_logs_dir}/{strain_name}_rRNA.fasta" # ngap_output/logs/barrnap/<strain>_rRNA.fasta

    # run barrnap on the command line
    barrnap_cmd = f"barrnap --quiet --kingdom bac --threads {num_threads} --outseq {outseq} {fasta}"
    subprocess.run(barrnap_cmd.split())

    # remove the fasta.fai from the ngap_output/assemblies directory
    fasta_fai = f"{fasta}.fai"
    rm_cmd = f"rm {fasta_fai}"
    subprocess.run(rm_cmd.split())

    # ngap_output/logs/barrnap/<strain>_rRNA.fasta
    return outseq


def extract_16S(strain_name, rrna_fasta, rrna_16s_logs):
    """
    Extract the 16S rRNA sequences from the full (16S, 18S, 23S) rRNA file generated by Barrnap.

        Params...
    (str) strain_name: Name of the strain
    (str) rrna_fasta: Path to the fasta file containing rRNA sequences
    (str) rrna_16s_logs: Path to "ngap_output/logs/16S"

        Returns...
    (tuple) (output_file, count): Path to the output fasta file containing only 16S rRNA sequences.
                                  Number of 16S rRNA sequences.
    """
    output_file = f"{rrna_16s_logs}/{strain_name}_16S.fasta" # ngap_output/logs/16S/<strain>_16S.fasta

    # read the input file (fasta file containing rRNA sequences - 16S, 18S, 23S)
    with open(rrna_fasta, 'r') as inf:
        lines = inf.readlines()
   
    # write to the output file (fasta file containing only 16S rRNA sequences)
    with open(output_file, 'w') as outf:
        count = 1
        for line in lines:
            if "16S" in line:
                # fasta content
                header = f">{strain_name}_16S_rRNA_{count}"
                seq = f"{lines[lines.index(line) + 1]}"

                # write to output file
                outf.write(f"{header}\n{seq}")
                count += 1

        count -= 1

    # ngap_output/logs/16S/<strain>_16S.fasta
    return output_file, count


def make_blast_input(strain_name, fasta_16s, blast_logs_dir):
    """
    Make a fasta file containing one 16S rRNA to be used for BLAST. 

        Params...
    (str) strain_name: Name of the strain.
    (str) fasta_16s: Path to the fasta file containing 16S rRNA sequences.
    (str) blast_logs_dir: Path to "ngap_output/logs/blast".

        Returns...
    (str) fasta_out: Path to the output file containing one 16S rRNA sequence.
    """
    # make a <strain> directory inside the blast logs directory
    blast_strain_dir = f"{blast_logs_dir}/{strain_name}" # ngap_output/logs/blast/<strain>
    os.mkdir(blast_strain_dir)

    # read the input file and extract the first 16S
    with open(fasta_16s, 'r') as inf:
        lines = inf.readlines()[:2]

    # write to output file 
    fasta_out = f"{blast_strain_dir}/{strain_name}_blast.fasta" # ngap_output/logs/blast/<strain>/<strain>_blast.fasta
    fasta_content = ''.join(lines) # file content
    with open(fasta_out, 'w') as outf:
        outf.write(fasta_content)

    # ngap_output/logs/blast/<strain>/<strain>_blast.fasta
    return fasta_out


def run_blast(strain_name, fasta, blast_logs_dir):
    """
    Run BLAST on a single 16S rRNA sequence.

        Params...
    (str) strain_name: Name of the strain.
    (str) fasta: Path to the fasta file containing one 16S rRNA sequence.
    (str) blast_logs_dir: Path to "ngap_output/logs/blast".

        Returns...
    (str) blast_out: Path to the blast output file.
    """
    print("\n# # # # # BLAST # # # # #\n")
    print(f"Performing BLAST search on {strain_name}\n")
    print(f"Input: {fasta}\nDatabase: blastdb/16S_ribosomal_RNA")

    # name the output file
    blast_strain_dir = f"{blast_logs_dir}/{strain_name}" # ngap_output/logs/blast/<strain>
    blast_out = f"{blast_strain_dir}/{strain_name}.out" # ngap_output/logs/blast/<strain>/<strain>.out

    # run blast on the command line
    blast_cmd = f"blastn -db blastdb/16S_ribosomal_RNA -query {fasta} -out {blast_out} -max_target_seqs 100 -sorthits 0"
    subprocess.run(blast_cmd.split())

    print(f"Output: {blast_out}")

    # ngap_output/logs/blast/<strain>/<strain>.out
    return blast_out


def parse_blast_out(blast_out):
    """
    Parse the BLAST output to get the top hit.

        Params...
    (str) blast_out: Path to the BLAST output file.
    
        Returns...
    (dict) top_hit_dict: Dictionary containing BLAST information.
    """
    # read the input file
    with open(blast_out, 'r') as inf:
        lines = inf.readlines()

    # get the top hit
    top_hit = ""
    for line in lines:
        if "Sequences producing significant alignments" in line:
            top_hit = lines[lines.index(line) + 2]

    # parse the top hit into a list and then into a dictionary
    top_hit_ls = top_hit.split()
    top_hit_dict = {
            "accession": top_hit_ls[0],
            "species": ' '.join(top_hit_ls[1:3]),
            "bit_score": top_hit_ls[-5],
            "max_identity": top_hit_ls[-1]
    }

    return top_hit_dict


def process_file_input(fastq_file, output_dir, medaka_model, num_threads):
    """
    Run the pipeline on a single fastq file.

        Params...
    (str) fastq_file: Path to the input fastq file.
    (str) output_dir: Path to the main output directory - default is "ngap_output/".
    (str) medaka_model: Medaka model.
    (str) num_threads: Number of threads to run the pipeline.

        Returns...
    (dict) output_dict: Dictionary containing the combined assembly and BLAST information.
    """
    # Preprocessing
    strain_name = os.path.basename(fastq_file).split('.')[0] 
    read_size = round(os.path.getsize(fastq_file) / 1000000, 1)

    # 1. Make the output directories
    out_dirs = make_output_dirs(output_dir)

    # 2. Run Flye to generate a draft assembly
    draft_fasta = run_flye(strain_name, fastq_file, out_dirs["flye_logs"], num_threads)

    # 3. Run Medaka for error correction and generate a consensus assembly
    consensus_fasta = run_medaka(strain_name, fastq_file, draft_fasta, medaka_model, out_dirs["medaka_logs"], out_dirs["assemblies"], num_threads)

    # 4. Run Barrnap to extract rRNA sequences from the genome assembly
    rrna_fasta = run_barrnap(strain_name, consensus_fasta, out_dirs["barrnap_logs"], num_threads)

    # 5. Extract the 16S rRNA from the Barrnap output
    rrna_16S, num_16S = extract_16S(strain_name, rrna_fasta, out_dirs["rrna_16s_logs"])

    # 6. Make a BLAST input fasta file from the 16S rRNA fasta file
    blast_input = make_blast_input(strain_name, rrna_16S, out_dirs["blast_logs"])

    # 7. Run BLAST - search for 16S sequence homology against the 16S ribosomal RNA BLAST database
    blast_output = run_blast(strain_name, blast_input, out_dirs["blast_logs"])

    # 8. Parse the BLAST output to get the relevant information
    parsed_blast_out = parse_blast_out(blast_output)

    # 9. Parse the Flye log to get some assembly statistics
    parsed_flye_log = parse_flye_log(strain_name, out_dirs["flye_logs"])

    # 10. Combine all of the relevant information and assembly statistics into one dictionary
    output_dict = {
            "strain": strain_name,
            "top_hit": parsed_blast_out["species"],
            "accession": parsed_blast_out["accession"],
            "bit_score": parsed_blast_out["bit_score"],
            "max_identity": parsed_blast_out["max_identity"],
            "read_size": read_size,
            "genome_size": parsed_flye_log["genome_size"],
            "contigs": parsed_flye_log["num_contigs"],
            "n50": parsed_flye_log["n50"],
            "mean_coverage": parsed_flye_log["mean_coverage"],
            "num_16S": num_16S
    }

    print(f"\n# # # # # FINISHED ASSEMBLING {strain_name} # # # # #\n")

    return output_dict


def write_one_row(output_dir, output_dict):
    """
    Generate a csv file and write one row containing assembly information.

        Params...
    (str) output_dir: Path to the main output directory - default is "ngap_output/".
    (dict) output_dict: Dictionary containing the combined assembly and BLAST information.

        Returns...
    (void)
    """
    file_name = f"{output_dir}/{output_dir}.csv" # ngap_output/ngap_output.csv
    with open(file_name, 'w') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=output_dict.keys())

        # write the csv columns
        writer.writeheader()
       
        # write the single row
        writer.writerow(output_dict)


def process_directory_input(fastq_dir, output_dir, medaka_model, num_threads):
    """
    Run the pipeline on a directory containing fastq files.

        Params...
    (str) fastq_dir: Path to the input directory containing fastq files.
    (str) output_dir: Path to the main output directory - default is "ngap_output/".
    (str) medaka_model: Medaka model.
    (str) num_threads: Number of threads to run the pipeline.

        Returns...
    (list) output_dict: List containing dictionaries of the combined assembly and BLAST information per strain.
    """
    # user will be able to input both "fastq_reads/" and "fastq_reads"
    fastq_dir = fastq_dir.strip('/')

    # list of dictionaries containing assembly information
    multi_output_info = []
    
    # list of fastq files
    fastq_list = os.listdir(fastq_dir)
    for fastq in fastq_list:
        input_read = f"{fastq_dir}/{fastq}"
        output_info = process_file_input(input_read, output_dir, medaka_model, num_threads)

        # append the assembly info
        multi_output_info.append(output_info)

    return multi_output_info


def write_multi_rows(output_dir, output_list):
    """
    Generate a csv file and write multiple rows containing assembly information for each strain.

        Params...
    (str) output_dir: Path to the main output directory - default is "ngap_output/".
    (list) output_list: List containing dictionaries of the combined assembly and BLAST information per strain.

        Returns...
    (void)
    """
    file_name = f"{output_dir}/{output_dir}.csv" # ngap_output/ngap_output.csv
    with open(file_name, 'w') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=output_list[0].keys())

        # write the csv columns
        writer.writeheader()

        # write the rows
        for output_dict in output_list:
            writer.writerow(output_dict)

# # # # # PARSER # # # # # 

MAX_CPU = multiprocessing.cpu_count()

parser = argparse.ArgumentParser(description='Assembly and correction of long ONT .fastq reads')
parser.add_argument('input', type=str, help='Input single fastq read file or directory containing fastq reads')
parser.add_argument('-t', '--threads', default=MAX_CPU, help='Number of threads to run NanoGAP (default: MAX_CPU)')
parser.add_argument('-o', '--outdir', default='ngap_output', help='Name of output directory (default: ngap_output)')
parser.add_argument('-m', '--model', default='r941_min_high_g360', help='Medaka model (default: r941_min_high_g360)')
args = parser.parse_args()

input_read = args.input # directory containing reads
num_threads = args.threads # number of cores/threads to be used by NanoGAP
output_dir = args.outdir # name of output dir
medaka_model = args.model # medaka model option

# # # # # PROGRAM # # # # # 
if __name__ == "__main__":
    
    if os.path.isfile(input_read):
        # run the pipeline on a single fastq file
        output_info = process_file_input(input_read, output_dir, medaka_model, num_threads)

        # write the output csv containing assembly information
        write_one_row(output_dir, output_info)

    elif os.path.isdir(input_read):
        # run the pipeline on a directory containing more than one fastq files
        output_list = process_directory_input(input_read, output_dir, medaka_model, num_threads)

        # write the output csv containing multiple rows of assembly information
        write_multi_rows(output_dir, output_list)

    else:
        print("Invalid input file/directory")
