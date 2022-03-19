import os
import pandas as pd
import sys

def main():
    """
    Accepts two arguments:
    1. a filename of a CSV file describing the selected constructs
    for which subsequent forward-folding simulations should be conducted.
    the only columns of that CSV that matter are the name (which will serve
    as a unique identifier for the forward folding simulation/designed variant)
    and `seq1` and `seq2`.
    2. a fasta file with two 'wildcard' regions for substitution. to permit
    modification using simple string replacement, these wildcard regions should
    be labeled [!SEQ1] and [!SEQ2], strings that should never appear in any 
    Rosetta-compatible fasta file.

    Run this script in the directory where you ran the design simulation, so
    the starting PDB is available. (Or a local working directory where you
    set up the design simulation in the first place.) This script will copy the
    starting PDB file into each of the desired new forward folding directories.
    It will also write the contents of your current flags file, although it will
    reduce the number of cycles to an appropriate number for forward folding.
    """

    # you can run this in your target working directory on the cluster, or you
    # can rsync this directory to the cluster. either way, it's pretty light.
    try:
        os.mkdir("indiv_simulations")
    except FileExistsError:
        pass
    
    if len(sys.argv) < 4:
        print('Usage: python parse_out_sequences.py [csv_file] [fasta_template] [starting_pdb]')
        quit()

    csv_file = sys.argv[1]
    fasta_template_file = sys.argv[2]
    starting_pdb = sys.argv[3]

    # Why parse the command line when you can just warn the user about mildly
    # suspicious decisions?
    if 'fasta' in csv_file:
        print("Warning: 'fasta' found in your CSV filename")
    if 'csv' in fasta_template_file:
        print("Warning: 'csv' found in your fasta filename")

    df = pd.read_csv(csv_file)
    fasta_contents, flags_lines = None, None
    with open(fasta_template_file) as f:
        fasta_contents = f.read()
    with open('flags') as f:
        flags_lines = f.readlines()

    for ii, row in df.iterrows():
        # os.chdir("indiv_simulations")
        seq1 = row.seq1
        seq2 = row.seq2
        name = row.name
        # name, sequence, seq1, seq2, score = line.strip().split(',')[:5]

        row_fasta_contents = fasta_contents.replace('[!SEQ1]', seq1.lower())
        row_fasta_contents = row_fasta_contents.replace('[!SEQ2]', seq2.lower())

        try:
            os.mkdir(f'indiv_simulations/{name}')
        except FileExistsError:
            pass

        with open(f"indiv_simulations/{name}/target.fasta", "w") as g:
            g.write()

        # Copies the current working directory's 
        os.system(f"cp {starting_pdb}.pdb ./indiv_simulations/{name}/")
        

        with open(f"indiv_simulations/{name}/flags", "w") as g:
            for line in flags_lines:
                if 'cycles' in line:
                    g.write('-cycles 500\n')
                else:
                    g.write(line)
        with open(f"indiv_simulations/{name}/README_SWM", "w") as g:
            g.write("stepwise @flags\n")

if __name__ == '__main__':
    main()
