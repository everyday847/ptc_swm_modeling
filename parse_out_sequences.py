import re
import json
import sys

def seq_map(target_len, grouped_master_seq_idxs, annotated):
    # remove all [] parts from annotated.
    annotated = annotated.replace(':NtermProteinFull', '')
    annotated = annotated.replace(':CtermProteinFull', '')
    annotated = annotated.replace(':Virtual_Protein_SideChain', '')
    annotated = annotated.replace(':Virtual_Phosphate', '')
    annotated = annotated.replace(':rna_cutpoint_lower', '')
    annotated = annotated.replace(':rna_cutpoint_upper', '')
    annotated = annotated.replace(':UpperRNA', '')
    annotated = annotated.replace('[HIS_D]', '')
    annotated = re.sub(r'\[...\]', '', annotated)

    if len(annotated) == target_len:
        tot = [''.join([annotated[i] for i in group]) for group in grouped_master_seq_idxs]
        # We want this to serve as a dict key, and you can't hash a list
        # but you can't do a 'tuple comprehension' above of course.
        return tuple(tot)
    else:
        return ''

def main():
    if len(sys.argv) < 2:
        print('Usage: python parse_out_sequences.py [silent_file]')
        quit()

    fn = sys.argv[1]
    seq_scores = {}
    missing_idx = None
    with open(fn) as f:
        score_current = None
        seq_current = None
        grouped_master_seq_idxs = []
        target_len = None
        missing = None
        for line in f:
            if line[:len('SEQUENCE:')] == 'SEQUENCE:':
                if target_len == None: 
                    master_seq = line.split()[1].strip()
                    target_len = len(master_seq)
                    for i, c in enumerate(master_seq):
                        if c in 'nbdhvwskmyr':
                            # Unnecessarily explicit logic, but: accumulate consecutive
                            # groups of string indices, so that we can rapidly parse out
                            # consecutive sequence groups from completed ANNOTATED_SEQUENCE
                            # lines of individual completed models.

                            # We only care about these 11 letters because acgu are fixed by
                            # definition and aren't part of the 'problem definition'
                            if len(grouped_master_seq_idxs) == 0:
                                grouped_master_seq_idxs.append([i])
                            elif i == max(grouped_master_seq_idxs[-1]) + 1:
                                grouped_master_seq_idxs[-1].append(i)
                            else:
                                grouped_master_seq_idxs.append([i])

            elif line[:len('SCORE:')] == 'SCORE:':
                if 'description' not in line:
                    # store score
                    score_current = float(line.split()[1])
                    missing = int(line.split()[missing_idx])
                else:
                    missing_idx = line.split().index('missing')
            elif line[:len('ANNOTATED_SEQUENCE:')] == 'ANNOTATED_SEQUENCE:':
                # store specific sequence
                seq_current = seq_map(target_len, grouped_master_seq_idxs, line.split()[1].strip())
            elif 'REMARK' in line:
                # use this as our trigger to reset
                if seq_current is None: continue
                if missing == 0:
                    if seq_current in seq_scores:
                        seq_scores[seq_current].append(score_current)
                    else:
                        seq_scores[seq_current] = [score_current]

                seq_current = None
                score_current = None

    # record the whole dict for posterity.
    with open(fn.replace('.out', '.json'), 'w') as g:
        g.write(json.dumps(seq_scores))

    # Here's an example of how you might write out the top 100 sequences
    # to a csv file suitable for use with set_up_forward_folding.py. We
    # sort based on the minimum sampled score

    # you can't sort a dict, but you can sort its items -- in this case, the 
    # sort-key is set up to be min(item[1]), i.e., the minimum value of the
    # score-list we've stored in this dictionary anyway.
    with open('top_100_seqs.csv', 'w') as f:
        f.write('name,seq1,seq2\n')
        for ii, (k, v) in enumerate(sorted(seq_scores.items(), key=lambda item: min(item[1]))[:100]):
            seq1 = k[0]
            seq2 = k[1]
            f.write(f'seq_{ii},{seq1},{seq2}\n')

if __name__ == '__main__':
    main()