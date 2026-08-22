import argparse
import pandas as pd
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--blast', required=True, help='BLAST results file (CSV)')
    parser.add_argument('--seqs', required=True, help='Input sequences file (FASTA)')
    parser.add_argument('--out_seqs', required=True, help='Output sequences file (FASTA)')
    parser.add_argument('--range', required=True, help='Specify the region to filter: protein_xy')
    parser.add_argument('--min_length', type=int, required=True, help='Minimum length for the specified range')
    parser.add_argument('--max_length', type=int, required=True, help='Maximum length for the specified range')
    parser.add_argument('--no-debug', action='store_false', dest='debug', help='Disable debug output (sequence & filter counts)')
    return parser.parse_args()

def main():
    args = parse_args()
    length_range = (args.min_length, args.max_length)  # Use the provided min and max lengths
    
    # Load BLAST results
    blast_results = pd.read_csv(args.blast, 
                                names=["qseqid","sseqid","pident","length","mismatch","gapopen",
                                        "qstart","qend","sstart","send","evalue","bitscore","qcovs"])
    if args.debug: print("\n*** Loading counts ***")
    if args.debug: print(f"Loaded {len(blast_results)} rows from BLAST file")

    # Load sequences
    sequences = list(SeqIO.parse(args.seqs, "fasta"))
    if args.debug: print(f"Loaded {len(sequences)} sequences from FASTA file")

    # Filter blast results and remove duplicates
    blast_results=blast_results.sort_values(['length','evalue'],ascending=[False,True]).drop_duplicates(subset='qseqid', keep='first')
    blast_results["diff_length_ref"]=abs(blast_results["send"]-blast_results["sstart"])
    blast_results.loc[:,["qseqid","diff_length_ref"]].to_csv(f"{args.range}/results/blast_{args.range}_length.csv",index=False,header=False) # needed later for adding to metadata
    selected_seqs = []
    if args.debug: print(f"After sort/deduplicate: {len(blast_results)} sequences")

    if args.debug: print(f"Filtering with length_range={length_range}\n")
    if args.debug: print("*** Filtering counts ***")

    # Process sequences based on the specified length range
    if args.range != "whole_genome":
            keep_mask = (blast_results.diff_length_ref >= length_range[0]) & (blast_results.diff_length_ref <= length_range[1])
            if args.debug: print(f"Sequences passing length filter: {keep_mask.sum()}")
            if args.debug: print(f"Sequences rejected by length filter: {(~keep_mask).sum()}")

            if (~keep_mask).sum() > 0:
                rejected = blast_results[~keep_mask]
                if args.debug: print(f"Rejected sequences' diff_length_ref range: {rejected['diff_length_ref'].min()} to {rejected['diff_length_ref'].max()}")

            blast_results = blast_results[keep_mask]

            if args.debug: print(f"After length filter: {len(blast_results)} sequences pass")

            if args.debug: print("\n*** Beginning secondary length check ***")
            for seq_record in sequences:
                if seq_record.id in blast_results.qseqid.unique():
                    # Get the relevant BLAST hit
                    blast_hit = blast_results[blast_results.qseqid == seq_record.id].iloc[0]
                    # Extract the aligned region
                    reg_seq = seq_record.seq[blast_hit["qstart"] - 1:blast_hit["qend"]]
                    # Create new sequence record
                    new_seq_record = seq_record[:0]  # Copy metadata
                    new_seq_record.seq = reg_seq
                    if length_range[0] <= len(reg_seq) <= length_range[1]:
                        selected_seqs.append(new_seq_record)
                    else:
                        if args.debug: print(f"\t{seq_record.id} REJECTED on second check - extracted length {len(reg_seq)}, outside range {length_range}")
                        if args.debug: print(f"\t{seq_record.id} - ref_length={blast_hit['diff_length_ref']}, query_extracted_length={len(reg_seq)}")

    else:
        for seq_record in sequences:
            seq_length = len(seq_record.seq)
            if length_range[0] <= seq_length <= length_range[1]:
                selected_seqs.append(seq_record)
        if args.debug: print(f"After whole genome length filter: {len(selected_seqs)} sequences pass\n")

    # Write selected sequences
    SeqIO.write(selected_seqs, args.out_seqs, "fasta")
    print(f"Filtering retained {len(selected_seqs)} out of {len(sequences)} sequences in this file")

if __name__ == "__main__":
    main()