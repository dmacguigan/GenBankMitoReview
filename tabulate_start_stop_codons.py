#!/usr/bin/env python3
"""
Script to read FASTA files containing mitochondrial genes and count
the first 3 bp (start codon) and last 3 bp (stop codon) of each sequence.
"""

import os
import sys
import glob
from collections import Counter

def print_help():
    """
    Print usage information and exit.
    """
    help_text = """
Usage: tabulate_start_stop_codons.py [OPTIONS]

Description:
  Reads all FASTA files in the current directory and tabulates the start codon
  (first 3 bp) and stop codon (last 3 bp, or remainder if length % 3 != 0)
  for each sequence. Sequences composed entirely of 'N' or '?' are skipped.

  Supported file extensions: .fasta, .fa, .fna

Options:
  -v, --verbose     Print per-sequence details (ID, length, start/stop codon)
  -h, --help        Show this help message and exit

Output:
  For each FASTA file, prints:
    - Total number of sequences (and any skipped)
    - Sequence length statistics (mean, min, max)
    - Start codon counts and percentages
    - Stop codon counts and percentages
    - [verbose] Per-sequence breakdown

Examples:
  # Run in a directory containing .fasta files
  python tabulate_start_stop_codons.py

  # Include per-sequence details
  python tabulate_start_stop_codons.py --verbose

Notes:
  - Sequences are converted to uppercase before processing.
  - The stop codon is the last 3 bp for sequences whose length is divisible
    by 3, or the last (length % 3) bp otherwise.
  - Only the first word of each FASTA header (after '>') is used as the ID.
"""
    print(help_text.strip())
    sys.exit(0)

def read_fasta(filepath):
    """
    Read a FASTA file and return a dictionary of sequences.
    
    Args:
        filepath: Path to the FASTA file
        
    Returns:
        Dictionary with sequence IDs as keys and sequences as values
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                # Start new sequence
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        
        # Save last sequence
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences

def is_invalid_sequence(sequence):
    """
    Check if a sequence is entirely composed of '?' or 'N'.
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        True if sequence is invalid, False otherwise
    """
    return all(base in ['?', 'N'] for base in sequence)

def extract_codons(sequence):
    """
    Extract the first 3 bp (start codon) and last 3 bp (stop codon) from a sequence.
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        Tuple of (start_codon, stop_codon)
    """
    seq_len = len(sequence)
    
    # Extract first 3 bp (or fewer if sequence is shorter)
    start_codon = sequence[:min(3, seq_len)]
    
    # Extract last 3 bp (or remaining bp if not divisible by 3)
    remainder = seq_len % 3
    if remainder == 0:
        stop_codon = sequence[-3:] if seq_len >= 3 else sequence
    else:
        stop_codon = sequence[-remainder:]
    
    return start_codon, stop_codon

def analyze_fasta_files(fasta_files):
    """
    Analyze start and stop codons for all sequences in the provided FASTA files.
    
    Args:
        fasta_files: List of FASTA file paths
        
    Returns:
        Dictionary with analysis results
    """
    results = {}
    
    for filepath in fasta_files:
        gene_name = os.path.basename(filepath)
        sequences = read_fasta(filepath)
        
        start_counter = Counter()
        stop_counter = Counter()
        sequence_details = []
        sequence_lengths = []
        skipped_sequences = []
        
        for seq_id, sequence in sequences.items():
            # Check if sequence is invalid
            if is_invalid_sequence(sequence):
                skipped_sequences.append(seq_id)
                print(f"WARNING: Skipping sequence '{seq_id}' in {gene_name} - "
                      f"entirely composed of '?' or 'N'")
                continue
            
            start_codon, stop_codon = extract_codons(sequence)
            
            start_counter[start_codon] += 1
            stop_counter[stop_codon] += 1
            sequence_lengths.append(len(sequence))
            
            sequence_details.append({
                'id': seq_id,
                'length': len(sequence),
                'start_codon': start_codon,
                'stop_codon': stop_codon
            })
        
        # Calculate length statistics
        length_stats = {}
        if sequence_lengths:
            length_stats = {
                'mean': sum(sequence_lengths) / len(sequence_lengths),
                'min': min(sequence_lengths),
                'max': max(sequence_lengths)
            }
        
        results[gene_name] = {
            'start_codons': start_counter,
            'stop_codons': stop_counter,
            'sequences': sequence_details,
            'total_sequences': len(sequence_details),
            'skipped_sequences': skipped_sequences,
            'length_stats': length_stats
        }
    
    return results

def print_summary(results, verbose=False):
    """
    Print a formatted summary of the codon analysis.
    
    Args:
        results: Dictionary with analysis results
        verbose: If True, print detailed information for each sequence
    """
    print("=" * 80)
    print("MITOCHONDRIAL GENE START/STOP CODON SUMMARY")
    print("=" * 80)
    
    # Per-file detailed results
    for gene_name, data in results.items():
        print(f"\n{'=' * 80}")
        print(f"File: {gene_name}")
        print(f"Number of sequences: {data['total_sequences']}")
        if data['skipped_sequences']:
            print(f"Skipped sequences: {len(data['skipped_sequences'])}")
        print(f"{'=' * 80}")
        
        # Length statistics
        if data['length_stats']:
            stats = data['length_stats']
            print(f"\n  Sequence Length Statistics:")
            print(f"    Mean: {stats['mean']:.1f} bp")
            print(f"    Min:  {stats['min']} bp")
            print(f"    Max:  {stats['max']} bp")
        
        print("\n  Start Codons (first 3 bp):")
        for codon, count in sorted(data['start_codons'].items(), key=lambda x: (-x[1], x[0])):
            percentage = (count / data['total_sequences']) * 100
            print(f"    {codon}: {count:4d} ({percentage:5.1f}%)")
        
        print("\n  Stop Codons (last 3 bp or remainder):")
        for codon, count in sorted(data['stop_codons'].items(), key=lambda x: (-x[1], x[0])):
            percentage = (count / data['total_sequences']) * 100
            print(f"    {codon}: {count:4d} ({percentage:5.1f}%)")
        
        # Only print sequence details if verbose mode is on
        if verbose:
            print("\n  Sequence Details:")
            for seq in data['sequences']:
                print(f"    {seq['id']:30s} Length: {seq['length']:6d}  Start: {seq['start_codon']:3s}  Stop: {seq['stop_codon']:3s}")

def main():
    """
    Main function to run the analysis.
    """
    # Check for help flag first
    if '-h' in sys.argv or '--help' in sys.argv:
        print_help()

    # Check for verbose flag
    verbose = '-v' in sys.argv or '--verbose' in sys.argv

    # Warn about any unrecognized flags
    known_flags = {'-v', '--verbose', '-h', '--help'}
    for arg in sys.argv[1:]:
        if arg.startswith('-') and arg not in known_flags:
            print(f"WARNING: Unrecognized option '{arg}' (ignored). Use -h for help.")

    # Get all FASTA files in the current directory
    fasta_files = glob.glob("*.fasta") + glob.glob("*.fa") + glob.glob("*.fna")
    
    if not fasta_files:
        print("No FASTA files found in the current directory.")
        print("Please ensure your files have .fasta, .fa, or .fna extensions.")
        print("Use -h or --help for usage information.")
        return
    
    print(f"Found {len(fasta_files)} FASTA file(s):")
    for f in fasta_files:
        print(f"  - {f}")
    print()
    
    # Analyze the files
    results = analyze_fasta_files(fasta_files)
    
    # Print summary
    print_summary(results, verbose=verbose)
    
    if not verbose:
        print("\n" + "=" * 80)
        print("Use -v or --verbose flag to see detailed sequence information")
        print("Use -h or --help for full usage information")
        print("=" * 80)

if __name__ == "__main__":
    main()
