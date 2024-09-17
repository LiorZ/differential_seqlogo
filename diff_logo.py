# diff_logo_with_substitution_matrix.py

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.transforms import Affine2D
from scipy.stats import gamma
from statsmodels.stats.multitest import multipletests
from Bio.SubsMat import MatrixInfo
import tqdm

# Define Alphabets

class Alphabet:
    def __init__(self, chars, cols, support_reverse_complement):
        self.chars = chars
        self.cols = cols
        self.size = len(chars)
        self.support_reverse_complement = support_reverse_complement
        self.char_to_index = {char: idx for idx, char in enumerate(chars)}

# DNA Alphabet
DNA = Alphabet(
    chars=['A', 'C', 'G', 'T'],
    cols=['green', 'blue', 'orange', 'red'],
    support_reverse_complement=True
)

# RNA Alphabet
RNA = Alphabet(
    chars=['A', 'C', 'G', 'U'],
    cols=['green', 'blue', 'orange', 'red'],
    support_reverse_complement=True
)

# Protein Alphabet (Standard 20 Amino Acids)
PROTEIN = Alphabet(
    chars=[
        'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G',
        'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
        'T', 'W', 'Y', 'V'
    ],
    cols=[
        'darkgreen', 'purple', 'red', 'red', 'green', 'red', 'red', 'darkgreen',
        'magenta', 'orange', 'orange', 'purple', 'orange', 'blue', 'pink', 'brown',
        'brown', 'blue', 'blue', 'orange'
    ],
    support_reverse_complement=False  # Reverse complement is not applicable for proteins
)

# Available substitution matrices in Biopython
available_matrices = {
    'BLOSUM45': MatrixInfo.blosum45,
    'BLOSUM50': MatrixInfo.blosum50,
    'BLOSUM62': MatrixInfo.blosum62,
    'BLOSUM80': MatrixInfo.blosum80,
    'BLOSUM90': MatrixInfo.blosum90,
}

# Functions for handling sequences and PWMs

def get_sequences_from_fasta_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    sequences = []
    seq = ''
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if seq != '':
                sequences.append(seq)
                seq = ''
        else:
            seq += line.upper()
    if seq != '':
        sequences.append(seq)
    return sequences

def get_counts_from_sequences(sequences, alphabet):
    alignment_length = len(sequences[0])
    counts = np.zeros((alignment_length, alphabet.size), dtype=int)
    for seq in sequences:
        for pos, char in enumerate(seq):
            if char in alphabet.char_to_index:
                idx = alphabet.char_to_index[char]
                counts[pos, idx] += 1
            else:
                pass  # Ignore unknown characters or handle as needed
    return counts

def counts_to_pwm(counts, pseudo_count=0.1):
    pwm = counts + pseudo_count
    pwm = pwm / pwm.sum(axis=1, keepdims=True)
    return pwm.T  # Transpose to match expected dimensions

def information_content(p):
    p = np.array(p)
    p_nonzero = p > 0
    ic = np.log2(len(p)) + np.sum(p[p_nonzero] * np.log2(p[p_nonzero]))
    return {'height': ic, 'ylab': 'Information Content [bits]'}

def shannon_divergence(p1, p2):
    p1 = np.array(p1)
    p2 = np.array(p2)
    if np.allclose(p1, p2):
        height = 0
    else:
        m = (p1 + p2) / 2
        p1_nonzero = p1 > 0
        p2_nonzero = p2 > 0
        term1 = 0.5 * np.sum(p1[p1_nonzero] * (np.log2(p1[p1_nonzero]) - np.log2(m[p1_nonzero])))
        term2 = 0.5 * np.sum(p2[p2_nonzero] * (np.log2(p2[p2_nonzero]) - np.log2(m[p2_nonzero])))
        height = term1 + term2
    return {'height': height, 'ylab': 'JS divergence'}

def create_substitution_matrix(alphabet, matrix_name):
    if matrix_name not in available_matrices:
        raise ValueError(f"Substitution matrix '{matrix_name}' is not available.")
    matrix_dict = available_matrices[matrix_name]
    size = len(alphabet.chars)
    sub_matrix = np.zeros((size, size))
    for i, aa1 in enumerate(alphabet.chars):
        for j, aa2 in enumerate(alphabet.chars):
            if (aa1, aa2) in matrix_dict:
                score = matrix_dict[(aa1, aa2)]
            elif (aa2, aa1) in matrix_dict:
                score = matrix_dict[(aa2, aa1)]
            else:
                score = 0  # Assign a default score if not found
            sub_matrix[i, j] = score
    # Normalize the matrix to obtain probabilities
    min_score = sub_matrix.min()
    adjusted_matrix = sub_matrix - min_score  # Make all scores positive
    sub_probs = adjusted_matrix / adjusted_matrix.sum(axis=1, keepdims=True)
    return sub_probs

def compute_position_divergences(pwm1, pwm2, alphabet, sub_probs=None):
    npos = pwm1.shape[1]
    divergences = np.zeros(npos)
    for j in range(npos):
        p1 = pwm1[:, j]
        p2 = pwm2[:, j]
        if alphabet == PROTEIN and sub_probs is not None:
            # Use substitution matrix probabilities
            divergence = substitution_divergence(p1, p2, sub_probs)
        else:
            d = shannon_divergence(p1, p2)
            divergence = d['height']
        divergences[j] = divergence
    return divergences

def substitution_divergence(p1, p2, sub_probs):
    p1 = np.array(p1)
    p2 = np.array(p2)
    # Compute expected probabilities under substitution matrix
    expected_p1 = sub_probs.dot(p1)
    expected_p2 = sub_probs.dot(p2)
    # Use Kullback-Leibler divergence as a measure
    kl_div1 = np.sum(p1 * np.log2(p1 / expected_p1 + 1e-10))
    kl_div2 = np.sum(p2 * np.log2(p2 / expected_p2 + 1e-10))
    divergence = (kl_div1 + kl_div2) / 2
    return divergence

def permutation_test(sequences1, sequences2, alphabet, n_permutations=1000, sub_probs=None):
    combined_sequences = sequences1 + sequences2
    n_samples1 = len(sequences1)
    n_total = len(combined_sequences)
    n_positions = len(sequences1[0])

    # Observed divergence
    counts1 = get_counts_from_sequences(sequences1, alphabet)
    counts2 = get_counts_from_sequences(sequences2, alphabet)
    pwm1 = counts_to_pwm(counts1)
    pwm2 = counts_to_pwm(counts2)
    observed_divergences = compute_position_divergences(pwm1, pwm2, alphabet, sub_probs=sub_probs)

    divergence_null = np.zeros((n_permutations, n_positions))
    for i in tqdm.tqdm(range(n_permutations)):
        permuted_indices = np.random.permutation(n_total)
        perm_sequences1 = [combined_sequences[idx] for idx in permuted_indices[:n_samples1]]
        perm_sequences2 = [combined_sequences[idx] for idx in permuted_indices[n_samples1:]]
        perm_counts1 = get_counts_from_sequences(perm_sequences1, alphabet)
        perm_counts2 = get_counts_from_sequences(perm_sequences2, alphabet)
        perm_pwm1 = counts_to_pwm(perm_counts1)
        perm_pwm2 = counts_to_pwm(perm_counts2)
        divergence_null[i, :] = compute_position_divergences(perm_pwm1, perm_pwm2, alphabet, sub_probs=sub_probs)
    p_values = np.mean(divergence_null >= observed_divergences, axis=0)
    return p_values

def add_letter(ax, letter, x_pos, y_pos, ht, wt, color):
    fp = FontProperties(family='DejaVu Sans', weight='bold')
    tp = TextPath((0, 0), letter, size=1, prop=fp)
    t = Affine2D().scale(wt, ht) + Affine2D().translate(x_pos, y_pos) + ax.transData
    patch = PathPatch(tp, transform=t, color=color, linewidth=0)
    ax.add_patch(patch)

def seq_logo(pwm, ax, alphabet=DNA, stack_height=information_content, base_distribution=lambda p: p):
    pwm = np.array(pwm)
    npos = pwm.shape[1]
    max_height = 0
    for j in range(npos):
        column = pwm[:, j]
        sh = stack_height(column)
        heights = base_distribution(column) * sh['height']
        letter_order = np.argsort(heights)
        y_pos = 0
        for i in letter_order[::-1]:
            ht = heights[i]
            if ht > 0:
                letter = alphabet.chars[i]
                color = alphabet.cols[i]
                add_letter(ax, letter, x_pos=j, y_pos=y_pos, ht=ht, wt=1, color=color)
                y_pos += ht
        if y_pos > max_height:
            max_height = y_pos
    ax.set_xlim(-0.5, npos - 0.5)
    ax.set_ylim(0, max_height + 0.5)
    ax.set_xticks(range(npos))
    ax.set_xticklabels(range(1, npos + 1))
    ax.set_ylabel(sh['ylab'])

def diff_logo(pwm1, pwm2, ax, alphabet=DNA, p_values=None, alpha=0.05, sub_probs=None):
    pwm1 = np.array(pwm1)
    pwm2 = np.array(pwm2)
    npos = pwm1.shape[1]
    max_height = 0
    min_height = 0
    for j in range(npos):
        p1 = pwm1[:, j]
        p2 = pwm2[:, j]
        if alphabet == PROTEIN and sub_probs is not None:
            sh = {'height': substitution_divergence(p1, p2, sub_probs)}
            diff = p1 - p2
            if np.sum(np.abs(diff)) != 0:
                heights = diff / np.sum(np.abs(diff)) * sh['height']
            else:
                heights = np.zeros_like(diff)
        else:
            sh = shannon_divergence(p1, p2)
            diff = p1 - p2
            if np.sum(np.abs(diff)) != 0:
                heights = diff / np.sum(np.abs(diff)) * sh['height']
            else:
                heights = np.zeros_like(diff)
        letter_order = np.argsort(np.abs(heights))
        y_pos_pos = 0
        y_pos_neg = 0
        for i in letter_order[::-1]:
            ht = heights[i]
            if ht > 0:
                y_pos = y_pos_pos
                y_pos_pos += ht
            else:
                y_pos = y_pos_neg
                y_pos_neg += ht
            letter = alphabet.chars[i]
            color = alphabet.cols[i]
            add_letter(ax, letter, x_pos=j, y_pos=y_pos, ht=ht, wt=1, color=color)
        max_height = max(max_height, y_pos_pos)
        min_height = min(min_height, y_pos_neg)
    ax.set_xlim(-0.5, npos - 0.5)
    ax.set_ylim(min_height - 0.5, max_height + 0.5)
    ax.axhline(0, color='black', linewidth=0.5)
    ax.set_xticks(range(npos))
    ax.set_xticklabels(range(1, npos + 1))
    ax.set_ylabel('Difference')
    ax.set_title('Difference Logo')
    # Add significance indicators
    if p_values is not None:
        for j in range(npos):
            if p_values[j] < alpha:
                ax.text(j, min_height - 0.1, '*', ha='center', va='bottom', fontsize=14, color='black')

def main():
    parser = argparse.ArgumentParser(description='Generate sequence logos and difference logo from two FASTA files.')
    parser.add_argument('fasta1', type=str, help='Path to the first FASTA file.')
    parser.add_argument('fasta2', type=str, help='Path to the second FASTA file.')
    parser.add_argument('--alphabet', type=str, default='DNA', choices=['DNA', 'RNA', 'PROTEIN'],
                        help='Alphabet to use: DNA, RNA, or PROTEIN. Default is DNA.')
    parser.add_argument('--matrix', type=str, default='BLOSUM62',
                        help='Substitution matrix to use for protein sequences. Default is BLOSUM62.')
    parser.add_argument('--n_permutations', type=int, default=1000,
                        help='Number of permutations for statistical significance testing. Default is 1000.')
    parser.add_argument('--pval-csv', type=str, default=None,
                        help='Path to save LOG(p-values) as a CSV file. Default is None. if None, p-values are not saved.')
    parser.add_argument('--fig-out', type=str, default=None,
                        help='Path to save the figure. Default is None. if None, the figure is displayed.')
    args = parser.parse_args()

    # Select the alphabet
    if args.alphabet.upper() == 'DNA':
        alphabet = DNA
        sub_probs = None
    elif args.alphabet.upper() == 'RNA':
        alphabet = RNA
        sub_probs = None
    elif args.alphabet.upper() == 'PROTEIN':
        alphabet = PROTEIN
        # Load the specified substitution matrix
        try:
            sub_probs = create_substitution_matrix(alphabet, args.matrix.upper())
            print(f"Using substitution matrix: {args.matrix.upper()}")
        except ValueError as e:
            print(e)
            sys.exit(1)
    else:
        print(f"Alphabet {args.alphabet} is not recognized. Using DNA alphabet.")
        alphabet = DNA
        sub_probs = None

    # Read sequences from FASTA files
    sequences1 = get_sequences_from_fasta_file(args.fasta1)
    sequences2 = get_sequences_from_fasta_file(args.fasta2)

    # Check if sequences are of equal length
    lengths1 = set(len(seq) for seq in sequences1)
    lengths2 = set(len(seq) for seq in sequences2)
    if len(lengths1) > 1 or len(lengths2) > 1:
        print("Sequences in each FASTA file must be of the same length.")
        sys.exit(1)
    if lengths1 != lengths2:
        print("Sequences in both FASTA files must be of the same length.")
        sys.exit(1)

    # Generate counts and PWMs
    counts1 = get_counts_from_sequences(sequences1, alphabet)
    counts2 = get_counts_from_sequences(sequences2, alphabet)
    pwm1 = counts_to_pwm(counts1)
    pwm2 = counts_to_pwm(counts2)

    # Compute divergence
    divergence = compute_position_divergences(pwm1, pwm2, alphabet, sub_probs=sub_probs).sum()
    print(f"Total divergence: {divergence}")

    # Compute p-values for each position using permutation test
    print("Computing p-values using permutation test...")
    p_values = permutation_test(sequences1, sequences2, alphabet, n_permutations=args.n_permutations, sub_probs=sub_probs)
    # Adjust p-values for multiple testing using Benjamini-Hochberg
    reject, pvals_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

    # Plot sequence logos and difference logo
    npos = pwm1.shape[1]
    fig_height = 6
    fig, axs = plt.subplots(3, 1, figsize=(npos, fig_height * 3))

    # Sequence logo for PWM1
    axs[0].set_title('Sequence Logo for PWM1')
    seq_logo(pwm1, axs[0], alphabet)
    axs[0].set_xlabel('Position')

    # Sequence logo for PWM2
    axs[1].set_title('Sequence Logo for PWM2')
    seq_logo(pwm2, axs[1], alphabet)
    axs[1].set_xlabel('Position')

    # Difference logo with p-values
    diff_logo(pwm1, pwm2, axs[2], alphabet, p_values=pvals_corrected, sub_probs=sub_probs)
    axs[2].set_xlabel('Position')

    plt.tight_layout()

    if args.fig_out is not None:
        plt.savefig(args.fig_out)
        print(f"Figure saved to {args.fig_out}")
    else:
        plt.show()

    #Saving the p-values as a CSV file, include column names, round to 2 decimal places:
    if args.pval_csv is not None:
        arr = []
        log_pvals = np.log(pvals_corrected)
        for i,p in enumerate(log_pvals):
            pvals_corrected[i] = np.format_float_scientific(p, precision=2)
            arr.append([i+1, pvals_corrected[i]])
        np.savetxt(args.pval_csv, arr, delimiter=",", fmt='%s', header="position,log-p-value", comments='')

        print(f"log(p-values) saved to {args.pval_csv}")


if __name__ == "__main__":
    main()

