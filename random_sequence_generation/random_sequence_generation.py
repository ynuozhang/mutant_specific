import numpy as np
import random
import matplotlib.pyplot as plt
from collections import Counter
from collections import defaultdict
import itertools
from Bio import SeqIO

# you will need biopython to be installed
# set random seed
np.random.seed(42)
random.seed(42)

def draw_aa_composition_barplots(inputseq, plotname, AA, sizes):
    count_aa_composition = []
    for i in inputseq:
        counts = Counter(i)
        counts_sorted = sorted(counts.items(), key=lambda pair:pair[1], reverse=True)
        count_aa_composition.append(counts_sorted)
    merged_list = list(itertools.chain(*count_aa_composition))
    d = defaultdict(int)
    for k, v in merged_list:
        d[k] += v
    aa_composition_check = sorted(d.items(), key=lambda pair:pair[1], reverse=True)
    aa_composition_percentage = [(instance, count / (min(sizes)*len(sizes))) for instance, count in aa_composition_check]
    sorted_original = sorted(AA.items(), key=lambda pair:pair[1], reverse=True)
    import pandas as pd
    AA_composition_df = pd.DataFrame(AA.items(), columns = ['aa', 'freq'])
    AA_composition_df.freq = AA_composition_df.freq/100
    generated_aa_composition_df = pd.DataFrame(aa_composition_percentage, columns = ['aa', 'freq'])
    merged_df = AA_composition_df.merge(generated_aa_composition_df, how='left', on='aa', suffixes=('_from_db', '_from_generation'))
    merged_df['freq_from_generation'] = merged_df['freq_from_generation'].fillna(0)
    merged_df = merged_df.sort_values(by=['freq_from_db'], ascending=False)
    # width, height
    plt.figure(figsize=(25,5))
    plt.suptitle(plotname)
    plt.subplot(1,3,1)
    plt.bar(list(zip(*aa_composition_check))[0], list(zip(*aa_composition_check))[1])
    plt.gca().set_title('Amino Acid Composition in Counts')
    plt.subplot(1,3,2)
    #plt.bar(list(zip(*aa_composition_percentage))[0], list(zip(*aa_composition_percentage))[1])
    plt.bar(merged_df['aa'], merged_df['freq_from_generation'])
    plt.gca().set_title('Amino Acid Composition in Percentage')
    plt.subplot(1,3,3)
    plt.bar(merged_df['aa'], merged_df['freq_from_db'])
    plt.gca().set_title('Amino Acid Composition from Database')
    plt.savefig('%s.png' % plotname)
    plt.show()
    from scipy.spatial.distance import cosine
    from scipy.stats import ks_2samp
    return print('The cosine similarity is %f \n'
                 'The ks test result is %s' % ((1 - cosine(merged_df.freq_from_db, merged_df.freq_from_generation)),
                                               ks_2samp(merged_df.freq_from_db, merged_df.freq_from_generation)))

AA_composition = {'A':9.05, 'R':5.83, 'N':3.79, 'D':5.46,'C':1.28,
                  'Q':3.81, 'E':6.23, 'G':7.27, 'H':2.22, 'I':5.54,
                  'L':9.87, 'K':4.93, 'M':2.34, 'F':3.88, 'P':4.98,
                  'S':6.80, 'T':5.55, 'W':1.30, 'Y':2.88, 'V':6.87,
                  'B':0.025, 'Z':0.025, 'X':0.02,
                  'O':0.025, 'U':0.025
                  }# NO J for esm1b vocabulary

fnames = '/data/rozen/home/e0833634/transformerDENV/data/DENV-E-VIPR.fasta'
sizes = [len(rec) for rec in SeqIO.parse(fnames, 'fasta')]
denv_seqs = [rec.seq for rec in SeqIO.parse(fnames, 'fasta')]
#len(sizes), min(sizes), max(sizes)

def generate_random_seqs(size_, aa, freq):
    return ''.join(np.random.choice(aa, size=size_, replace=True, p=freq))
# random shaffle the lengths of seqs and generate aa seqs based on database composition
random_seqs = []
for size in sizes:
    random_seqs.append(generate_random_seqs(size, np.array(list(AA_composition.keys())), np.array(list(AA_composition.values()))/100))
#random_seqs

# input: 1, list of seqs, 2, name for the plot, 3, AA composition dict from the online database
draw_aa_composition_barplots(denv_seqs, 'DENV_E_aa_composition', AA_composition, sizes)
draw_aa_composition_barplots(random_seqs, 'Random_seqs_aa_composition', AA_composition, sizes)

# write random sequences into fasta
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
random_seqs_copy = random_seqs
for index in range(len(random_seqs_copy)):
    random_seqs_copy[index] = SeqRecord(Seq(random_seqs_copy[index]), id='seq'+str(index))
#random_seqs_copy

SeqIO.write(random_seqs_copy, 'random_seqs_obey_db_distribution.fasta', 'fasta')






