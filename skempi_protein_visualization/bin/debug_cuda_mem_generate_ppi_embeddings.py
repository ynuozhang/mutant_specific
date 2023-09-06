import os, sys
import torch
import esm
import pandas as pd
import pickle
import numpy as np
import datetime
from torch.utils.data import TensorDataset
from esm import Alphabet, FastaBatchedDataset
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'max_split_size_mb:32'
import gc
gc.collect()
torch.cuda.empty_cache()
torch.backends.cudnn.benchmark = False
torch.backends.cudnn.deterministic = True

def tprint(string):
    string = str(string)
    sys.stdout.write(str(datetime.datetime.now()) + ' | ')
    sys.stdout.write(string + '\n')
    sys.stdout.flush()

df = pd.read_csv('./skempi_protein_visualization/data/high_confidence_ppi.csv', index_col=0)
columns_to_keep = ['UniProtID_1', 'UniProtID_2', 'symbol_1', 'symbol_2', 'seq_1', 'seq_2',
                   'Experimental System', 'Throughput', 'len_1', 'len_2', 'Pubmed ID']
df = df[columns_to_keep]
df.drop_duplicates(inplace=True) # no duplicate

seqs_wt1 = df.seq_1.values.tolist()
seqs_wt2 = df.seq_2.values.tolist()
seqs_wt1 = set(seqs_wt1)
seqs_wt2 = set(seqs_wt2)
seqs_labeled_wt1 = []
count = 0
for seq in seqs_wt1:
    seqs_labeled_wt1.append(tuple((str('seq' + str(count)), seq)))
    count += 1
seqs_labeled_wt2 = []
count = 0
for seq in seqs_wt2:
    seqs_labeled_wt2.append(tuple((str('seq' + str(count)), seq)))
    count += 1


model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # disables dropout for deterministic results

batch_size = 1000
dataset = FastaBatchedDataset(list(zip(*seqs_labeled_wt1))[0], list(zip(*seqs_labeled_wt1))[1])
batches = dataset.get_batch_indices(batch_size, extra_toks_per_seq=1)
data_loader = torch.utils.data.DataLoader(dataset, collate_fn=Alphabet.from_architecture("roberta_large").get_batch_converter(),
            batch_sampler=batches, pin_memory=True)
dataset_seq2 = FastaBatchedDataset(list(zip(*seqs_labeled_wt2))[0], list(zip(*seqs_labeled_wt2))[1])
batches_seq2 = dataset_seq2.get_batch_indices(batch_size, extra_toks_per_seq=1)
data_loader_seq2 = torch.utils.data.DataLoader(dataset_seq2, collate_fn=Alphabet.from_architecture("roberta_large").get_batch_converter(), batch_sampler=batches_seq2, pin_memory=True)

if torch.cuda.is_available():
    model = model.cuda()
    tprint('Transferred model to GPU')

representation_store_dict = {}
for batch_idx, (labels, strs, toks) in enumerate(data_loader):
    if torch.cuda.is_available():
        toks = toks.to(device='cuda', non_blocking=True)
    with torch.no_grad():
        results = model(toks, repr_layers = [33], return_contacts = True)['representations'][33]
    results_cpu = results.to(device='cpu')
    del results
    tprint('Batch ID: '+str(batch_idx)+str(labels)+str(strs))
    tprint(torch.cuda.memory_allocated())
    tprint(torch.cuda.memory_stats())
    for i, str_ in enumerate(strs):
        # only select representations relate to the sequence
        # rest of the sequences are paddings, check notebook
        # create dictionary {sequence: embeddings}
        representation_store_dict[str_] = results_cpu[i, 1: (len(strs[i])+1)].numpy()
sequence_embeddings = {key: np.mean(value, axis=0, keepdims=True) for key, value in representation_store_dict.items()}

def update_embeddings(row, embedding_dict):
    """
    add embeddings to the metadata column.
    cannot do the reverse, because due to mislabel, several different protein names share the same sequences
    but as long as sequences are correct, so will the embeddings
    """
    for key, value in embedding_dict.items():
        if row == key:
            return value
df['wild_seq_1_embeddings'] = df['wild_seq_1'].apply(update_embeddings, embedding_dict=sequence_embeddings)
df['wild_seq_2_embeddings'] = df['wild_seq_2'].apply(update_embeddings, embedding_dict=sequence_embeddings_seq2)

df.to_hdf('../outputs/dataframes/PPIs_embeddings_meta.hdf', key='df', mode='w')

path = '../outputs/variables/'
if not os.path.exists(path):
    os.makedirs(path)
with open(path+'PPI_seq1_embeddings_full.pk', 'wb') as f:
    pickle.dump(representation_store_dict, f)