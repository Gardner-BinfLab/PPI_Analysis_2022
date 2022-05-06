import os
import sys
import time
import argparse
from argparse import RawDescriptionHelpFormatter
from itertools import product
from functools import partial
import pandas as pd
import numpy as np
import pickle
import multiprocessing
from multiprocessing import Pool
import threading
from threading import Semaphore
from datetime import datetime
from sklearn.ensemble import RandomForestClassifier
# from sklearn.metrics import roc_auc_score, average_precision_score
from protlearn.features import aaindex1
from protlearn.preprocessing import remove_unnatural



def valid_file(param):
    base, ext = os.path.splitext(param)
    if ext.lower() not in ('.txt', '.tsv', '.csv', '.fasta', '.fa', '.fas', '.faa', '.pickle', '.pkl'):
        raise argparse.ArgumentTypeError('File must have a fasta or csv extension')
    return param



def check_arg(args=None):
    parser = argparse.ArgumentParser(formatter_class=RawDescriptionHelpFormatter, description='Prediction of protein-protein interactions using a random forest classifier trained using AAindex. \n\nexamples:\n\tpython rf.py -a ../test/test_file1.fasta -p ../test/test_pairs.txt -c ../clf/lazypair_clfs.pickle\n\tpython rf.py -a ../test/test_file1.fasta -b ../test/test_file2.fasta -p p -c ../clfs/lazypair_clfs.pickle\n\tpython rf.py -a ../test/test_file1.fasta -p m -c ../clfs/lazypair_clfs.pickle\n\tpython rf.py -a ../test/test_file1.fasta -p s -c ../clfs/lazypair_clfs.pickle')
    
    parser.add_argument('-a', '--seqfile1',
                        type=valid_file,
                        metavar='STR',
                        help='Sequences in fasta format, or tsv/csv format without header.', required='True')
    parser.add_argument('-b', '--seqfile2',
                        type=valid_file,
                        metavar='STR',
                        help='Sequences in fasta format, or tsv/csv format without header. Use it when -p p')
    parser.add_argument('-c', '--classifier',
                        type=valid_file,
                        metavar='STR',
                        help='Classifier file in pickle format', required='True')
    parser.add_argument('-p', '--pair',
                        metavar='STR',
                        help='A file in tsv/csv format (no header) that contains pairs of protein sequence IDs. Use a lowercase p to predict interactions between protein sequences in files from -a and -b using a single process. Use m or s to predict interactions between all protein pairs, using multiprocessing or a single process.', required='True')
    parser.add_argument('-t', '--processes',
                        type=int,
                        metavar='INT',
                        help='Number of processes to spawn. Default = half of the number of CPU.')

    results = parser.parse_args(args)
    return (results.seqfile1, results.seqfile2, results.classifier, results.pair, results.processes)



screen_lock = Semaphore(value=1)
_stop_timer = threading.Event() #global var for thread
def time_count():
    starttime = datetime.now()
    while not _stop_timer.isSet():
        screen_lock.acquire()
        time_message = '\r' + str(datetime.now() - starttime)[:-7]
        sys.stdout.write(time_message)
        sys.stdout.flush()
        screen_lock.release()
        time.sleep(0.01)



def print_time():
    timerthread = threading.Thread(target = time_count, args = ())
    timerthread.start()

    
    
def progress(iteration, total, message=None):
    '''Simple progressbar
    '''
    if message is None:
        message = ''
    bars_string = int(float(iteration) / float(total) * 50.)
    print("\r|%-50s| %d%% (%s/%s) %s "% ('█'*bars_string+ "░" * \
                                     (50 - bars_string), float(iteration)/\
                                     float(total) * 100, iteration, total, \
                                     message), end='\r', flush=True)

    if iteration == total:
        print('\nCompleted!')
        

        
def fasta_reader(file):
    '''Converts .fasta to a pandas dataframe with accession as index
    and sequence in a column 'sequence'
    '''
    fasta_df = pd.read_csv(file, sep='>', lineterminator='>', header=None)
    fasta_df[['Accession', 'Sequence']] = fasta_df[0].str.split('\n', 1, \
                                        expand=True)
    fasta_df['Accession'] = fasta_df['Accession'].str.split('\s').apply(lambda x: x[0])
    fasta_df['Sequence'] = fasta_df['Sequence'].replace('\n', '', regex=True).\
                            astype(str).str.upper().replace('U', 'C').str.replace('\r','')
    total_seq = fasta_df.shape[0]
    fasta_df.drop(0, axis=1, inplace=True)
    fasta_df = fasta_df[(fasta_df.Sequence != '') & (fasta_df.Sequence != 'NONE') & (fasta_df['Sequence'].str.isalpha())]
    fasta_df['Sequence_'] = fasta_df.Sequence.apply(lambda x: len(remove_unnatural(x)))
    fasta_df = fasta_df[fasta_df.Sequence_==1].drop('Sequence_', axis=1)
    fasta_df = fasta_df[(~fasta_df.Sequence.str.contains('X')) & (~fasta_df.Sequence.str.contains('Z'))]
    final_df = fasta_df.dropna()
    remained_seq = final_df.shape[0]
#     if total_seq != remained_seq:
#         print("{} sequences were removed due to inconsistencies in"
#                       "provided file.".format(total_seq-remained_seq))
    return final_df


def predict(clfs, seq, pr, ids, v):
    if type(v) is list:
#         print('Running a single process...')
        v = v

    elif type(v) is str:
#         print('Running multiprocessing...')
        v = [v] 
        pr = pd.DataFrame(product(ids, v))
        pr.columns = ['Protein1','Protein2']
    
    df = pd.merge(pd.merge(pr, seq.rename(columns={'Accession':'Protein1'}), on='Protein1'), seq.rename(columns={'Accession':'Protein2'}), on='Protein2')
    df['AAindex_mean'] = df[['AAindex_x', 'AAindex_y']].values.tolist()
    df['AAindex_mean'] = df['AAindex_mean'].apply(lambda x: np.mean(x, axis=0))
    df['Pairs'] = df[['Protein1','Protein2']].values.tolist()
    df['Pairs'] = df.Pairs.apply(lambda x: ' '.join(sorted(x)))
    preds = pd.DataFrame([clf.predict_proba(df.AAindex_mean.tolist())[:,1] for clf in clfs]).T
    preds['Median'] = preds.median(axis = 1, skipna = True)
    df = pd.concat([df[['Pairs']],preds], axis=1)
    
    return df.to_string(header=False, index=False)
#     filename = o + '.out'
#     pr.to_csv(filename, sep='\t', index=False, encoding='utf-8')
        

def main():
#     startTime = datetime.now()

    clfs = pd.read_pickle(c)
    columns = ['Protein1', 'Protein2', \
               'STRING:full', 'STRING:binding', 'STRING:ptmod', \
               'STRING:activation', 'STRING:reaction', 'STRING:inhibition', \
               'STRING:catalysis', 'STRING:expression', \
               'signor:phosphorylation', 'signor:binding', \
               'signor:transcriptional_regulation', \
               'signor:dephosphorylation', 'signor:cleavage', \
               'signor:ubiquitination', 'signor:relocalization', \
               'signor:guanine_nucleotide_exchange_factor', \
               'signor:gtpase-activating_protein', \
               'signor:post_transcriptional_regulation', 'Median']
    
    base1, ext1 = os.path.splitext(a)
    if ext1.lower() in ('.fasta', '.fa', '.faa', '.fas', '.fna'):
        s1 = fasta_reader(a).reset_index()
    elif ext1.lower() in ('.txt', '.tsv'):
        s1 = pd.read_csv(a, sep='\t', header=None)
    else:
        print('Please provide file in .fasta, .fa, .faa, .fas or .fna extension!')
        sys.exit()
        

           
        
    if os.path.isfile(p) is True:
        base3, ext3 = os.path.splitext(p)
        if bool(ext3):
            if ext3.lower() in ('.txt', '.tsv'):
                pr = pd.read_csv(p, sep='\t', header=None)
                pr.columns = ['Protein1','Protein2']
            elif ext3.lower() in ('.csv'):
                pr = pd.read_csv(p, header=None)
                pr.columns = ['Protein1','Protein2']
            else:
                print('Please provide file in .txt, .tsv or .csv extension!')
                sys.exit()
                    
            pr['Pairs'] = pr[['Protein1','Protein2']].values.tolist()
            pr['Pairs'] = pr.Pairs.apply(lambda x: ' '.join(sorted(x)))

            set_d = set((' '.join(pr.Pairs)).split(' '))
            s1 = pd.merge(s1, pd.DataFrame(set_d).rename(columns={0:'Accession'}), on='Accession')
            aaind1, _ = aaindex1(s1.Sequence.tolist())
            s1['AAindex'] = list(aaind1)
            ids = s1.Accession.tolist()        
            print(' '.join(columns))
            print(predict(clfs, s1, pr, ids, ids))
        else:
            print('Please provide file in .txt, .tsv or .csv extension!')
            sys.exit()

        
    elif p=='p' and os.path.isfile(b) is True:
        base2, ext2 = os.path.splitext(b)
        if ext2.lower() in ('.fasta', '.fa', '.faa', '.fas', '.fna'):
            s2 = fasta_reader(b).reset_index()
        elif ext2.lower() in ('.txt', '.tsv'):
            s2 = pd.read_csv(b, sep='\t', header=None)
        else:
            print('Please provide file in .fasta, .fa, .faa, .fas or .fna extension!')
            sys.exit()

        aaind1, _ = aaindex1(s1.Sequence.tolist())
        s1['AAindex'] = list(aaind1)
        ids = s1.Accession.tolist()
        
        aaind2, _ = aaindex1(s2.Sequence.tolist())
        s2['AAindex'] = list(aaind2)
        id2 = s2.Accession.tolist()
        
        seq = pd.concat([s1,s2])
        
        pr = pd.DataFrame(product(ids, id2))
        pr.columns = ['Protein1','Protein2']
        print(' '.join(columns))
        print(predict(clfs, seq, pr, ids, id2))        

        
    elif p=='m':
        if t>0:
            aaind1, _ = aaindex1(s1.Sequence.tolist())
            s1['AAindex'] = list(aaind1)
            ids = s1.Accession.tolist() 
        
            pools = Pool(t)
            results = []
            pred = partial(predict, clfs, s1, False, ids)
            print(' '.join(columns))

            for result in pools.imap(pred, ids):
                print(result)
            pools.close()
            pools.join()

            
    elif p=='s':
        print(' '.join(['Protein1','Protein2'] + columns))
        aaind1, _ = aaindex1(s1.Sequence.tolist())
        s1['AAindex'] = list(aaind1)
        ids = s1.Accession.tolist()
        
        for i, v in enumerate(ids):
            for j, k in enumerate(ids):
                if i == len(ids):
                    break
                else:

                    aaind_mean = list(np.mean([s1[s1.Accession==v].AAindex.tolist(), s1[s1.Accession==k].AAindex.tolist()], axis=0))
                    pairs = sorted([v,k])
                    preds = pd.DataFrame([clf.predict_proba(aaind_mean)[:,1] for clf in clfs]).T
                    print(pairs[0],pairs[1],preds.to_string(header=False, index=False))
                    
    
    else:
        print('-p takes either a file or a lowercase p, m or s.')
        sys.exit()

#     print('\nPrediction completed in', datetime.now() - startTime, flush = True)    



if __name__ == "__main__":
    a, b, c, p, t = check_arg(sys.argv[1:])
    if t is None:
        t = os.cpu_count()//2
    main()
