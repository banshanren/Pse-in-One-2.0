"""
Created on Sat May 13 15:35:42 2016
@version:0.2.1./pyc
@author: Fule Liu, Nackel, luo
"""

import sys
import re
import time, math
import multiprocessing

from numpy import array
from itertools import combinations, combinations_with_replacement, permutations, product
import numpy as np

from util import frequency
from util import get_data
from data import index_list
import const


#===========================Kmer===================================================
def make_kmer_list(k, alphabet):
    if k < 0:
        print("Error, k must be an inter and larger than 0.")

    kmers = []
    for i in range(1, k + 1):
        if len(kmers) == 0:
            kmers = list(alphabet)
        else:
            new_kmers = []
            for kmer in kmers:
                for c in alphabet:
                    new_kmers.append(kmer + c)
            kmers = new_kmers

    return kmers


def find_revcomp(sequence, revcomp_dictionary):
    # Save time by storing reverse complements in a hash.
    if sequence in revcomp_dictionary:
        return revcomp_dictionary[sequence]

    # Make a reversed version of the string.
    rev_sequence = list(sequence)
    rev_sequence.reverse()
    rev_sequence = ''.join(rev_sequence)

    return_value = ""
    for letter in rev_sequence:
        if letter == "A":
            return_value += "T"
        elif letter == "C":
            return_value += "G"
        elif letter == "G":
            return_value += "C"
        elif letter == "T":
            return_value += "A"
        elif letter == "N":
            return_value += "N"
        else:
            error_info = ("Unknown DNA character (%s)\n" % letter)
            sys.exit(error_info)
 
    # Store this value for future use.
    revcomp_dictionary[sequence] = return_value
 
    return return_value
 
 
def _cmp(a, b):
    return (a > b) - (a < b)
 
 
def make_revcomp_kmer_list(kmer_list):
    revcomp_dictionary = {}
    new_kmer_list = [kmer for kmer in kmer_list if _cmp(kmer, find_revcomp(kmer, revcomp_dictionary)) <= 0]
    return new_kmer_list
 
 
def make_kmer_vector(k, alphabet, filename, revcomp=False):
    """Generate kmer vector."""
    with open(filename) as f:
        seq_list = get_data(f, alphabet=alphabet)

        if revcomp and re.search(r'[^acgtACGT]', ''.join(alphabet)) is not None:
            sys.exit("Error, Only DNA sequence can be reverse compliment.")
 
        vector = []
        kmer_list = make_kmer_list(k, alphabet)
        for seq in seq_list:
            count_sum = 0
 
            # Generate the kmer frequency dict.
            kmer_count = {}
            for kmer in kmer_list:
                temp_count = frequency(seq, kmer)
                if not revcomp:
                    if kmer not in kmer_count:
                        kmer_count[kmer] = 0
                    kmer_count[kmer] += temp_count
                else:
                    rev_kmer = find_revcomp(kmer, {})
                    if kmer <= rev_kmer:
                        if kmer not in kmer_count:
                            kmer_count[kmer] = 0
                        kmer_count[kmer] += temp_count
                    else:
                        if rev_kmer not in kmer_count:
                            kmer_count[rev_kmer] = 0
                        kmer_count[rev_kmer] += temp_count
 
                count_sum += temp_count
 
            # Normalize.
            if not revcomp:
                count_vec = [kmer_count[kmer] for kmer in kmer_list]
            else:
                revc_kmer_list = make_revcomp_kmer_list(kmer_list)
                count_vec = [kmer_count[kmer] for kmer in revc_kmer_list]
            count_vec = [round(float(e)/count_sum, 8) for e in count_vec]

            vector.append(count_vec)

    return vector
#==============================================================================

#===================IDKmer=====================================================
def make_index_upto_k_revcomp(k):
    """Generate the index for revcomp and from 1 to k."""
    sum = 0
    index = [0]
    for i in range(1, k + 1):
        if i % 2 == 0:
            sum += (math.pow(2, 2 * i - 1) + math.pow(2, i - 1))
            index.append(int(sum))
        else:
            sum += math.pow(2, 2 * i - 1)
            index.append(int(sum))

    return index


def make_index_upto_k(k):
    """Generate the index from 1 to k."""
    sum = 0
    index = [0]
    for i in range(1, k + 1):
        sum += math.pow(4, i)
        index.append(int(sum))

    return index


def make_index(k):
    """Generate the index just for k."""
    index = [0, int(math.pow(4, k))]

    return index


def make_kmer_vector_ID(seq_list, kmer_list, rev_kmer_list, k, upto, revcomp, normalize):
    """Generate kmer vector."""

    # Generate the alphabet index.
    if upto:
        index = make_index_upto_k(k)
        sum = [0] * k
        len_k = k
    else:
        index = make_index(k)
        sum = [0]
        len_k = 1

    vector = []
    for seq in seq_list:
        kmer_count = {}
        # Generate the kmer frequency vector.
        for i in range(len_k):
            sum[i] = 0
            for j in range(index[i], index[i + 1]):
                kmer = kmer_list[j]
                temp_count = frequency(seq, kmer)
                # print temp_count
                if revcomp:
                    rev_kmer = find_revcomp(kmer, {})
                    if kmer <= rev_kmer:
                        if kmer not in kmer_count:
                            kmer_count[kmer] = 0
                        kmer_count[kmer] += temp_count
                    else:
                        if rev_kmer not in kmer_count:
                            kmer_count[rev_kmer] = 0
                        kmer_count[rev_kmer] += temp_count
                else:
                    if kmer not in kmer_count:
                        kmer_count[kmer] = 0
                    kmer_count[kmer] += temp_count
                sum[i] += temp_count

        # Store the kmer frequency vector.
        if revcomp:
            temp_vec = [kmer_count[kmer] for kmer in rev_kmer_list]
        else:
            temp_vec = [kmer_count[kmer] for kmer in kmer_list]

        # Normalize.
        if normalize:
            i = 0
            if not upto:
                temp_vec = [round(float(e)/sum[i], 3) for e in temp_vec]
            if upto:
                if revcomp:
                    upto_index = make_index_upto_k_revcomp(k)
                else:
                    upto_index = make_index_upto_k(k)
                j = 0
                for e in temp_vec:
                    if j >= upto_index[i + 1]:
                        i += 1
                    temp_vec[j] = round(float(e) / sum[i], 3)
                    j += 1

        vector.append(temp_vec)
    # if 0 != len(rev_kmer_list):
    #     print "The kmer is", rev_kmer_list
    # else:
    #     print "The kmer is", kmer_list
    return vector


def diversity(vec):
    """Calculate diversity.

    :param vec: kmer vec
    :return: Diversity(X)
    """
    m_sum = sum(vec)
    from math import log
    return m_sum*log(m_sum, 2) - sum([e*log(e, 2) for e in vec if e != 0])


def id_x_s(vec_x, vec_s, diversity_s):
    """Calculate ID(X, S)

    :param vec_x: kmer X
    :param vec_s: kmer S
    :return: ID(X, S) = Diversity(X + S) - Diversity(X) - Diversity(S)
    """
    # print 'vec_x', vec_x
    # print 'vec_s', vec_s
    vec_x_s = [sum(e) for e in zip(vec_x, vec_s)]
    # print 'vec_x_s', vec_x_s
    # print diversity(vec_x_s), diversity(vec_x), diversity_s
    return diversity(vec_x_s) - diversity(vec_x) - diversity_s


def check_nac_para(k, normalize=False, upto=False, alphabet='ACGT'):
    """Check the nac parameter's validation.
    """
    try:
        if not isinstance(k, int) or k <= 0:
            raise ValueError("Error, parameter k must be an integer and larger than 0.")
        elif not isinstance(normalize, bool):
            raise ValueError("Error, parameter normalize must be bool type.")
        elif not isinstance(upto, bool):
            raise ValueError("Error, parameter upto must be bool type.")
        elif alphabet != 'ACGT':
            raise ValueError("Error, parameter alphabet must be 'ACGT'.")
    except ValueError:
        raise


class IDkmer():
    def __init__(self, k=6, upto=True, alphabet='ACGT'):
        """
        :param k: int, the k value of kmer, it should be larger than 0.
        :param upto: bool, whether to generate 1-kmer, 2-kmer, ..., k-mer.
        :param alphabet: string.
        """
        self.k = k
        self.upto = upto
        self.alphabet = alphabet
        check_nac_para(k=self.k, upto=self.upto, alphabet=self.alphabet)

    def make_idkmer_vec(self, data, hs, non_hs):
        """Make IDKmer vector.

        :param data: Need to processed FASTA file.
        :param hs: Positive FASTA file.
        :param non_hs: Negative FASTA file.
        """

        rev_kmer_list, upto, revcomp, normalize = [], False, False, False

        pos_s_list = get_data(hs, self.alphabet)
        neg_s_list = get_data(non_hs, self.alphabet)
        # print self.k
        if self.upto is False:
            k_list = [self.k]
        else:
            k_list = list(range(1, self.k+1))

        # print 'k_list =', k_list

        # Get all kmer ID from 1-kmer to 6-kmer.
        # Calculate standard source S vector.
        pos_s_vec, neg_s_vec = [], []
        diversity_pos_s, diversity_neg_s = [], []
        for k in k_list:
            kmer_list = make_kmer_list(k, self.alphabet)

            temp_pos_s_vec = make_kmer_vector_ID(pos_s_list, kmer_list, rev_kmer_list, k, upto, revcomp, normalize)
            temp_neg_s_vec = make_kmer_vector_ID(neg_s_list, kmer_list, rev_kmer_list, k, upto, revcomp, normalize)

            temp_pos_s_vec = [sum(e) for e in zip(*[e for e in temp_pos_s_vec])]
            temp_neg_s_vec = [sum(e) for e in zip(*[e for e in temp_neg_s_vec])]

            pos_s_vec.append(temp_pos_s_vec)
            neg_s_vec.append(temp_neg_s_vec)

            diversity_pos_s.append(diversity(temp_pos_s_vec))
            diversity_neg_s.append(diversity(temp_neg_s_vec))

        # Calculate Diversity(X) and ID(X, S).
        sequence_list = get_data(data, self.alphabet)
        vec = []

        for seq in sequence_list:
            # print seq
            temp_vec = []
            for k in k_list:
                kmer_list = make_kmer_list(k, self.alphabet)
                seq_list = [seq]
                kmer_vec = make_kmer_vector_ID(seq_list, kmer_list, rev_kmer_list, k, upto, revcomp, normalize)
                # print 'k', k
                # print 'kmer_vec', kmer_vec

                # print diversity_pos_s
                if upto is False:
                    k = 1

                # print 'pos_vec', pos_s_vec
                # print 'neg_vec', neg_s_vec
                # print 'diversity_pos_s', diversity_pos_s

                temp_vec.append(round(id_x_s(kmer_vec[0], pos_s_vec[k-1], diversity_pos_s[k-1]), 3))
                temp_vec.append(round(id_x_s(kmer_vec[0], neg_s_vec[k-1], diversity_neg_s[k-1]), 3))

            vec.append(temp_vec)

        return vec


def idkmer(k, filename, pos_src_name, neg_src_name):
    idkmer = IDkmer(k=k)
    vec = idkmer.make_idkmer_vec(open(filename), open(pos_src_name), open(neg_src_name))
    return vec





#==============================================================================

#==================getting (k,m)-mismatch profile==============================
def getMismatchProfileMatrix(filename, alphabet, k, m):
    alphabet = list(alphabet)
    p = len(alphabet)
    with open(filename) as f:
        seq_list = get_data(f, alphabet=alphabet)
        kmerdict = getKmerDict(alphabet, k)
        features = []
        if m==0 and m < k:
            for sequence in seq_list:
                vector=getSpectrumProfileVector(sequence, kmerdict, p, k)
                features.append(vector)
        elif m > 0 and m < k:
            for sequence in seq_list:
                vector=getMismatchProfileVector(sequence, alphabet, kmerdict, p, k)
                features.append(vector)   
    return array(features)
        
def getKmerDict(alphabet, k):
    kmerlst = []
    partkmers = list(combinations_with_replacement(alphabet, k))
    for element in partkmers:
        elelst = set(permutations(element, k))
        strlst = [''.join(ele) for ele in elelst]
        kmerlst += strlst
    kmerlst = np.sort(kmerlst)
    kmerdict = {kmerlst[i]:i for i in range(len(kmerlst))}
    return kmerdict


def getSpectrumProfileVector(sequence, kmerdict, p, k):    
    vector = np.zeros((1, p**k))
    n = len(sequence)
    for i in range(n-k+1):
        subsequence=sequence[i:i+k]
        position=kmerdict.get(subsequence)
        vector[0,position] += 1
    return list(vector[0])


def getMismatchProfileVector(sequence, alphabet, kmerdict, p, k): 
    n = len(sequence)
    vector = np.zeros((1, p**k))
    for i in range(n-k+1):
        subsequence = sequence[i:i+k]
        position = kmerdict.get(subsequence)
        vector[0, position]+=1
        for j in range(k):
            substitution = subsequence
            for letter in list(set(alphabet)^set(subsequence[j])):
                substitution = list(substitution)
                substitution[j] = letter
                substitution = ''.join(substitution)
                position = kmerdict.get(substitution)
                vector[0,position] += 1    
    return list(vector[0])

#==============================================================================

#=================getting (k, delta)-subsequence profile=======================
def getSubsequenceProfileByParallel(filename, alphabet, k, delta):
    alphabet = list(alphabet)
    with open(filename) as f:
        seq_list = get_data(f, alphabet=alphabet)
        cpu_num = multiprocessing.cpu_count()   
        batches = constructPartitions(seq_list, cpu_num)
        pool = multiprocessing.Pool(cpu_num)
        results = []
        
        for batch in batches:
            temp=pool.apply_async(getSubsequenceProfile, (batch, alphabet, k, delta))
            results.append(temp)
        pool.close()
        pool.join()
        i = 1
        for temp in results:
            temp_X = temp.get()
            if len(temp_X) != 0:
                if i == 1:
                    X = temp_X
                else:
                    X = np.vstack((X,temp_X))
                i += 1
        return X
    
def constructPartitions(seq_list, cpu_num):
    seqs_num = len(seq_list)
    batch_num = seqs_num//cpu_num
    batches = []
    for i in range(cpu_num-1):
        batch = seq_list[i*batch_num:(i+1)*batch_num]
        batches.append(batch)
    batch=seq_list[(cpu_num-1)*batch_num:]
    batches.append(batch)
    return batches
    
def getSubsequenceProfile(seq_list, alphabet, k, delta):
    kmerdict = getKmerDict(alphabet, k)
    X=[]
    for sequence in seq_list:
        vector = getSubsequenceProfileVector(sequence, kmerdict, k, delta)
        X.append(vector)
    X=array(X)    
    return X

def getSubsequenceProfileVector(sequence, kmerdict, k, delta):      
    vector = np.zeros((1,len(kmerdict)))
    sequence = array(list(sequence))
    n = len(sequence)
    index_lst = list(combinations(range(n), k))
    for subseq_index in index_lst:
        subseq_index = list(subseq_index)
        subsequence = sequence[subseq_index]
        position = kmerdict.get(''.join(subsequence))     
        subseq_length = subseq_index[-1] - subseq_index[0] + 1
        subseq_score = 1 if subseq_length == k else delta**subseq_length    
        vector[0,position] += subseq_score
    #return list(vector[0])
    return [round(f, 4) for f in list(vector[0])]
#==============================================================================
#==============================DR for Protein==================================
def dr_method(inputfile, max_dis):
    if int(max_dis) > 0:
        aa_pairs = make_kmer_list(2, index_list.PROTEIN)
    aa_list = list(index_list.PROTEIN)
    vector_list = []
    with open(inputfile, 'r') as f:
        for line in f:
            vector = []
            if line.strip().startswith('>'):
                continue
            else:
                line = list(line)
                len_line = len(line)
                for i in xrange(max_dis + 1):
                    if i == 0:
                        temp = [line.count(j) for j in aa_list]
                        vector.extend(temp)
                    else:
                        new_line = []
                        for index, elem in enumerate(line):
                            if (index + i) < len_line:
                                new_line.append(line[index] + line[index+i])
                        temp = [new_line.count(j) for j in aa_pairs]
                        vector.extend(temp)
            vector_list.append(vector)
    return vector_list

#==============================================================================
#================================Distance Pair=================================

def get_pseaacdis_dict(raaas_lst, k):
    '''
    :param raaas_lst: A list of reduced amino acid alphabet scheme.
    :param k: The length of pseudo amino acid composition.
    return: a pseaac_dis pattern dictionary.
    '''
    pseaacdis_lst = []
    pseaacdis_dict = {}
    if k == 2:
        part_pseaa = list(combinations_with_replacement(raaas_lst, 2))
        for element in part_pseaa:
            elelst = set(permutations(element, 2))
            pseaacdis_lst += elelst
        pseaacdis_lst.sort()

        for i in range(len(pseaacdis_lst)):
            a, b = pseaacdis_lst[i]
            for j in product(list(a), list(b)):
                pseaacdis_dict[j] = i
        return pseaacdis_dict
    elif k == 1:
        pseaacdis_lst = raaas_lst
        for i in range(len(pseaacdis_lst)):
            for j in list(pseaacdis_lst[i]):
                pseaacdis_dict[j] = i
        return pseaacdis_dict
    else:
        return False


def get_pseaacdis_vector_d(sequence, raaas_lst, distance):
    if distance == 0:
        pseaacdis_dict = get_pseaacdis_dict(raaas_lst, 1)
        sequence = list(sequence)
        vector = np.zeros((1, len(set(pseaacdis_dict.itervalues()))))
        for i in sequence:
            position = pseaacdis_dict.get(i)
            vector[0, position] += 1
        return [round(f,3) for f in vector[0] / sum(vector[0])]
    elif distance > 0:
        pseaacdis_dict = get_pseaacdis_dict(raaas_lst, 2)
        sequence = list(sequence)
        vector = np.zeros((1, len(set(pseaacdis_dict.itervalues()))))
        for i in range(len(sequence) - distance):
            a, b = sequence[i], sequence[i + distance]
            position = pseaacdis_dict.get((a, b))
            vector[0, position] += 1
        return [round(f, 3) for f in vector[0] / sum(vector[0])]
    else:
        return False


def get_pseaacdis_vector(sequence, raaas_lst, max_distance):
    if max_distance >= 0:
        for i in range(max_distance + 1):
            vector_tmp = get_pseaacdis_vector_d(sequence, raaas_lst, i)
            if i == 0:
                vector = vector_tmp
            else:
                vector = np.concatenate((vector_tmp, vector))
        return vector
    else:
        return False


def get_pseaacdis_matrix(filename, reduce_alphabet_scheme, max_distance, alphabet):
    with open(filename) as f:
        seqs = get_data(f, alphabet)
    features = []
    for seq in seqs:
        vector = get_pseaacdis_vector(seq, reduce_alphabet_scheme, max_distance)
        features.append(vector)
    return features

#==============================================================================


def main(args):
    # Set revcomp parameter.
    if args.r != 1:
        args.r = False
    elif args.r == 1 and args.alphabet != 'DNA':
        print("Error, the -r parameter can only be used in DNA.")
    elif args.r == 1 and args.alphabet == 'DNA':
        args.r = True

    # Set alphabet parameter.
    if args.alphabet == 'DNA':
        args.alphabet = index_list.DNA
    elif args.alphabet == 'RNA':
        args.alphabet = index_list.RNA
    elif args.alphabet == 'Protein':
        args.alphabet = index_list.PROTEIN

    if args.method.upper() == 'KMER':
        if args.k is None:
            print "parameters k is required. The default value of k is 2."
            args.k = 2
        if args.r is None:
            print "parameters r is required. The default value of r is 0."
            args.r = 0
        res = make_kmer_vector(k=args.k, alphabet=args.alphabet, filename=args.inputfile, revcomp=args.r)
    elif args.method.upper() == 'IDKMER':
        if args.k is None:
            print "parameters k is required. The default value of k is 6."
            args.k = 6
        if args.ps is None or args.ns is None:
            print 'The positive  and the negative source files are required.'
            return False
        res = idkmer(k=args.k, filename=args.inputfile, pos_src_name=args.ps, neg_src_name=args.ns)
    elif args.method.upper() == "MISMATCH":
        if args.k is None:
            print "parameters k is required. The default value of k is 3."
            args.k = 3
        if args.m is None:
            print "parameters m is required. The default value of m is 1."
            args.m = 1
        if args.m >= args.k:
            print "parameters m should be less than parameter k."
        else:
            res = getMismatchProfileMatrix(args.inputfile, args.alphabet, args.k, args.m)
    elif args.method.upper() == "SUBSEQUENCE":
        if args.delta is None:
            print "parameters delta is required. The default value of delta is 1."
            args.delta = 1
        elif args.delta > 1 or args.delta < 0:
            print "delta should be greater than or equal to 0 and less than or equal to 1."
        if args.k is None:
            print "parameters k is required. The default value of k is 3."
            args.k = 3
        res = getSubsequenceProfileByParallel(filename=args.inputfile, alphabet=args.alphabet, k=args.k, delta=args.delta)
    elif args.method.upper() == 'DR':
        if args.alphabet != index_list.PROTEIN:
            print 'DR method is only available for Protein.'
            return False
        elif args.max_dis < 0 or args.max_dis > 10:
            print 'The max distance can not be negative integer and should be smaller than 11.'
            return False
        else:
            res = dr_method(inputfile=args.inputfile, max_dis=args.max_dis)
            print res
    elif args.method.upper() == 'DP':
        if args.alphabet != index_list.PROTEIN:
            print 'Distance Pair method is only available for Protein.'
            return False
        elif args.max_dis < 0 or args.max_dis > 10:
            print 'The max distance can not be negative integer and should be smaller than 11.'
            return False
        else:
            if args.cp == 'cp_13':
                reduce_alphabet_scheme = const.cp_13
            elif args.cp == 'cp_14':
                reduce_alphabet_scheme = const.cp_14
            elif args.cp == 'cp_19':
                reduce_alphabet_scheme = const.cp_19
            elif args.cp == 'cp_20':
                reduce_alphabet_scheme = const.cp_20
            res = get_pseaacdis_matrix(filename=args.inputfile, reduce_alphabet_scheme=reduce_alphabet_scheme, 
                max_distance=args.max_dis, alphabet=args.alphabet)

    else:
        print("Method error!")

    # Write correspond res file.
    if args.f == 'svm':
        if args.multi == 0 and args.l is None:
            args.l = '+1'
        elif args.multi == 0 and (args.l != '+1' and args.l != '-1'):
            print "For binary classification, the label should be either '+1' or '-1'."
            return False
        elif args.multi == 1 and args.l is None:
            args.l = '0'
        elif args.multi == 1 and args.l is not None:
            try:
                label = int(args.l)
            except ValueError:
                print 'The labels should be integer.'
                return False
        from util import write_libsvm
        write_libsvm(res, [args.l] * len(res), args.outputfile)
    elif args.f == 'tab':
        from util import write_tab
        write_tab(res, args.outputfile)
    elif args.f == 'csv':
        from util import write_csv
        write_csv(res, args.outputfile)


if __name__ == '__main__':
    import argparse
    from argparse import RawTextHelpFormatter

    parse = argparse.ArgumentParser(description="This is kmer module for generate nucleic acid compositio vector.",
                                    formatter_class=RawTextHelpFormatter)
    parse.add_argument('inputfile',
                       help="The input file in FASTA format.")
    parse.add_argument('outputfile',
                       help="The output file stored results.")
    parse.add_argument('alphabet', choices=['DNA', 'RNA', 'Protein'],
                       help="The sequence type.")
    parse.add_argument('method', type=str,
                        help="The method name of nucleic acid composition. {Kmer,mismatch,subsequence}")
    parse.add_argument('-k', type=int,
                       help="For Kmer, IDKmer, mismatch, subsequence methods. The k value of kmer.")
    parse.add_argument('-m', type=int, default=1,
                       help="For mismatch. The max value inexact matching. (m<k)")
    parse.add_argument('-delta', type=float, default=1,
                       help="For subsequence method. The value of penalized factor. (0<=delta<=1)")
    parse.add_argument('-r', type=int, choices=[1, 0],
                       help="Whether consider the reverse complement or not.\n"
                            "1 means True, 0 means False. (default = 0)")
    parse.add_argument('-f', default='tab', choices=['tab', 'svm', 'csv'],
                       help="The output format (default = tab).\n"
                            "tab -- Simple format, delimited by TAB.\n"
                            "svm -- The libSVM training data format.\n"
                            "csv -- The format that can be loaded into a spreadsheet program.")
    parse.add_argument('-l',
                       help="The libSVM output file label.\n"
                       "For binary classification problem, the labels can only be '+1' or '-1'.\n"
                       "For multiclass classification problem, the labels can be set as an integer.")
    parse.add_argument('-ps',
                       help="The input positive source file in FASTA format for IDKmer.")
    parse.add_argument('-ns',
                       help="The input negative source file in FASTA format for IDKmer.")
    parse.add_argument('-max_dis', type=int, default=3,
                       help="For DR and Distance Pair methods. The max distance value of DR and Distance Pair.")
    parse.add_argument('-cp', default='cp_14', choices=['cp_13', 'cp_14', 'cp_19', 'cp_20'],
                       help="For Distance Pair method. The reduced alphabet scheme. Choose one of the four:\n"
                       "cp_13, cp_14, cp_19, cp_20 ")
    parse.add_argument('-multi', type=int, default=0, choices=[0, 1],
                       help="Whether binary classification or multiclass classification.\n"
                       "0: binary classification, default value.\n"
                       "1: multiclass classification.")

    args = parse.parse_args()

    print("Calculating...")
    start_time = time.time()
    main(args)
    print("Done.")
    print("Used time: %ss" % (time.time() - start_time))

#===========================test misnatch==========================================
    #matrix = getMismatchProfileMatrix("test_s.txt", index_list.RNA, 3,0)
#    f=open("test_s.txt")
#    alphabet = index_list.DNA
#    k = 3
#    delta = 1
#    matrix_subseq=getSubsequenceProfileByParallel("test_s.txt", alphabet, k, delta)
#==============================================================================