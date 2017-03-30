#!/usr/bin/env python
# -*- coding: utf-8 -*-


METHODS_ALL = ['Kmer', 'RevKmer', 'PseKNC',
               'PseDNC', 'PC-PseDNC-General', 'SC-PseDNC-General', 'PC-PseTNC-General', 'SC-PseTNC-General',
               'PC-PseAAC', 'PC-PseAAC-General', 'SC-PseAAC', 'SC-PseAAC-General',
               'DAC', 'DCC', 'DACC', 'TAC', 'TCC', 'TACC', 'AC', 'CC', 'ACC', 'IDKmer', 'Mismatch', 'Subsequence', 'DR', 'DP']

METHODS_DNA = ['Kmer', 'RevKmer', 'PseKNC', 'PseDNC',
               'PC-PseDNC-General', 'SC-PseDNC-General', 'PC-PseTNC-General', 'SC-PseTNC-General',
               'DAC', 'DCC', 'DACC', 'TAC', 'TCC', 'TACC', 'IDKmer', 'Mismatch', 'Subsequence', 'MAC', 'GAC', 'NMBAC']
METHODS_DNA_KMER = ['Kmer', 'RevKmer', 'IDKmer', 'IDKmer', 'Mismatch', 'Subsequence']
METHODS_DNA_ACC = ['DAC', 'DCC', 'DACC', 'TAC', 'TCC', 'TACC', 'MAC', 'GAC', 'NMBAC']
METHODS_DNA_PSE = ['PseDNC', 'PseKNC',
                   'PC-PseDNC-General', 'SC-PseDNC-General', 'PC-PseTNC-General', 'SC-PseTNC-General',]

METHODS_RNA = ['Kmer', 'Mismatch', 'Subsequence', 'PC-PseDNC-General', 'SC-PseDNC-General', 'DAC', 'DCC', 'DACC', 'MAC', 'GAC', 'NMBAC']
METHODS_RNA_KMER = ['Kmer', 'Mismatch', 'Subsequence']
METHODS_RNA_ACC = ['DAC', 'DCC', 'DACC', 'MAC', 'GAC', 'NMBAC']
METHODS_RNA_PSE = ['PC-PseDNC-General', 'SC-PseDNC-General']

METHODS_PROTEIN = ['Kmer', 'PC-PseAAC', 'PC-PseAAC-General', 'SC-PseAAC', 'SC-PseAAC-General', 'AC', 'CC', 'ACC', 'DR', 'DP']
METHODS_PROTEIN_KMER = ['Kmer', 'DR', 'DP']
METHODS_PROTEIN_ACC = ['AC', 'CC', 'ACC']
METHODS_PROTEIN_PSE = ['PC-PseAAC', 'PC-PseAAC-General', 'SC-PseAAC', 'SC-PseAAC-General']

KMER_FILENAME = 'kmer.py'
ACC_FILENAME = 'acc.py'
PSE_FILENAME = 'pse.py'

METHODS_AC = ['DAC', 'TAC', 'AC']
METHODS_CC = ['DCC', 'TCC', 'CC']
METHODS_ACC = ['DACC', 'TACC', 'ACC']

K_2_DNA_METHODS = ['PseDNC', 'PC-PseDNC-General', 'SC-PseDNC-General', 'DAC', 'DCC', 'DACC']
K_3_DNA_METHODS = ['PC-PseTNC-General', 'SC-PseTNC-General', 'TAC', 'TCC', 'TACC']

THETA_1_METHODS = ['PseDNC', 'PC-PseDNC-General', 'PC-PseTNC-General', 'PC-PseAAC', 'PC-PseAAC-General']
THETA_2_METHODS = ['SC-PseDNC-General', 'SC-PseTNC-General', 'SC-PseAAC', 'SC-PseAAC-General']

DI_INDS_6_DNA = ['Rise', 'Roll', 'Shift', 'Slide', 'Tilt', 'Twist']
TRI_INDS_DNA = ['Dnase I', 'Bendability (DNAse)']

DI_INDS_RNA = ['Rise (RNA)', 'Roll (RNA)', 'Shift (RNA)', 'Slide (RNA)', 'Tilt (RNA)', 'Twist (RNA)']
INDS_3_PROTEIN = ['Hydrophobicity', 'Hydrophilicity', 'Mass']

#reduced alphabet for Distance Pair
cp_13 = ['MF', 'IL', 'V', 'A', 'C', 'WYQHP', 'G', 'T', 'S', 'N', 'RK', 'D', 'E']
cp_14 = ['IMV', 'L', 'F', 'WY', 'G', 'P', 'C', 'A', 'S', 'T', 'N', 'HRKQ', 'E', 'D']
cp_19 = ['P', 'G', 'E', 'K', 'R', 'Q', 'D', 'S', 'N', 'T', 'H', 'C', 'I', 'V', 'W', 'YF', 'A', 'L', 'M']
cp_20 = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


#indices for MAC GAC NMBAC

ALL_DI_DNA_IND = ['Base stacking', 'Protein induced deformability', 'B-DNA twist',
                                 'Dinucleotide GC Content', 'A-philicity', 'Propeller twist',
                                 'Duplex stability-free energy', 'Duplex stability-disrupt energy', 'DNA denaturation',
                                 'Bending stiffness', 'Protein DNA twist', 'Stabilising energy of Z-DNA',
                                 'Aida_BA_transition', 'Breslauer_dG', 'Breslauer_dH', 'Breslauer_dS',
                                 'Electron_interaction', 'Hartman_trans_free_energy', 'Helix-Coil_transition',
                                 'Ivanov_BA_transition', 'Lisser_BZ_transition', 'Polar_interaction', 'SantaLucia_dG',
                                 'SantaLucia_dH', 'SantaLucia_dS', 'Sarai_flexibility', 'Stability', 'Stacking_energy',
                                 'Sugimoto_dG', 'Sugimoto_dH', 'Sugimoto_dS', 'Watson-Crick_interaction', 'Twist',
                                 'Tilt', 'Roll', 'Shift', 'Slide', 'Rise', 'Stacking energy', 'Bend', 'Tip',
                                 'Inclination', 'Major Groove Width', 'Major Groove Depth', 'Major Groove Size',
                                 'Major Groove Distance', 'Minor Groove Width', 'Minor Groove Depth',
                                 'Minor Groove Size', 'Minor Groove Distance', 'Persistance Length',
                                 'Melting Temperature', 'Mobility to bend towards major groove',
                                 'Mobility to bend towards minor groove', 'Propeller Twist', 'Clash Strength',
                                 'Enthalpy', 'Free energy', 'Twist_twist', 'Tilt_tilt', 'Roll_roll', 'Twist_tilt',
                                 'Twist_roll', 'Tilt_roll', 'Shift_shift', 'Slide_slide', 'Rise_rise', 'Shift_slide',
                                 'Shift_rise', 'Slide_rise', 'Twist_shift', 'Twist_slide', 'Twist_rise', 'Tilt_shift',
                                 'Tilt_slide', 'Tilt_rise', 'Roll_shift', 'Roll_slide', 'Roll_rise', 'Slide stiffness',
                                 'Shift stiffness', 'Roll stiffness', 'Rise stiffness', 'Tilt stiffness',
                                 'Twist stiffness', 'Wedge', 'Direction', 'Flexibility_slide', 'Flexibility_shift',
                                 'Entropy']
DEFAULT_DI_DNA_IND = ['Twist', 'Tilt', 'Roll', 'Shift', 'Slide', 'Rise']

ALL_TRI_DNA_IND = ['Bendability-DNAse', 'Bendability-consensus', 'Trinucleotide GC Content',
                                  'Nucleosome positioning', 'Consensus_roll', 'Consensus_Rigid', 'Dnase I',
                                  'Dnase I-Rigid', 'MW-Daltons', 'MW-kg', 'Nucleosome', 'Nucleosome-Rigid']
DEFAULT_TRI_DNA_IND = ['Nucleosome positioning', 'Dnase I']

ALL_RNA_IND = ['Shift', 'Slide', 'Rise', 'Tilt', 'Roll', 'Twist', 'Stacking energy', 'Enthalpy', 'Entropy', 'Free energy', 'Hydrophilicity']
DEFAULT_RNA_IND = ['Shift', 'Slide', 'Rise', 'Tilt', 'Roll', 'Twist']


# For Pse-in-One-Analysis
# Directory of temporary files.
TEMP_DIR = './data/temp/'
TEMP_INDEPENDENT_DIR = './data/temp/independent/'
GEN_FILE_PATH = '/data/gen_files/'
FINAL_RESULTS_PATH = '/data/final_results/'

# Alpabets of DNA, RNA, PROTEIN.
DNA = "ACGT"
RNA = "ACGU"
PROTEIN = "ACDEFGHIKLMNPQRSTVWY"
