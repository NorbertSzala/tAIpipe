#!/usr/bin/python3

'''
Program reads table from tRNAscan-SE output and count previously complemented codons. Output is in format needed in gtAI (R).

Note: input Anti Codon is in DNA (as seq in input tRNAscan-SE)
'''

from collections import OrderedDict
path = '/home/norbert_s/tRNAscan-SE/tRNAscan-SE_postia_output'
# path = '/home/norbert/IBB/skrypty/tRNAscan-SE_postia_output'

def count_tRNA(path):
    translate = {'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'} #dict to translate nts
    not_in_translate = [] #to count how many nucleotides are different than nt in translate
    anti_codons = [] #list of anticodons
    codons_count = {} #dict in format codon : count
    codons = [] #translated anticodons similiar with anti_codon
    codons_order = '''
    1    TTT
    2    TTC
    3    TTA
    4    TTG
    5    TCT
    6    TCC
    7    TCA
    8    TCG
    9    TAT
    10    TAC
    11    TAA
    12    TAG
    13    TGT
    14    TGC
    15    TGA
    16    TGG
    17    CTT
    18    CTC
    19    CTA
    20    CTG
    21    CCT
    22    CCC
    23    CCA
    24    CCG
    25    CAT
    26    CAC
    27    CAA
    28    CAG
    29    CGT
    30    CGC
    31    CGA
    32    CGG
    33    ATT
    34    ATC
    35    ATA
    36    ATG
    37    ACT
    38    ACC
    39    ACA
    40    ACG
    41    AAT
    42    AAC
    43    AAA
    44    AAG
    45    AGT
    46    AGC
    47    AGA
    48    AGG
    49    GTT
    50    GTC
    51    GTA
    52    GTG
    53    GCT
    54    GCC
    55    GCA
    56    GCG
    57    GAT
    58    GAC
    59    GAA
    60    GAG
    61    GGT
    62    GGC
    63    GGA
    64    GGG ''' #in TCAG order

    codons_order = [line.split()[1] for line in codons_order.strip().splitlines()]  #creates list of codons in proper order as a pattern

    with open(path, 'r') as file:
        for line in file:
            anti_codons.append(line.upper().strip().split('\t')[5])

        anti_codons = anti_codons[3:] #deleting header and other rubbish

        for anti_codon in anti_codons: #translating anticodons to codons
            codon = ''
            for nt in anti_codon:
                if nt not in translate:
                    not_in_translate.append(anti_codon)
                    continue  # Goes to next anticodon

                else:
                    codon += translate[nt]
            if len(codon) == 3:
                codons.append(codon[::-1])

        for codon in codons:    #adding amount of codons to dictionary
            if not codon in codons_count:
                codons_count[codon] = 0
            else:
                codons_count[codon] += 1

    for codon in codons_order: # adding missing codons, which are absent in input - total amount of codons should equeals 64
        if codon not in codons_count:
            codons_count[codon] = 0

    sorted_codons_count = dict(OrderedDict((codon, codons_count[codon]) for codon in codons_order if codon in codons_count)) #sort codons in order



    # with open('/home/norbert/IBB/skrypty/count_codons.txt', 'w') as file:
    with open('/home/norbert_s/gtAI/count_codons.txt', 'w') as file:
        file.writelines([str(sorted_codons_count.get(codon, 0)) + '\n' for codon in codons_order])


count_tRNA(path)

