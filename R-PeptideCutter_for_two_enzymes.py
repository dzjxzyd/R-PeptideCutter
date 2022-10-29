# author Zhenjiao Du zhenjiao@ksu.edu 
# Grain science & Industry, Kansas State University
# primary intention of the simulated tools
#   1. collect all the unknown enzyme cleavage cites and create a one-stop tools
#   2. no web tools take the hydrolyzing order into consideration

#### a total of 51 enzymes and chemicals
### cleavage posiiton rule
### p4 p3 p2 p1 p1' p2'
### the number in the cleavage was the position number of the last amino acid position (p1 residue position) of the peptide


# the following enzyme cleavage cite are available at
#   https://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html#exceptions
# for the last 16 enzymes, they are available another paper published at JAFC
#       the related tool are available at http://hazralab.iitr.ac.in/new-documentation-01.html


#### input the  enzymes and chemicals functions

# Arg-C proteinase
# the R amino acid was retained by the N terminal amino acid
def Arg_C(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='R':
                cleavage.append(i)
    return cleavage

# ASP-D
# the D amino acid was retanied by the C terminal peptides
def Asp_D(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i+1]=='D':
                cleavage.append(i)
    return cleavage

# Asp-N Endopeptidase + N-terminal Glu:
# the D and E amino acid was retanied by the C terminal peptides
# Asp(D) or Glu(E) in position P1'
def Asp_N(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i+1]=='D' or seq[i+1] == 'E':
                cleavage.append(i)
    return cleavage

# BNPs endopeptidase
def BNPs(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='W':
                cleavage.append(i)
    return cleavage

#Caspase 1
def Caspase1(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-4:
            if seq[i]=='F'or seq[i] == 'W'or seq[i] == 'Y'or seq[i] == 'L':
                if seq[i+2]=='H'or seq[i+2] == 'A'or seq[i+2] == 'T':
                    if seq[i+3]=='D':
                        if seq[i+4]!='P'or seq[i+4]!='E'or seq[i+4]!='D'or seq[i+4]!='Q'or seq[i+4]!='K'or seq[i+4]!='R':
                            cleavage.append(i+3)
    return cleavage


#Caspase 2
def Caspase2(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-4:
            if seq[i]=='D':
                if seq[i+1]=='V':
                    if seq[i+2]=='A':
                        if seq[i+3]=='D':
                            if seq[i+4]!='P'or seq[i+4]!='E'or seq[i+4]!='D'or seq[i+4]!='Q'or seq[i+4]!='K'or seq[i+4]!='R':
                                cleavage.append(i+3)
    return cleavage
#Caspase 3
def Caspase3(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-4:
            if seq[i]=='D':
                if seq[i+1]=='M':
                    if seq[i+2]=='Q':
                        if seq[i+3]=='D':
                            if seq[i+4]!='P'or seq[i+4]!='E'or seq[i+4]!='D'or seq[i+4]!='Q'or seq[i+4]!='K'or seq[i+4]!='R':
                                cleavage.append(i+3)
    return cleavage
#Caspase 4
def Caspase4(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-4:
            if seq[i]=='L' :
                if seq[i+1]=='E':
                    if seq[i+2]=='V':
                        if seq[i+3]=='D':
                            if seq[i+4]!='P'or seq[i+4]!='E'or seq[i+4]!='D'or seq[i+4]!='Q'or seq[i+4]!='K'or seq[i+4]!='R':
                                cleavage.append(i+3)
    return cleavage
#Caspase 5
def Caspase5(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-3:
            if seq[i]=='L' or seq[i]== 'W'  :
                if seq[i+1]=='E':
                    if seq[i+2]=='H':
                        if seq[i+3]=='D':
                            cleavage.append(i+3)
    return cleavage
#Caspase 6
def Caspase6(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-4:
            if seq[i]=='V' :
                if seq[i+1]=='E':
                    if seq[i+2]=='H'or seq[i]== 'I'  :
                        if seq[i+3]=='D':
                            if seq[i+4]!='P'or seq[i+4]!='E'or seq[i+4]!='D'or seq[i+4]!='Q'or seq[i+4]!='K'or seq[i+4]!='R':
                                cleavage.append(i+3)
    return cleavage
#Caspase 7
def Caspase7(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-4:
            if seq[i]=='D' :
                if seq[i+1]=='E':
                    if seq[i+2]=='V'  :
                        if seq[i+3]=='D':
                            if seq[i+4]!='P'or seq[i+4]!='E'or seq[i+4]!='D'or seq[i+4]!='Q'or seq[i+4]!='K'or seq[i+4]!='R':
                                cleavage.append(i+3)
    return cleavage
#Caspase 8
def Caspase8(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-4:
            if seq[i]=='I' or seq[i]=='L':
                if seq[i+1]=='E':
                    if seq[i+2]=='T':
                        if seq[i+3]=='D':
                            if seq[i+4]!='P'or seq[i+4]!='E'or seq[i+4]!='D'or seq[i+4]!='Q'or seq[i+4]!='K'or seq[i+4]!='R':
                                cleavage.append(i+3)
    return cleavage
#Caspase 9
def Caspase9(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-3:
            if seq[i]=='L' :
                if seq[i+1]=='E':
                    if seq[i+2]=='H':
                        if seq[i+3]=='D':
                            cleavage.append(i+3)
    return cleavage
#Caspase 10
def Caspase10(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-3:
            if seq[i]=='I' :
                if seq[i+1]=='E':
                    if seq[i+2]=='A':
                        if seq[i+3]=='D':
                            cleavage.append(i+3)
    return cleavage


#Chymotrypsin-high specificity(C-term to [FYW], not before P)
def chymo_high(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='F' or seq[i] == 'Y':
                if seq[i+1] != 'P':
                    cleavage.append(i)
            if seq[i]=='W' :
                if seq[i+1] != 'P' and seq[i] != 'M':
                    cleavage.append(i)
    return cleavage

#Chymotrypsin-low specificity (C-term to [FYWML], not before P)
def chymo_low(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='F' or seq[i] == 'L' or seq[i] == 'Y':
                if seq[i+1] != 'P':
                    cleavage.append(i)
            if seq[i]=='W' :
                if seq[i+1] != 'P' and seq[i] != 'M':
                    cleavage.append(i)
            if seq[i]=='M' :
                if seq[i+1] != 'P' and seq[i] != 'Y':
                    cleavage.append(i)
            if seq[i]=='H' :
                if seq[i+1] != 'D' and seq[i] != 'M' and seq[i] != 'P' and seq[i] != 'W':
                    cleavage.append(i)
    return cleavage


#Clostripain (Clostridiopeptidase B)
def Clostripain(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='R':
                cleavage.append(i)
    return cleavage
#CNBr
def CNBr(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='M':
                cleavage.append(i)
    return cleavage
#Enterokinase
def Enterokinase(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-3:
            if seq[i]=='D' or seq[i]== 'E':
                if seq[i+1]=='D' or seq[i+1]== 'E':
                    if seq[i+2]=='D' or seq[i+2]== 'E':
                        if seq[i+3]=='K':
                            cleavage.append(i+3)
    return cleavage

#Factor_xa
#A,F,G,I,L,T,V or M
def Factor_xa(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-3:
            if seq[i]=='A' or seq[i]== 'F'or seq[i]== 'G'or seq[i]== 'I'or seq[i]== 'L'or seq[i]== 'T'or seq[i]== 'V'or seq[i]== 'M':
                if seq[i+1]=='D' or seq[i+1]== 'E':
                    if seq[i+2]=='G':
                        if seq[i+3]=='R':
                            cleavage.append(i+3)
    return cleavage
#Formic_acid
def Formic_acid(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='D':
                cleavage.append(i)
    return cleavage
#Glutamyl_endopeptidase
def Glutamyl_endopeptidase(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='E':
                cleavage.append(i)
    return cleavage
#GranzymeB
def GranzymeB(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-3:
            if seq[i]=='I':
                if seq[i+1]== 'E':
                    if seq[i+2]=='P':
                        if seq[i+3]=='D':
                            cleavage.append(i+3)
    return cleavage
#Hydroxylamine
def Hydroxylamine(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='N' and seq[i+1] == 'G':
                cleavage.append(i)
    return cleavage
#Iodosobenzoic_acid
def Iodosobenzoic_acid(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='W':
                cleavage.append(i)
    return cleavage
#LysC
def LysC(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='K':
                cleavage.append(i)
    return cleavage
#Neutrophil_elastase
def Neutrophil_elastase(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='A' or seq[i]== 'V':
                cleavage.append(i)
    return cleavage
#NTCB (2-nitro-5-thiocyanobenzoic acid)
def NTCB(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i+1]=='C':
                cleavage.append(i)
    return cleavage
#Pepsin (pH1.3)
def Pepsin_1(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i+1]== 'L' or seq[i+1]== 'F':
                cleavage.append(i)
                if i>=2 and i<seq_len-1:
                    if seq[i-1]=='P'or seq[i-2] == 'H' or  seq[i-2] == 'K' or  seq[i-2] == 'R' or seq[i-2]=='P':
                        cleavage.remove(i)
                if i ==1:
                    if seq[i-1]=='P':
                        cleavage.remove(i)
                if i==seq_len-1:
                    if seq[i-1]=='R':
                        cleavage.remove(i)

            if seq[i]== 'L' or seq[i]== 'F':
                if i not in cleavage: # when meet ASLLAL the double L will cause the repeat postion extraction
                    cleavage.append(i)
                if i>=2 and i<seq_len-1:
                    if seq[i-1]=='P'or seq[i-2] == 'H' or  seq[i-2] == 'K' or  seq[i-2] == 'R' or seq[i-2]=='P':
                        cleavage.remove(i)
                if i ==1:
                    if seq[i-1]=='P':
                        cleavage.remove(i)
                if i==seq_len-1:
                    if seq[i-1]=='R':
                        cleavage.remove(i)

    return cleavage

#Pepsin (pH>2)
def Pepsin_2(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i+1]== 'L' or seq[i+1]== 'F' or seq[i+1]== 'Y' or seq[i+1]== 'W':
                cleavage.append(i)
                if i>=2 and i<seq_len-1:
                    if seq[i-1]=='P'or seq[i-2] == 'H' or  seq[i-2] == 'K' or  seq[i-2] == 'R' or seq[i-2]=='P':
                        cleavage.remove(i)
                if i ==1:
                    if seq[i-1]=='P':
                        cleavage.remove(i)
                if i==seq_len-1:
                    if seq[i-1]=='R':
                        cleavage.remove(i)
            if seq[i]== 'L' or seq[i]== 'F'or seq[i]== 'Y' or seq[i]== 'W':
                if i not in cleavage: # when meet ASLLAL the double L will cause the repeat postion extraction
                    cleavage.append(i)
                if i>=2 and i<seq_len-1:
                    if seq[i-1]=='P'or seq[i-2] == 'H' or  seq[i-2] == 'K' or  seq[i-2] == 'R' or seq[i-2]=='P':
                        cleavage.remove(i)
                if i ==1:
                    if seq[i-1]=='P':
                        cleavage.remove(i)
                if i==seq_len-1:
                    if seq[i-1]=='R':
                        cleavage.remove(i)

    return cleavage
#Proline_endopeptidase
def Proline_endopeptidase(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-2:
            if seq[i]=='H' and seq[i]== 'K'and seq[i]== 'R':
                if seq[i+1]=='P' :
                    if seq[i+2]!='P':
                        cleavage.append(i+1)
    return cleavage
#Proteinase_K
#A,E,F,I,L,T,V,W or Y
def Proteinase_K(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='A' or seq[i]=='E' or seq[i]== 'F'or seq[i]== 'I'or seq[i]== 'L'or seq[i]== 'T'or seq[i]== 'V'or seq[i]== 'W'or seq[i]== 'Y':
                cleavage.append(i)
    return cleavage
#Staphylococcal_peptidase_I
def Staphylococcal_peptidase_I(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if i==0:
                if seq[i]== 'E':
                    cleavage.append(i)
            if seq[i]!='E':
                if seq[i+1]== 'E':
                    cleavage.append(i+1)
    return cleavage

#Thermolysin
def Thermolysin(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]!='E'and seq[i]!= 'D':
                if seq[i+1]== 'A' or seq[i+1] == 'F'or seq[i+1] == 'I'or seq[i+1] == 'L'or seq[i+1] == 'M'or seq[i+1] == 'V':
                    cleavage.append(i)
    return cleavage
#Thrombin
def Thrombin(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-2:
            if seq[i]=='G':
                if seq[i+1] == 'R':
                    if seq[i+2] =='G':
                        cleavage.append(i+1)
        if i < seq_len-4:
            if seq[i]== 'A' or seq[i] == 'F'or seq[i] == 'G'or seq[i] == 'I'or seq[i] == 'L'or seq[i] == 'T'or seq[i] == 'V' or seq[i] == 'M':
                if seq[i+1]== 'A' or seq[i+1] == 'F'or seq[i+1] == 'G' or seq[i+1] == 'I'or seq[i+1] == 'L'or seq[i+1] == 'T'or seq[i+1] == 'V' or seq[i+1] == 'W':
                    if seq[i+2]=='P':
                        if seq[i+3]=='R':
                            if seq[i+4]!='D' and seq[i+4] != 'E':
                                cleavage.append(i+3)
                            if i+5 == seq_len-1: # position index is from 0 instead of 1
                                if seq[i+5]=='D' and seq[i+5] == 'E':
                                    cleavage.remove(i+3)
    return cleavage

#Trypsin (please note the exceptions!)
def Trypsin(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='K' or seq[i] =='R':
                if seq[i+1] != 'P':
                    cleavage.append(i)
        if i < seq_len-2:
            if seq[i] == 'W':
                if seq[i+1] == 'K':
                    if seq[i+2] =='P':
                        cleavage.append(i+1)
        if i < seq_len-2:
            if seq[i] =='M':
                if seq[i+1] == 'R':
                    if seq[i+2] =='P':
                        cleavage.append(i+1)
        if i < seq_len-2:
            if seq[i] =='E':
                if seq[i+1] == 'R':
                    if seq[i+2] =='P':
                        cleavage.append(i+1)
        #situation will block the cleavage
        if i < seq_len-2:
            if seq[i]=='D' or seq[i]=='C':
                if seq[i+1] == 'K':
                    if seq[i+2] == 'D':
                        if i+1 in cleavage:
                            cleavage.remove(i+1)
        if i < seq_len-2:
            if seq[i]=='C':
                if seq[i+1] == 'K':
                    if seq[i+2] == 'H' or seq[i+2] == 'Y':
                        if i+1 in cleavage:
                            cleavage.remove(i+1)
        if i < seq_len-2:
            if seq[i]=='C':
                if seq[i+1] == 'R':
                    if seq[i+2] == 'K':
                        if i+1 in cleavage:
                            cleavage.remove(i+1)
        if i < seq_len-2:
            if seq[i]=='R':
                if seq[i+1] == 'R':
                    if seq[i+2] == 'H' or seq[i+2] =='R':
                        if i+1 in cleavage:
                            cleavage.remove(i+1)
    return cleavage


# only pepsin and trypsin are following the cleavage cites data from PeptideCutter
# obtanined from http://hazralab.iitr.ac.in/new-documentation-01.html
# Elastase_1
def Elastase_1(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='A':
                cleavage.append(i)
    return cleavage

# Elastase_2 (Neutrophil_elastase)
def Elastase_2(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='A' or seq[i]=='V':
                cleavage.append(i)
    return cleavage

# chymotrypsinogen_B1
def chymotrypsinogen_B1(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='F' or seq[i]=='L'or seq[i]=='W'or seq[i]=='Y':
                cleavage.append(i)
    return cleavage
# chymotrypsinogen_C
def chymotrypsinogen_C(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='F' or seq[i]=='M' or seq[i]=='L'or seq[i]=='W'or seq[i]=='Y' or seq[i]=='N'or seq[i]=='Q':
                cleavage.append(i)
    return cleavage

# pancreatic_enteropeptidase_E
def pancreatic_enteropeptidase_E(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='A':
                cleavage.append(i)
    return cleavage

# Enteropeptidase   ///this protein need to take the p1' into consideration
def Enteropeptidase(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='K' and seq[i+1] =='I':
                cleavage.append(i)
    return cleavage
# prostasin  ///highly specificity
def prostasin(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-7:
            if seq[i]=='R':
                if seq[i+1]=='K':
                    if seq[i+2]=='R':
                        if seq[i+3]=='K':
                            if seq[i+4]=='I':
                                if seq[i+5]=='S':
                                    if seq[i+6]=='G':
                                        if seq[i+7]=='K':
                                            cleavage.append(i+3)
    return cleavage

#Gastricsin
def Gastricsin(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='Y':
                cleavage.append(i)
    return cleavage

# Fruit_Bromelain
def Fruit_Bromelain(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-2:
            if seq[i]=='F':
                if seq[i+1]=='V':
                    if seq[i+2]=='R':
                        cleavage.append(i+2)
    return cleavage
#Stem_Bromelain
def Stem_Bromelain(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='R':
                if seq[i+1]=='R':
                    cleavage.append(i+1)
    return cleavage
# Ananain
def Ananain(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-2:
            if seq[i]=='F':
                if seq[i+1]=='V':
                    if seq[i+2]=='R':
                        cleavage.append(i+2)
    return cleavage
# Papaya_proteinase
def Papaya_proteinase(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='G':
                cleavage.append(i)
    return cleavage

# Chymopapain
def Chymopapain(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]== 'A' or seq[i] == 'V'or seq[i] == 'L'or seq[i] == 'I'or seq[i] == 'F'or seq[i] == 'W':
                if seq[i+1] == 'R' or seq[i+1] == 'L':
                    cleavage.append(i+1)
    return cleavage

# Chymosin
def Chymosin(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]=='S' or seq[i]=='F':
                if seq[i+1]=='M'or seq[i]=='A':
                    cleavage.append(i)
    return cleavage

#Caricain
#same to the chymopapain
def Caricain(seq,seq_len):
    cleavage=[]
    for i in range(seq_len):
        if i < seq_len-1:
            if seq[i]== 'A' or seq[i] == 'V'or seq[i] == 'L'or seq[i] == 'I'or seq[i] == 'F'or seq[i] == 'Y' or seq[i] == 'W':
                if seq[i+1] == 'R' or seq[i+1] == 'L':
                    cleavage.append(i+1)
    return cleavage

#### auto_input the protien sequence
import pandas as pd
import numpy as np
import collections

import os
os.getcwd()
os.chdir('/Users/zhenjiaodu/Desktop/wheat bran study/reference proteomics/sequence info collection 52 sequence')
os.getcwd()
import glob
# read and storage all the fasta file name in my folder
fasta_name = glob.glob("*.fasta")
sum_accession_number=[]
sum_seq=[]

for file_name in fasta_name:
    #print(file_name)
    fasta = open(file_name, mode='r')
    temp =[] # create a list for the fasta content temporaty storage

    for line in fasta:
        temp.append(line.strip()) # delete all the /n in the original fasta file

    temp_seq=str() # collect the seq info
    temp_name=str() #collect the seq name
    for i in range(len(temp)):
        if i == 0:
            temp_name=temp[0].split('|')[1] # get the first line and extract the second values (accession_number)
        if i >= 1:
            temp_seq= temp_seq+temp[i]

    sum_accession_number.append(temp_name)
    sum_seq.append(temp_seq)

# running the 51 enzymes and chemicals below
#   then running the options
#       options contanined the 47 functions
options= [Arg_C,Asp_D,Asp_N,BNPs,Caspase1,Caspase2,Caspase3,Caspase4,Caspase5,Caspase6,Caspase7,Caspase8,Caspase9,Caspase10,chymo_high,chymo_low,Clostripain,CNBr,Enterokinase,Factor_xa,Formic_acid,Glutamyl_endopeptidase,GranzymeB,Hydroxylamine,Iodosobenzoic_acid,LysC,Neutrophil_elastase,NTCB,Pepsin_1,Pepsin_2,Proline_endopeptidase,Proteinase_K,Staphylococcal_peptidase_I,Thermolysin,Thrombin,Trypsin,Elastase_1,Elastase_2, chymotrypsinogen_B1, chymotrypsinogen_C, pancreatic_enteropeptidase_E, Enteropeptidase,prostasin,Gastricsin,Fruit_Bromelain,Stem_Bromelain,Ananain,Papaya_proteinase,Chymopapain,Chymosin, Caricain]
# assign the function name to option_names for further file generating
options_names= ['Arg_C','Asp_D','Asp_N','BNPs','Caspase1','Caspase2','Caspase3','Caspase4','Caspase5','Caspase6','Caspase7','Caspase8','Caspase9','Caspase10','chymo_high','chymo_low','Clostripain','CNBr','Enterokinase','Factor_xa','Formic_acid','Glutamyl_endopeptidase','GranzymeB','Hydroxylamine','Iodosobenzoic_acid','LysC','Neutrophil_elastase','NTCB','Pepsin_1','Pepsin_2','Proline_endopeptidase','Proteinase_K','Staphylococcal_peptidase_I','Thermolysin','Thrombin','Trypsin','Elastase_1','Elastase_2', 'chymotrypsinogen_B1', 'chymotrypsinogen_C', 'pancreatic_enteropeptidase_E', 'Enteropeptidase','prostasin','Gastricsin','Fruit_Bromelain','Stem_Bromelain','Ananain','Papaya_proteinase','Chymopapain','Chymosin', 'Caricain']

#desired peptide length
len_peptide = 2
# Create a Pandas Excel writer using XlsxWriter as the engine. for the data collection at different sheet_name
writer = pd.ExcelWriter('two enzyme results.xlsx', engine='xlsxwriter')

for name in range(len(sum_seq)):
    seq = sum_seq[name]
    seq_len = len(seq)
    accession_number = sum_accession_number[name]
    summary = np.empty([150,51*51*3],dtype=object) #approximate the final size of data
    times = 0
    for j in range(len(options)):
        two_cleavage = []
        two_cleavage = options[j](seq,seq_len)

        two_cleavage.insert(0,0)
        two_cleavage.insert(len(two_cleavage),seq_len-1) # only to make it easy when extrac peptides below (we need to add 1 extract from the P1' instead of P1)
        for m in range(len(options)):
            two_pep_seq = [] # used to ccollect peptide sequences from the second hydrolysis
            second_cleavage=[]
            for i in range(len(two_cleavage)-1):
                # extract the sequence from the original protein sequences
                # the cleavage position was at the P1, so we should extract from P1'
                if i ==0:
                    # i ==0 means the sequence from zero instead of the 1 (different from the other position in the sequence (need to add 1))
                    pep_seq_temp = seq[two_cleavage[i]:two_cleavage[i+1]+1]
                    pep_seq_temp
                    len_pep_seq_temp = len(pep_seq_temp)
                else:
                    pep_seq_temp = seq[two_cleavage[i]+1:two_cleavage[i+1]+1]
                    pep_seq_temp
                    len_pep_seq_temp = len(pep_seq_temp)
                # find the cleavage P1 cite position
                temp_list= options[m](pep_seq_temp,len_pep_seq_temp)
                # the the second peptide sequence the index is not from 0(in the original sequence)
                # we need to add the initial index two_cleavage[i]
                # why there is a additional plus 1, because the temp_list is from 0 adn two_cleavage is also from 0
                new_list = [x+two_cleavage[i]+1 for x in temp_list]

                #asign the initial position of the peptide sequences
                new_list.insert(0,two_cleavage[i])
                # assing the last index of the peptide sequences
                new_list.append(two_cleavage[i+1])
                second_cleavage = second_cleavage+new_list
                new_list
            second_cleavage
            second_cleavage.insert(0,0)
            second_cleavage.insert(len(second_cleavage),seq_len-1)
            second_cleavage
            for n in range(len(second_cleavage)-1):
                if n ==0:
                    second_pep_seq_temp = seq[second_cleavage[n]:second_cleavage[n+1]+1] # same reason as above
                else:
                    second_pep_seq_temp = seq[second_cleavage[n]+1:second_cleavage[n+1]+1] # same reason as above
                #print(second_pep_seq_temp)
                if len(second_pep_seq_temp) == len_peptide: # check the peptide length
                    two_pep_seq.append(second_pep_seq_temp)
            two_pep_seq
            # delete the duplicated
            temp_list=[]
            temp_count=[]
            for i in two_pep_seq:
                if i not in temp_list:
                    temp_list.append(i) # collect the unduplicate peptide info
                    temp_count.append(1) #count as 1 for the occurance of a peptide
                else:
                    index_repeat=temp_list.index(i) # find the position of the repeated peptide sequences
                    temp_count[index_repeat]=temp_count[index_repeat]+1


            if len(temp_list) != 0:
                #temp_list
                # match the theoretical peptide sequence info and their number generated from the proteins under the enzymatic hydrolysis
                # match
                results = collections.defaultdict(list)
                for i in range(len(temp_list)):
                    results[i].append(temp_list[i])
                    results[i].append(temp_count[i])
                #results export
                df = pd.DataFrame(results)
                df = df.T
                df = df.sort_values(by = df.columns[1],ascending=None) # column[1] asend based on the count number
                ky = df.to_numpy()

                #assign the well arrange data into the summary numpy
                for ky_row in range(ky.shape[0]+1):
                    for ky_column in range(2):
                        if ky_row==0: # assign the number of the column
                            if ky_column==0:
                                summary[ky_row,times*2+ky_column]=options_names[j]+'+'+options_names[m]
                                options_names[j]
                            if ky_column==1:
                                summary[ky_row,times*2+ky_column]='theoretical peptide amounts'
                        if ky_row>0:
                            summary[ky_row,times*2+ky_column]=ky[ky_row-1,ky_column]  # column-1, so index back to the original data
                times = times +1


    summary= pd.DataFrame(summary)
    summary_filename= accession_number
    summary.to_excel(writer, sheet_name=summary_filename)
writer.save()
