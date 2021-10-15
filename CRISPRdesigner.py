# source file of CDS for set of applicable
h_sequences = open("CDS sequences.txt", 'r').read().splitlines() #split into list - easier to handle

# flat files for forward and reverse genomic sequence, all in upper case
whole_genome = open("whole_genome.txt","r").read()
whole_genome_rc = open("whole_genome_rc.txt","r").read()

all_lengths = []

# clear output file
open("output_file.txt","w").close()

def reverse_complement(to_rc):
    to_c = to_rc[::-1]
    rc = str()
    for base in to_c:
        if base == "A":
            rc = rc+"T"
        if base == "T":
            rc = rc+"A"
        if base == "G":
            rc = rc+"C"
        if base == "C":
            rc = rc+"G"
        if base == "a":
            rc = rc+"t"
        if base == "t":
            rc = rc+"a"
        if base == "g":
            rc = rc+"c"
        if base == "c":
            rc = rc+"g"
    return rc

def Tm_oligo (sequence):
    As = sequence.count("A")
    Ts = sequence.count("T")
    Gs = sequence.count("G")
    Cs = sequence.count("C")
    return 64.9 + 41*(Gs+Cs-16.4)/(As+Ts+Gs+Cs)


# codon usage dictionary
Dic_codon_usage = {'TTT': 0.59, 'TTC': 0.41, 'TTA': 0.28, 'TTG': 0.29, 'CTT': 0.13, 'CTC': 0.06, 'CTA': 0.14, 'CTG': 0.11, 'ATT': 0.46, 'ATC': 0.26, 'ATA': 0.27, 'GTT': 0.39, 'GTC': 0.21, 'GTA': 0.21, 'GTG': 0.19, 'TCT': 0.26, 'TCC': 0.16, 'TCA': 0.21, 'TCG': 0.1, 'AGT': 0.16, 'AGC': 0.11, 'CCT': 0.31, 'CCC': 0.15, 'CCA': 0.41, 'CCG': 0.12, 'ACT': 0.35, 'ACC': 0.22, 'ACA': 0.3, 'ACG': 0.13, 'GCT': 0.38, 'GCC': 0.22, 'GCA': 0.29, 'GCG': 0.11, 'TAT': 0.56, 'TAC': 0.44, 'CAT': 0.64, 'CAC': 0.36, 'CAA': 0.69, 'CAG': 0.31, 'AAT': 0.59, 'AAC': 0.41, 'AAA': 0.58, 'AAG': 0.42, 'GAT': 0.65, 'GAC': 0.35, 'GAA': 0.71, 'GAG': 0.29, 'TGT': 0.63, 'TGC': 0.37, 'CGT': 0.15, 'CGC': 0.06, 'CGA': 0.07, 'CGG': 0.04, 'AGA': 0.48, 'AGG': 0.21, 'GGT': 0.47, 'GGC': 0.19, 'GGA': 0.22, 'GGG': 0.12, 'ATG': 1.0, 'TGG': 1.0, 'TAA': 0.61, 'TAG': 0.09, 'TGA': 0.3}

# edits 'edits.txt' in format: gene name' 'CDS position of edited codon' 'replacing codon
with open('edits.txt') as h_edits:
    line_edit_counter = 0 # for printing edit number in output file
    for line_edit in h_edits:
        line_edit_counter = line_edit_counter +1
        edit_line = line_edit.split()
        h_sequences_index=int(0)
        for sequences_line in h_sequences:
            h_sequences_index = h_sequences_index+1
            #introduce desired edit into each sequence
            if edit_line[0] in sequences_line:
                ori_sequence = (h_sequences[h_sequences_index])
                codon_replaced_sequence = ori_sequence[:int(edit_line[1])-1]+edit_line[2].lower()+ori_sequence[int(edit_line[1])+2:]

        edit = edit_line[0]+"_"+edit_line[1]

        targetname = "targets_"+edit_line[0]
        filename = "%s.txt" % targetname
        print("\nGene: "+edit_line[0]+" Position: "+edit_line[1]+" Edit codon to: "+edit_line[2]+"\n")

        with open(filename) as h_targets:
            overlap_target = "none"
            closest_target = "NA"
            min_pam_dist = 300

            for line_target in h_targets:
                line_target = line_target.split()

                # remove header
                try:
                    int(line_target[0])
                except:
                    continue
                # remove targets with TTTT:
                if line_target[6] == '1':
                    continue
                # test uniqueness of target, skip possible restriction site entries:
                try:
                    hit_20mer = int(line_target[7])
                except:
                    try:
                        hit_20mer = int(line_target[8])
                    except:
                        try:
                            hit_20mer = int(line_target[9])
                        except:
                            print("!!! "+edit_line[0]+"_targets: It seems there are more than 2 restriction sites in this target entry !!!")
                if hit_20mer != 1:
                    continue

                # try to find (core) target with overlap with codon edit 12mer and PAM

                if line_target[2] == '+':
                    core_target = line_target[3][8:23]
                elif line_target[2] == '-':
                    core_target = line_target[3][0:15]

                if codon_replaced_sequence.find(core_target)==-1:
                    overlap_target = line_target
                    chosen_target = line_target
                    edited_sequence = codon_replaced_sequence
                    print("Chosen target: "+str(chosen_target))
                    print("\n***Intended edit results in removal of chosen target***\n")
                    #print("Core target: "+core_target)
                    chosen_core_target = core_target
                    break

                if overlap_target == "none":
                # find target with PAM closest to edit when no target with overlap to intended codon edit
                    if line_target[2] == '+':
                        pam_dist = abs(int(line_target[1])-int(edit_line[1]))
                    elif line_target[2] == '-':
                        pam_dist = abs(int(line_target[0])-int(edit_line[1]))
                    if pam_dist < min_pam_dist:
                        min_pam_dist = pam_dist
                        closest_target = line_target

            # introduce edit to remove target when no target with overlap to intended edit
            if overlap_target == "none" :
                if closest_target == "NA":
                    print("!!!Error - no near target found!!!")
                chosen_target = closest_target
                if closest_target[2] == '+':
                    core_target = closest_target[3][8:23]
                elif closest_target[2] == '-':
                    core_target = closest_target[3][0:15]
                print("Chosen target: "+str(chosen_target))
                print("\n***Introducing additional edits to remove target in donor***\n")
                #print("Core target: "+core_target)
                chosen_core_target = core_target

                # Splits line into chunks of 3 codons each
                codons = [codon_replaced_sequence[i:i+3] for i in range(0, len(codon_replaced_sequence), 3)]
                # find codon number to be altered
                change_codon = int(codon_replaced_sequence.find(core_target)/3+1)

                target_edits = 0 # iteration variable
                codons_to_be_changed = 2
                # introducing chosen number of edits chosen
                while target_edits < codons_to_be_changed:
                    # Assigns all the codons which can code for amino acids
                    PHE = ['TTT', 'TTC']
                    LEU = ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG']
                    ILE = ["ATT", "ATC", "ATA"]
                    MET = ["ATG"]
                    VAL = ["GTT", "GTC", "GTA", "GTG"]
                    SER = ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"]
                    PRO = ["CCT", "CCC", "CCA", "CCG"]
                    THR = ["ACT", "ACC", "ACA", "ACG"]
                    ALA = ["GCT", "GCC", "GCA", "GCG"]
                    TYR = ["TAT", "TAC"]
                    HIS = ["CAT", "CAC"]
                    GLN = ["CAA", "CAG"]
                    ASN = ["AAT", "AAC"]
                    LYS = ["AAA", "AAG"]
                    ASP = ["GAT", "GAC"]
                    GLU = ["GAA", "GAG"]
                    CYS = ["TGT", "TGC"]
                    TRP = ["TGG"]
                    ARG = ["CGT", "CGC", "CGA", "CGG", "AGG", "AGA"]
                    GLY = ["GGT", "GGC", "GGA", "GGG"]
                    # Full list of all possible AAs
                    AA = [PHE,LEU,ILE,MET,VAL,SER,PRO,THR,ALA,TYR,HIS,GLN,ASN,LYS,ASP,GLU,CYS,TRP,ARG,GLY]

                    for sublist in AA:
                        if codons[change_codon] in sublist:
                            if len(sublist) == 1 : # if no synonymous codons go to next
                                change_codon = change_codon+1
                                continue
                            else:
                                synonymous_codon = codons[change_codon] #synonymous_codon here is to be changed
                                frq_to_match = Dic_codon_usage[codons[change_codon]]
                                synonymous_codons = sublist
                                synonymous_codons.remove(codons[change_codon])

                                min_frq_diff = 1
                                for entry in synonymous_codons:
                                    frq_of_syn = Dic_codon_usage[entry]
                                    frq_diff = abs(float(frq_to_match) - float(frq_of_syn))
                                    if frq_diff < min_frq_diff:
                                        min_frq_diff = frq_diff
                                        synonymous_codon = entry

                                while synonymous_codon == codons[change_codon]:
                                    print("!!! Error !!! Codon unchanged !!!")

                                print("Codon changed from: "+codons[change_codon])

                                # make only differences lower case, not whole new codon
                                synonymous_codon_diff = ""
                                for i in range(len(synonymous_codon)):
                                        if synonymous_codon[i] != codons[change_codon][i]:
                                            synonymous_codon_diff = synonymous_codon_diff + synonymous_codon[i].lower()
                                        elif synonymous_codon[i] == codons[change_codon][i]:
                                            synonymous_codon_diff = synonymous_codon_diff + synonymous_codon[i]

                                codons[change_codon]=synonymous_codon_diff
                                print("Codon changed to:   "+codons[change_codon])
                                target_edits = target_edits+1
                                change_codon = change_codon+1
                                break

                edited_sequence = (''.join(codons))

        #print("\n"+"Edited CDS:\n"+edited_sequence)

        ################## Find to be edited region on genome, get genomic flanking regions to donor ##############

        # string length added on each side of edit for search in genome
        len_gsearch = 10
        # string length added from genome to query sequence on each side
        add_genomic = 150 - len_gsearch

        # find edits from left side
        left_sideindex = -1
        for base in edited_sequence:
            if base.isupper():
                left_sideindex = left_sideindex +1
            if base.islower():
                break
        # find edits from right side
        right_sideindex = len(edited_sequence)
        for base in edited_sequence[::-1]:
            if base.isupper():
                right_sideindex = right_sideindex -1
            if base.islower():
                break

        print("")
        around_toedit = ori_sequence[left_sideindex+1-len_gsearch:right_sideindex+len_gsearch]
        print("Query region to edit:   "+around_toedit)
        edit_s = edited_sequence[left_sideindex+1:right_sideindex]
        around_edit = edited_sequence[left_sideindex+1-len_gsearch:right_sideindex+len_gsearch]
        print("Query region with edit: "+around_edit)

        if whole_genome.find(around_toedit) != -1:
            print("***To be edited region found in fw genome sequence***")
            genomic_pos = whole_genome.find(around_toedit)
            donor = whole_genome[genomic_pos-add_genomic:genomic_pos]+around_edit+whole_genome[genomic_pos+len(around_edit):genomic_pos+len(around_edit)+add_genomic]
        elif whole_genome_rc.find(around_toedit) != -1:
            print("***To be edited region found in rv genome sequence***")
            rc_genomic_pos = whole_genome_rc.find(around_toedit)
            donor = whole_genome_rc[rc_genomic_pos-add_genomic:rc_genomic_pos]+around_edit+whole_genome_rc[rc_genomic_pos+len(around_edit):rc_genomic_pos+len(around_edit)+add_genomic]
        else:
            print("!!!Could not find to be edited region query in genome!!!")
            donor = "!!! none !!!"

        print("")
        # create oligo sequences for guide cloning
        if chosen_target[2] == '+':
            target_noPAM = chosen_target[3][0:-3]
        elif chosen_target[2] == '-':
            target_noPAM = reverse_complement(chosen_target[3])[0:-3]
        print("Target without PAM: "+target_noPAM)
        forward_guide_oligo = "gactt"+target_noPAM
        reverse_guide_oligo = "aaac"+reverse_complement(target_noPAM)+"aa"

        print("Forward guide oligo: "+forward_guide_oligo)
        print("Reverse guide oligo: "+reverse_guide_oligo)
        print("")
        print("Donor: "+donor+"\n")

        ################# make primers for donor / genomic amplification

        min_Tm = 55.0

        forward_primer = donor[0]
        reverse_primer = reverse_complement(donor)[0]


        for base in donor[1:]:
            forward_primer += base
            if Tm_oligo(forward_primer ) < min_Tm:
                continue
            else:
                if base == "A" or base == "T":
                    continue

                else:
                    break
        print("Forward primer: "+forward_primer+"\tTm: "+str(round(Tm_oligo(forward_primer),2))+"\tLength: "+str(len(forward_primer)))

        for base in reverse_complement(donor)[1:]:
            reverse_primer += base
            if Tm_oligo(reverse_primer ) < 54.5:
                continue
            else:
                if base == "A" or base == "T":
                    continue
                else:
                    break
        print("Reverse primer: "+reverse_primer+"\tTm: "+str(round(Tm_oligo(reverse_primer),2))+"\tLength: "+str(len(reverse_primer)))

        output_line = open("output_file.txt","a")
        output_line.write(str(line_edit_counter)+"\t"+edit_line[0]+"\t"+edit_line[1]+"\t"+edit_line[2]+"\t"+forward_guide_oligo+"\t"+reverse_guide_oligo+"\t"+donor+"\t"+forward_primer+"\t"+reverse_primer+"\n")
        output_line.close()

        print("\n****************************************************")
