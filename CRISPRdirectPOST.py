import requests
import time


###################################
h_sequences = open("CDS sequences.txt").read().splitlines() #split into list - easier to handle

with open('edits.txt') as h_edits:
    for line_edit in h_edits:
        edit_pair = line_edit.split()

        h_sequences_index=int(0)
        for h_sequences_line in h_sequences:
            h_sequences_index = h_sequences_index+1
            if edit_pair[0] == h_sequences_line:
                seq = (h_sequences[h_sequences_index])

                targetname = "targets_"+edit_pair[0]
                parameters = {"userseq":seq, "db":"sacCer3", "format":"txt"}
                print(parameters)
                print("****************************************************")

                # get Cas9 target sites from CRISPRdirect
                r = requests.post("http://crispr.dbcls.jp/", data = parameters)
                print(r.text)

                filename = "%s.txt" % targetname
                print(filename+"\n")
                print("****************************************************")
                f = open(filename,"w")
                f.write(r.text)
                f.close()

                time.sleep(5)
