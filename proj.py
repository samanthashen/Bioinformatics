# Write a web-server application using TurboGears to 
# carry out an in silico enzymatic digest of a user-provided protein sequence.
# Users should be able to specify min and max length, 
# min and max molecular weight, # of missed cleavages, and specific enzyme.
# Output should be a table of peptides, with their length, molecular weight, 
## of missed cleavages, and amino-acids to left and right of each peptide 
#in the protein sequence.


def proj(prot, enzyme, maxmisscleave, minlenfprot, maxlenfprot, minweightfprot, maxweightfprot):
    def findcutpos(char, tempprotein, cutpos, terminalside):
        pos = 0
        templen = len(tempprotein)
        #tempprotein is used instead of prot so originial prot is not edited
        #doesnt look at last character in prot because it will be auto added at end of this code
        for letter in tempprotein[:-1]:
            #char is a list of cleave proteins (see get cuts func)
            if letter == char:
                temppos = pos
                #find position of cleave protein
                pos = tempprotein.find(char)+1
                #find protein sequence related to cleavage site
                tempprotein = tempprotein[pos:]
                #if statement to catch first cutpos 0 that is auto added
                if len(cutpos)>0:
                    #find cutside depending on terminalside entered
                    if terminalside == "C":
                        pos = pos + cutpos[-1]
                    else:
                        pos = pos-1 + cutpos[-1]
                #this is where exceptions are taken into account
                #check if there is an exception, if not append cutpos in else
                if redict[enzyme][1] != '':
                    #get triplet surrounding cutpos
                    tripletest = protein[pos+temppos-2:pos+temppos+1]
                    #special rules for trpysin (specific)
                    trypsin = ["CKY", "DKD", "CKH", "CKD", "KKR"
                                "RRH", "RRR", "CRK", "DRD", "RRF", "KRR"]
                    #if the cut pos is after first aa:
                    if len(tripletest)<3:
                        if protein[pos]!="P":
                            if pos not in cutpos:
                                cutpos.append(pos)
                    #checks if aa after cutpos is P, then wil
                    else:
                        if tripletest[2] != "P":
                            if enzyme == "Trypsin (higher specificity)":
                                if tripletest not in trypsin:
                                    if pos not in cutpos:
                                        cutpos.append(pos+temppos)
                            else:
                                if pos not in cutpos:
                                    cutpos.append(pos+temppos)
                else:
                    if pos not in cutpos:
                        cutpos.append(pos)
            
        #if seq doesnt end with char, add final seq num to list
        if len(tempprotein)>0:
            cutpos.append(templen)

    #func to parse cleave pos from expasy
    def rulelist(query, letter, oppositeq, enzcleave):
        
        if cleave.find(query)!=-1:
            terminuspos = cleave.find(query)
            if cleave.find(oppositeq)!=-1:
                cterm = ''.join(cleave[terminuspos+19:cleave.find("C-terminal")-1])
            else:
                cterm = ''.join(cleave[terminuspos+19:])
            if cterm.find("or")!=-1:
                enzcleave[x] = [cterm[0], cterm[5]]
            elif cterm.find("and")!=-1:
                enzcleave[x] = [cterm[0], cterm[3], cterm[9]]
            else:
                if len(cterm)<3:
                    cterm = cterm.replace(',', '')
                enzcleave[x] = cterm.split(", ")

    def getcuts(letter, enzcleave):
        cutpos = [0]
        #try/except to account for no C or N rule
        try:
            for x in enzcleave[enzyme]:
                findcutpos(x, protein, cutpos, letter)
                for y in cutpos:
                    if y not in finallis:
                        finallis.append(y)
                cutpos = [0]
        except KeyError:
            pass

    #table of molecular weights
    mw = {'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07,
          'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08,
          'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10,
          'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06 }

    #before there was error handling in root.py there was this,
    #other characters beside sym would be ignored and incorp
    #not sure if this method is better or not b/c it will accept
    #spaces as per fasta format
    def valid_symbol(sym):
        if sym in 'ACDEFGHIKLMNPQRSTVWY':
            return True
        else:
            return False

    def molWt(seq):
        molwt = 18
        for x in seq:
            if valid_symbol(x) == True:
                molwt = molwt + mw[x]
            else:
                continue
        return molwt

    #handle lowercase
    protein = prot.upper()
    prot = prot.upper()
    
    #make enzymes into a dict
    redict={}
    handle = open('expasydigest.txt',encoding ="ISO-8859-1")
    #use loop to access row
    for line in handle:
        ls = line.split('\t')
        #create value in key for dict by separating values
        #via tabs. 
        enzymes = ls[0]
        cleavage = ls[1]
        exceptions = ls[2].strip('\n')
        redict[enzymes]=[cleavage,exceptions]
    handle.close()

    enzcleavec = {}
    enzcleaven = {}
    #get rules for all enzymes
    for x in redict:
        cleave = redict[x][0]
        if cleave[-1] == '"':
            cleave = cleave[:-1]
        rulelist("C-terminal", "C", ", N-terminal", enzcleavec)
        rulelist("N-terminal", "N", ", C-terminal", enzcleaven)

    #implement getcuts that will also findcutpos
    finallis = []
    getcuts("C", enzcleavec)
    getcuts("N", enzcleaven)

    #sort cutpos by size from smallest to lowest
    finallis = sorted(finallis)
    #print(finallis)

    #print("**************")
    #print(len(redict[enzyme][1]))
    #print(redict[enzyme][1])
    #print("**************")

    #get fragments
    fragments = {}
    for i in range(len(finallis)):
        for j in range(i+1, len(finallis)):
            gi = finallis[i]
            gj = finallis[j]
            if gi == 0 and gj >1:
                fragments[(gi,gj)] = (protein[gj-2], protein[gi:gj], protein[gj:gj+1])
            else:    
                fragments[(gi,gj)] = (protein[gj-2:gj-1], protein[gi:gj], protein[gj:gj+1])

    #print(fragexcept)

    #dict of missed cleavages
    missedcleavage = {}
    miss = 0
    for key in fragments:
        for x in finallis:
            if x in range(key[0]+1,key[1]):
                miss = miss + 1
        missedcleavage[key] = miss
        miss = 0

    table = []
    leftend = {}
    rightend = {}
    #creating table for turbogears
    for x in fragments:
        #min/max/cutpostext is range of fragment in human
        minend = x[0]+1
        maxend = x[1]
        cutpostext = str(minend)+"-"+str(maxend)
        #filters cleavages we entered
        if missedcleavage[x] <= int(maxmisscleave):
            #filters lengths/weights we specifed
            if int(minlenfprot)<len(fragments[x][1])<int(maxlenfprot):
                if int(minweightfprot)<molWt(fragments[x][1])<int(maxweightfprot):
                    #get left and right aa to fragment
                    if minend ==1:
                        leftend[x] = ''
                    else:
                        leftend[x] = protein[minend-2]
                    if maxend == len(protein):
                        rightend[x] = ''
                    else:
                        rightend[x] = protein[maxend]
                    #add everything to table
                    table.append([str(molWt(fragments[x][1])), cutpostext, str(missedcleavage[x]), len(fragments[x][1]), fragments[x][1], leftend[x], rightend[x]])
    return table
