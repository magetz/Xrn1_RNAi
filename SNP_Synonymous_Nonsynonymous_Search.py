import sys, csv, re

Mut_SNP_Gene_in = sys.argv[1]
Mut_SNP_BCF_Contrast_Call_in = sys.argv[2]
Scastellii_ORF_in = sys.argv[3]
Scastellii_Genome_BED_in = sys.argv[4]
Scastellii_Synonymous_SNPs_out = sys.argv[5]
Scastellii_Nonsynonymous_SNPs_out = sys.argv[6]
Scastellii_Nonsynonymous_Frameshift_SNPs_out = sys.argv[7]
Scastellii_Nonsense_SNPs_out = sys.argv[8]
SNPS_Table_out = sys.argv[9]

Mut_SNP_Gene_in_open = csv.reader(open(Mut_SNP_Gene_in,'U'),delimiter="\t")
Mut_SNP_BCF_Contrast_Call_in_open = csv.reader(open(Mut_SNP_BCF_Contrast_Call_in,'U'),delimiter="\t")
Scastellii_ORF_in_open = csv.reader(open(Scastellii_ORF_in,'U'),delimiter="\t")
Scastellii_Genome_BED_in_open = csv.reader(open(Scastellii_Genome_BED_in,'U'),delimiter="\t")
Scastellii_Synonymous_SNPs_out_open = open(Scastellii_Synonymous_SNPs_out,'w')
Scastellii_Nonsynonymous_SNPs_out_open = open(Scastellii_Nonsynonymous_SNPs_out,'w')
Scastellii_Nonsynonymous_Frameshift_SNPs_out_open = open(Scastellii_Nonsynonymous_Frameshift_SNPs_out,'w')
Scastellii_Nonsense_SNPs_out_open = open(Scastellii_Nonsense_SNPs_out, 'w')
SNPS_Table_out_open = open(SNPS_Table_out,'a')

geneticCode = {'GGG':'G', 'GGA':'G', 'GGC':'G', 'GGT':'G',
               'GAG':'E', 'GAA':'E', 'GAC':'D', 'GAT':'D',
               'GCG':'A', 'GCA':'A', 'GCC':'A', 'GCT':'A',
               'GTG':'V', 'GTA':'V', 'GTC':'V', 'GTT':'V',
               'AGG':'R', 'AGA':'R', 'AGC':'S', 'AGT':'S',
               'AAG':'K', 'AAA':'K', 'AAC':'N', 'AAT':'N',
               'ACG':'T', 'ACA':'T', 'ACC':'T', 'ACT':'T',
               'ATG':'M', 'ATA':'I', 'ATC':'I', 'ATT':'I',
               'CGG':'R', 'CGA':'R', 'CGC':'R', 'CGT':'R',
               'CAG':'Q', 'CAA':'Q', 'CAC':'H', 'CAT':'H',
               'CCG':'P', 'CCA':'P', 'CCC':'P', 'CCT':'P',
               'CTG':'L', 'CTA':'L', 'CTC':'L', 'CTT':'L',
               'TGG':'W', 'TGA':'*', 'TGC':'C', 'TGT':'C',
               'TAG':'*', 'TAA':'*', 'TAC':'Y', 'TAT':'Y',
               'TCG':'S', 'TCA':'S', 'TCC':'S', 'TCT':'S',
               'TTG':'L', 'TTA':'L', 'TTC':'F', 'TTT':'F'}

Mut_SNP_Gene_List = []
Mut_SNP_BCF_Contrast_Call_List = []
Scastellii_ORF_List = []
Scastellii_Genome_BED_List = []
Scastellii_WT_ORF_Mut_List = []
Unique_Scastellii_WT_ORF_Name_Mut_List = []
Unique_Scastellii_WT_ORF_Mut_List = []
Scastellii_Mut_ORF_Mut_List = []

for line in Mut_SNP_Gene_in_open:
    Mut_SNP_Gene_List.append(line)

for line in Mut_SNP_BCF_Contrast_Call_in_open:
    Mut_SNP_BCF_Contrast_Call_List.append(line)

for line in Scastellii_ORF_in_open:
    Scastellii_ORF_List.append(line)

for line in Scastellii_Genome_BED_in_open:
    Scastellii_Genome_BED_List.append(line)

for line in Mut_SNP_Gene_List:
    for row in Mut_SNP_BCF_Contrast_Call_List:
        if re.search(line[0], row[0]) and line[1]==row[1]:
            line.append(row[3])

for row in Scastellii_ORF_List:
    for line in Mut_SNP_Gene_List:
        if line[3] == row[0]:
            Scastellii_WT_ORF_Mut_List.append(row[0])

ORF_Set = set(Scastellii_WT_ORF_Mut_List)
Unique_Scastellii_WT_ORF_Name_Mut_List = list(ORF_Set)

for row in Unique_Scastellii_WT_ORF_Name_Mut_List:
    for line in Scastellii_ORF_List:
        if line[0] == row:
            Unique_Scastellii_WT_ORF_Mut_List.append(line)
            

for row in Unique_Scastellii_WT_ORF_Mut_List:
    for line in Mut_SNP_Gene_List:
        Mut_ORF = list(row[1])
        if line[3] == row[0]:
            for column in Scastellii_Genome_BED_List:
                if row[0] == column[3] and column[5] == '+' and len(line[2]) == len(line[4]): #Does point mutations on positive strand
                    Mut_ORF[int(line[1]) - int(column[1]) - int(1)] = line[2]
                    Mut_ORF = "".join(Mut_ORF)
                    Mut_ORF_Num_Seq = [row[0],line[2],line[1],line[4],Mut_ORF]
                    Scastellii_Mut_ORF_Mut_List.append(Mut_ORF_Num_Seq)
                elif row[0] == column[3] and column[5] == '-' and len(line[2]) == len(line[4]): #Does point mutations on negative strand
                    if line[2] == 'A':
                        Mut_ORF[int(column[2]) - int(line[1])] = 'T'
                        Mut_ORF = "".join(Mut_ORF)
                        Mut_ORF_Num_Seq = [row[0],line[2],line[1],line[4],Mut_ORF]
                        Scastellii_Mut_ORF_Mut_List.append(Mut_ORF_Num_Seq)
                    if line[2] == 'T':
                        Mut_ORF[int(column[2]) - int(line[1])] = 'A'
                        Mut_ORF = "".join(Mut_ORF)
                        Mut_ORF_Num_Seq = [row[0],line[2],line[1],line[4],Mut_ORF]
                        Scastellii_Mut_ORF_Mut_List.append(Mut_ORF_Num_Seq)
                    if line[2] == 'C':
                        Mut_ORF[int(column[2]) - int(line[1])] = 'G'
                        Mut_ORF = "".join(Mut_ORF)
                        Mut_ORF_Num_Seq = [row[0],line[2],line[1],line[4],Mut_ORF]
                        Scastellii_Mut_ORF_Mut_List.append(Mut_ORF_Num_Seq)
                    if line[2] == 'G':
                        Mut_ORF[int(column[2]) - int(line[1])] = 'C'
                        Mut_ORF = "".join(Mut_ORF)
                        Mut_ORF_Num_Seq = [row[0],line[2],line[1],line[4],Mut_ORF]
                        Scastellii_Mut_ORF_Mut_List.append(Mut_ORF_Num_Seq)
                elif row[0] == column[3] and len(line[2]) != len(line[4]):
                        Mut_ORF_Num_Seq = [row[0],line[2],line[1],line[4],"INDEL"]
                        Scastellii_Mut_ORF_Mut_List.append(Mut_ORF_Num_Seq)

Scastellii_WT_Prot_Mut_List = []

for item in Unique_Scastellii_WT_ORF_Mut_List:
    protein = ''
    for i in range(0, len(item[1]), 3):
        if re.search("N", item[1][i:i+3]):
            protein += "?"
        elif len(item[1][i:i+3]) < 3:
            protein += "*"
        else:
            codon = item[1][i:i+3]
            protein += geneticCode[codon]   #uses dictionary structure
    Scastellii_WT_Prot_Mut_List.append([item[0],protein])

Scastellii_Mut_Prot_Mut_List = []

for item in Scastellii_Mut_ORF_Mut_List:
    protein = ''
    for i in range(0, len(item[4]), 3):
        if item[4] == "INDEL":
            protein = "INDEL"
        elif re.search("N", item[4][i:i+3]):
            protein += "?"
        elif len(item[4][i:i+3]) < 3:
            protein += "*"
        else:
            codon = item[4][i:i+3]
            protein += geneticCode[codon]   #uses dictionary structure
    Scastellii_Mut_Prot_Mut_List.append([item[0],item[1],item[2],item[3],protein])

Synonymous_Nonsynonymous_Mut_List = []

Total = 0
Frameshift = 0
Synonymous = 0
Nonsynonymous = 0
Nonsense = 0
##
##print Scastellii_WT_Prot_Mut_List
##print Scastellii_Mut_Prot_Mut_List

for line in Scastellii_WT_Prot_Mut_List:
    for row in Scastellii_Mut_Prot_Mut_List:
        if line[0] == row[0]:
            if row[4] == "INDEL":
                Synonymous_Nonsynonymous_Mut_List.append([row[0],row[1],row[2],row[3],row[4],"Frameshift"])
                Scastellii_Nonsynonymous_Frameshift_SNPs_out_open.write(row[0]+"\t"+row[1]+"\t"+row[2]+"\t"+row[3]+"\t"+"Frameshift"+"\n")
                Frameshift += 1
                Total += 1
            elif line[1] == row[4]: 
                Synonymous_Nonsynonymous_Mut_List.append([row[0],row[1],row[2],row[3],row[4],"Synonymous"])
                Scastellii_Synonymous_SNPs_out_open.write(row[0]+"\t"+row[1]+"\t"+row[2]+"\t"+row[3]+"\n")
                Synonymous += 1
                Total += 1
            elif line[1] != row[4]:                                                       
                Synonymous_Nonsynonymous_Mut_List.append([row[0],row[1],row[2],row[3],row[4],"Nonsynonymous"])
                Scastellii_Nonsynonymous_Frameshift_SNPs_out_open.write(row[0]+"\t"+row[1]+"\t"+row[2]+"\t"+row[3]+"\t"+"Nonsynonymous"+"\n")
                Nonsynonymous += 1
                Total += 1

                
SNPS_Table_out_open.write("Exonic SNP Characterization-"+'\n')
SNPS_Table_out_open.write('Total Exonic SNPs:'+" "+str(Total)+'\n')
SNPS_Table_out_open.write('Frameshift SNPs:'+" "+str(Frameshift)+'\n')
SNPS_Table_out_open.write('Synonymous SNPs:'+" "+str(Synonymous)+'\n')
SNPS_Table_out_open.write('Nonsynonymous SNPs:'+" "+str(Nonsynonymous)+'\n')


for line in Synonymous_Nonsynonymous_Mut_List:
    for row in Scastellii_WT_Prot_Mut_List:
        if line[0] == row[0] and line[5] == "Nonsynonymous":
            for i in range(0,len(row[1])):
                if line[4][i] != row[1][i]:
                    position = i + 1
                    Mut_AA = line[4][i]
                    WT_AA = row[1][i]
                    Scastellii_Nonsynonymous_SNPs_out_open.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[5]+"\t"+str(position)+"\t"+WT_AA+"\t"+Mut_AA+"\n")
                    if Mut_AA == "*":
                        Scastellii_Nonsense_SNPs_out_open.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[5]+"\t"+str(position)+"\t"+WT_AA+"\t"+Mut_AA+"\n")
                        Nonsense += 1
                        
SNPS_Table_out_open.write('Nonsense SNPs:'+" "+str(Nonsense)+'\n')
SNPS_Table_out_open.write('\n'+'\n')
            
Scastellii_Synonymous_SNPs_out_open.close()
Scastellii_Nonsynonymous_SNPs_out_open.close()
Scastellii_Nonsense_SNPs_out_open.close()
SNPS_Table_out_open.close()
