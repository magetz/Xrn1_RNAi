import sys, csv, copy, re

Mut_SNP_Gene_in = sys.argv[1]
Scastellii_ORF_in = sys.argv[2]
Exonic_SNPS_out = sys.argv[3]
Intronic_SNPS_out = sys.argv[4]
SNPS_Table_out = sys.argv[5]

Mut_SNP_Gene_in_open = csv.reader(open(Mut_SNP_Gene_in,'U'),delimiter="\t")
Scastellii_ORF_in_open = csv.reader(open(Scastellii_ORF_in,'U'),delimiter="\t")
Exonic_SNPS_out_open = open(Exonic_SNPS_out,'w')
Intronic_SNPS_out_open = open(Intronic_SNPS_out,'w')
SNPS_Table_out_open = open(SNPS_Table_out,'a')

Mut_SNP_Gene_List = []
Intronic_SNPS_List = []

for line in Mut_SNP_Gene_in_open:
    Mut_SNP_Gene_List.append(line)
    Intronic_SNPS_List.append(line)

##print Mut_SNP_Gene_List
##print Intronic_SNPS_List

Scastellii_ORF_List = []
for line in Scastellii_ORF_in_open:
    Scastellii_ORF_List.append(line)

All_SNPS = 0
Exonic_SNPS = 0

Exonic_SNPS_List = []

for line in Mut_SNP_Gene_List:
    All_SNPS+=1
    for row in Scastellii_ORF_List:
        if line[3] == row[0] and re.search(",", row[7]): # Genes with multiple exons
            Exon = re.findall(r'\d+',row[7])
            ExonCount = (len(Exon)/2)
            for i in range(0, int(ExonCount)):
                if line[1] in range(int(Exon[2*i]), int(Exon[2*i+1])):
                    Exonic_SNPS_out_open.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\n")
                    Exonic_SNPS_List.append(line)
                    Exonic_SNPS+=1
        if line[3] == row[0] and not re.search(",", row[7]): # Genes with single exons
            Exonic_SNPS_out_open.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\n")
            Exonic_SNPS_List.append(line)
            Exonic_SNPS+=1

SNPS_Table_out_open.write("Exonic/Intronic Characterization-"+'\n')
SNPS_Table_out_open.write("Total Genic SNPS:"+" "+str(All_SNPS)+'\n')            
SNPS_Table_out_open.write("Exonic SNPS:"+" "+str(Exonic_SNPS)+'\n')

for line in Exonic_SNPS_List:
    for row in Intronic_SNPS_List:
        if line == row:
            Intronic_SNPS_List.remove(line)

for line in Intronic_SNPS_List:
    Intronic_SNPS_out_open.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\n")

SNPS_Table_out_open.write("Intronic SNPS:"+" "+str(len(Intronic_SNPS_List))+'\n')
SNPS_Table_out_open.write('\n'+'\n')

Exonic_SNPS_out_open.close()
Intronic_SNPS_out_open.close()
SNPS_Table_out_open.close()
