import sys, csv, copy

Mutant_SNPS_in = sys.argv[1]
Gene_Annotations_in = sys.argv[2]
Genic_SNPs_out = sys.argv[3]
Intergenic_SNPs_out = sys.argv[4]
SNPS_Table_out = sys.argv[5]

Mutant_SNPS_in_open = csv.reader(open(Mutant_SNPS_in,'U'),delimiter="\t")
Gene_Annotations_in_open = csv.reader(open(Gene_Annotations_in,'U'),delimiter="\t")
Genic_SNPs_out_open = open(Genic_SNPs_out,'w')
Intergenic_SNPs_out_open = open(Intergenic_SNPs_out,'w')
SNPS_Table_out_open = open(SNPS_Table_out,'a')

Mutant_SNPS_list = [] #If you do not do this, it only iterates through the list once. Property of csv reader.
for line in Mutant_SNPS_in_open:
    Mutant_SNPS_list.append(line)

Intergenic_SNPS_list = [] #Making initial intergenic SNPs list the same as the mutant SNP list.
for line in Mutant_SNPS_list:
    Intergenic_SNPS_list.append(line)
    
Gene_Annotations_list = [] #If you do not do this, it only iterates through the list once. Property of csv reader. 
for line in Gene_Annotations_in_open:
    Gene_Annotations_list.append(line)

All_SNPs = 0
Genic_SNPs = 0
Intergenic_SNPs_below = 0
Intergenic_SNPs_above = 0

Genic_SNPS_list = []

for line in Mutant_SNPS_list:
    All_SNPs+=1
    for row in Gene_Annotations_list:
        if line[0] == row[5]:
            if int(line[1]) > int(row[2]) and int(line[1]) < int(row[3]):
                    Genic_SNPs_out_open.write(line[0]+"\t"+line[1]+"\t"+line[3]+"\t"+row[0]+"\n")
                    Genic_SNPS_list.append(line)
                    Genic_SNPs+=1

SNPS_Table_out_open.write("Genic/Intergenic Characterization-"+'\n')
SNPS_Table_out_open.write("Total SNPS:"+" "+str(All_SNPs)+'\n')            
SNPS_Table_out_open.write("Genic SNPS:"+" "+str(Genic_SNPs)+'\n')

for line in Genic_SNPS_list:
    for row in Intergenic_SNPS_list:
        if line == row:
            Intergenic_SNPS_list.remove(line)

SNPS_Table_out_open.write("Intergenic SNPS:"+" "+str(len(Intergenic_SNPS_list))+'\n')


for line in Intergenic_SNPS_list:
    for row in reversed(Gene_Annotations_list):
        if line[0] == row[5]:
            if int(line[1]) > int(row[2]) and int(line[1]) > int(row[3]):
                Intergenic_SNPs_out_open.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+row[0]+"\n")
                Intergenic_SNPs_below+=1
                break
    for row in (Gene_Annotations_list):
        if line[0] == row[5]:
            if int(line[1]) < int(row[2]) and int(line[1]) < int(row[3]):
                Intergenic_SNPs_out_open.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+row[0]+"\n")
                Intergenic_SNPs_above+=1
                break
            
## The length of the intergenic SNPs list doesn't correspond perfectly.

SNPS_Table_out_open.write("Intergenic SNPS Mapped from Beginning to End:"+" "+str(Intergenic_SNPs_below)+'\n')
SNPS_Table_out_open.write("Intergenic SNPS Mapped from End to Beginning:"+" "+str(Intergenic_SNPs_above)+'\n')
SNPS_Table_out_open.write('\n'+'\n')
    
Genic_SNPs_out_open.close()
Intergenic_SNPs_out_open.close()
SNPS_Table_out_open.close()
