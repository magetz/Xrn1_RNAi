import sys, csv, re

Mutant_SNPS_in = sys.argv[1]
Parental_SNPS_in = sys.argv[2]
Final_SNPS_out = sys.argv[3]
Final_Homozygous_SNPS_out = sys.argv[4]
Final_Heterozygous_SNPS_out = sys.argv[5]
SNPS_Table_out = sys.argv[6]

Mutant_SNPS_in_open = csv.reader(open(Mutant_SNPS_in,'U'),delimiter="\t")
Parental_SNPS_in_open = csv.reader(open(Parental_SNPS_in,'U'),delimiter="\t")
Final_SNPS_out_open = open(Final_SNPS_out,'w')
Final_Homozygous_SNPS_out_open = open(Final_Homozygous_SNPS_out,'w')
Final_Heterozygous_SNPS_out_open = open(Final_Heterozygous_SNPS_out,'w')
SNPS_Table_out_open = open(SNPS_Table_out,'w')

Parental_SNP_list = []
Mutant_SNP_list = []
Mutant_Homozygous_SNP_list = []
Mutant_Heterozygous_SNP_list = []
Homozygous_Indel_List = []
Homozygous_G_to_A_list = []
Homozygous_C_to_T_list = []
Heterozygous_Indel_List = []
Heterozygous_G_to_A_list = []
Heterozygous_C_to_T_list = []

Parental_SNPs = 0
Mutant_SNPs = 0
Mutant_Homozygous_SNPs = 0
Mutant_Heterozygous_SNPs = 0
Unique_SNPs = 0
Unique_Homozygous_SNPs = 0
Unique_Heterozygous_SNPs = 0

Homo_Indels = 0
Homo_G_to_A = 0
Homo_C_to_T = 0
Homo_Others = 0

Hetero_Indels = 0
Hetero_G_to_A = 0
Hetero_C_to_T = 0
Hetero_Others = 0

Unique_Homo_Indels = 0
Unique_Homo_G_to_A = 0
Unique_Homo_C_to_T = 0
Unique_Homo_Others = 0

Unique_Hetero_Indels = 0
Unique_Hetero_G_to_A = 0
Unique_Hetero_C_to_T = 0
Unique_Hetero_Others = 0

for row in Parental_SNPS_in_open:
     Parental_SNP_list.append(row[0]+" "+row[1]+" "+row[3]+" "+row[4])
     Parental_SNPs += 1

for row in Mutant_SNPS_in_open:
     Mutant_SNP_list.append(row[0]+" "+row[1]+" "+row[3]+" "+row[4])
     Mutant_SNPs += 1
     if re.search('1/1',row[9]):
          Mutant_Homozygous_SNP_list.append(row[0]+" "+row[1]+" "+row[3]+" "+row[4])
          Mutant_Homozygous_SNPs += 1
          if re.search('INDEL',row[7]):
               Homozygous_Indel_List.append(row[0]+" "+row[1]+" "+row[3]+" "+row[4])
               Homo_Indels+=1
          if re.search('G',row[3]) and re.search('A',row[4]) and not re.search('INDEL',row[7]):
               Homozygous_G_to_A_list.append(row[0]+" "+row[1]+" "+row[3]+" "+row[4])
               Homo_G_to_A+=1
          if re.search('C',row[3]) and re.search('T',row[4]) and not re.search('INDEL',row[7]):
               Homozygous_C_to_T_list.append(row[0]+" "+row[1]+" "+row[3]+" "+row[4])
               Homo_C_to_T+=1
     if re.search('0/1',row[9]):
          Mutant_Heterozygous_SNP_list.append(row[0]+" "+row[1]+" "+row[3]+" "+row[4])
          Mutant_Heterozygous_SNPs += 1
          if re.search('INDEL',row[7]):
               Heterozygous_Indel_List.append(row[0]+" "+row[1]+" "+row[3]+" "+row[4])
               Hetero_Indels+=1
          if re.search('G',row[3]) and re.search('A',row[4]) and not re.search('INDEL',row[7]):
               Heterozygous_G_to_A_list.append(row[0]+" "+row[1]+" "+row[3]+" "+row[4])
               Hetero_G_to_A+=1
          if re.search('C',row[3]) and re.search('T',row[4]) and not re.search('INDEL',row[7]):
               Heterozygous_C_to_T_list.append(row[0]+" "+row[1]+" "+row[3]+" "+row[4])
               Hetero_C_to_T+=1
               
for line in Mutant_SNP_list:
    if line not in Parental_SNP_list:
        Final_SNPS_out_open.write(line +"\n")
        Unique_SNPs+=1

for line in Mutant_Homozygous_SNP_list:
    if line not in Parental_SNP_list:
        Final_Homozygous_SNPS_out_open.write(line +"\n")
        Unique_Homozygous_SNPs+=1
        
for line in Mutant_Heterozygous_SNP_list:
    if line not in Parental_SNP_list:
        Final_Heterozygous_SNPS_out_open.write(line +"\n")
        Unique_Heterozygous_SNPs+=1

for line in Homozygous_Indel_List:
    if line not in Parental_SNP_list:
         Unique_Homo_Indels+=1

for line in Homozygous_G_to_A_list:
    if line not in Parental_SNP_list:
         Unique_Homo_G_to_A+=1

for line in Homozygous_C_to_T_list:
    if line not in Parental_SNP_list:
         Unique_Homo_C_to_T+=1
         
for line in Heterozygous_Indel_List:
    if line not in Parental_SNP_list:
         Unique_Hetero_Indels+=1

for line in Heterozygous_G_to_A_list:
    if line not in Parental_SNP_list:
         Unique_Hetero_G_to_A+=1

for line in Heterozygous_C_to_T_list:
    if line not in Parental_SNP_list:
         Unique_Hetero_C_to_T+=1
         
SNPS_Table_out_open.write("SNP Characterization-"+'\n')         
SNPS_Table_out_open.write("Parental SNPS:"+" "+str(Parental_SNPs)+'\n')
SNPS_Table_out_open.write("Mutant All SNPS:"+" "+str(Mutant_SNPs)+'\n')
SNPS_Table_out_open.write("Mutant All Homozygous SNPS:"+" "+str(Mutant_Homozygous_SNPs)+'\n')
SNPS_Table_out_open.write("Mutant All Heterozygous SNPS:"+" "+str(Mutant_Heterozygous_SNPs)+'\n')
SNPS_Table_out_open.write("Mutant Unique SNPS:"+" "+str(Unique_SNPs)+'\n')
SNPS_Table_out_open.write("Mutant Unique Homozygous SNPS:"+" "+str(Unique_Homozygous_SNPs)+'\n')
SNPS_Table_out_open.write("Mutant Unique Heterozygous SNPS:"+" "+str(Unique_Heterozygous_SNPs)+'\n'+'\n')

SNPS_Table_out_open.write("All Homozygous Indels:"+" "+str(Homo_Indels)+'\n')
SNPS_Table_out_open.write("All Homozygous G to A:"+" "+str(Homo_G_to_A)+'\n')
SNPS_Table_out_open.write("All Homozygous C to T:"+" "+str(Homo_C_to_T)+'\n')
SNPS_Table_out_open.write("All Homozygous Other:"+" "+str(Mutant_Homozygous_SNPs - Homo_Indels - Homo_G_to_A - Homo_C_to_T)+'\n')
SNPS_Table_out_open.write("All Homozygous Expected from EMS:"+" "+str(Homo_G_to_A + Homo_C_to_T)+'\n')
SNPS_Table_out_open.write("All Homozygous Not Expected from EMS:"+" "+str(Homo_Indels + (Mutant_Homozygous_SNPs - Homo_Indels - Homo_G_to_A - Homo_C_to_T))+'\n'+'\n')

SNPS_Table_out_open.write("All Heterozygous Indels:"+" "+str(Hetero_Indels)+'\n')
SNPS_Table_out_open.write("All Heterozygous G to A:"+" "+str(Hetero_G_to_A)+'\n')
SNPS_Table_out_open.write("All Heterozygous C to T:"+" "+str(Hetero_C_to_T)+'\n')
SNPS_Table_out_open.write("All Heterozygous Other:"+" "+str(Mutant_Heterozygous_SNPs - Hetero_Indels - Hetero_G_to_A - Hetero_C_to_T)+'\n')
SNPS_Table_out_open.write("All Heterozygous Expected from EMS:"+" "+str(Hetero_G_to_A + Hetero_C_to_T)+'\n')
SNPS_Table_out_open.write("All Heterozygous Not Expected from EMS:"+" "+str(Hetero_Indels + (Mutant_Heterozygous_SNPs - Hetero_Indels - Hetero_G_to_A - Hetero_C_to_T))+'\n'+'\n')

SNPS_Table_out_open.write("Unique Homozygous Indels:"+" "+str(Unique_Homo_Indels)+'\n')
SNPS_Table_out_open.write("Unique Homozygous G to A:"+" "+str(Unique_Homo_G_to_A)+'\n')
SNPS_Table_out_open.write("Unique Homozygous C to T:"+" "+str(Unique_Homo_C_to_T)+'\n')
SNPS_Table_out_open.write("Unique Homozygous Other:"+" "+str(Unique_Homozygous_SNPs - Unique_Homo_Indels - Unique_Homo_G_to_A - Unique_Homo_C_to_T)+'\n')
SNPS_Table_out_open.write("Unique Homozygous Expected from EMS:"+" "+str(Unique_Homo_G_to_A + Unique_Homo_C_to_T)+'\n')
SNPS_Table_out_open.write("Unique Homozygous Not Expected from EMS:"+" "+str(Unique_Homo_Indels + (Unique_Homozygous_SNPs - Unique_Homo_Indels - Unique_Homo_G_to_A - Unique_Homo_C_to_T))+'\n'+'\n')

SNPS_Table_out_open.write("Unique Heterozygous Indels:"+" "+str(Unique_Hetero_Indels)+'\n')
SNPS_Table_out_open.write("Unique Heterozygous G to A:"+" "+str(Unique_Hetero_G_to_A)+'\n')
SNPS_Table_out_open.write("Unique Heterozygous C to T:"+" "+str(Unique_Hetero_C_to_T)+'\n')
SNPS_Table_out_open.write("Unique Heterozygous Other:"+" "+str(Unique_Heterozygous_SNPs - Unique_Hetero_Indels - Unique_Hetero_G_to_A - Unique_Hetero_C_to_T)+'\n')
SNPS_Table_out_open.write("Unique Heterozygous Expected from EMS:"+" "+str(Unique_Hetero_G_to_A + Unique_Hetero_C_to_T)+'\n')
SNPS_Table_out_open.write("Unique Heterozygous Not Expected from EMS:"+" "+str(Unique_Hetero_Indels + (Unique_Heterozygous_SNPs - Unique_Hetero_Indels - Unique_Hetero_G_to_A - Unique_Hetero_C_to_T))+'\n'+'\n')
SNPS_Table_out_open.write('\n'+'\n')

Final_SNPS_out_open.close()
Final_Homozygous_SNPS_out_open.close()
Final_Heterozygous_SNPS_out_open.close()
SNPS_Table_out_open.close()
