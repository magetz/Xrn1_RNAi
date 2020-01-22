import sys, csv, re

Scas_Mut_Gene_in = sys.argv[1]
Scas_Gene_List_in = sys.argv[2]
Scer_Gene_List_in = sys.argv[3]
Homology_List_in = sys.argv[4]
Scas_Scer_Gene_Annotation_out = sys.argv[5]
Scas_Scer_No_Ortholog_out = sys.argv[6]
SNPS_Table_out = sys.argv[7]

Scas_Mut_Gene_in_open = csv.reader(open(Scas_Mut_Gene_in,'U'),delimiter="\t")
Scas_Gene_List_in_open = csv.reader(open(Scas_Gene_List_in,'U'),delimiter="\t")
Scer_Gene_List_in_open = csv.reader(open(Scer_Gene_List_in,'U'),delimiter="\t")
Homology_List_in_open = csv.reader(open(Homology_List_in,'U'),delimiter="\t")
Scas_Scer_Gene_Annotation_out_open = open(Scas_Scer_Gene_Annotation_out,'w')
Scas_Scer_No_Ortholog_out_open = open(Scas_Scer_No_Ortholog_out,'w')
SNPS_Table_out_open = open(SNPS_Table_out,'a')

Scas_Mut_Gene_List = [] #If you do not do this, it only iterates through the list once. Property of csv reader.
for line in Scas_Mut_Gene_in_open:
    Scas_Mut_Gene_List.append(line)

Scas_Gene_List = [] #If you do not do this, it only iterates through the list once. Property of csv reader.
for line in Scas_Gene_List_in_open:
    Scas_Gene_List.append(line)

Scer_Gene_List = [] #If you do not do this, it only iterates through the list once. Property of csv reader.
for line in Scer_Gene_List_in_open:
    Scer_Gene_List.append(line)

Homology_List = []
for line in Homology_List_in_open:
    Homology_List.append(line)

Scas_Mut_Ortholog_List = []    
Scas_Mut_No_Ortholog_List = []

SNPS_Considered = 0

for line in Scas_Mut_Gene_List:
    for row in Scas_Gene_List:
        if line[0] == row[0]:
            Scas_Mut_Ortholog_List.append(row[0])
            Scas_Mut_No_Ortholog_List.append(row[0])
            SNPS_Considered+=1

SNPS_Table_out_open.write("Nonsense Ortholog Search-"+'\n')        
SNPS_Table_out_open.write("SNPs Considered:"+" "+str(SNPS_Considered)+'\n') 

Scas_Mut_Ortholog_Set = set(Scas_Mut_Ortholog_List)
Scas_Mut_No_Ortholog_Set = set(Scas_Mut_No_Ortholog_List)

Scer_Orthologs = len(Scas_Mut_Ortholog_Set)

SNPS_Table_out_open.write("Genes Considered:"+" "+str(Scer_Orthologs)+'\n')

Scas_Mut_Ortholog_List_2 = list(Scas_Mut_Ortholog_Set)
Scas_Mut_No_Ortholog_List_2 = list(Scas_Mut_No_Ortholog_Set)

Scas_Mut_Extended_Ortholog_List = []

for line in Scas_Mut_Ortholog_List_2:
    for row in Scas_Gene_List:
        if line == row[0]:
            Scas_Mut_Extended_Ortholog_List.append(row)

Scas_Ortholog_List = []
Temp_Scas_Ortholog_List = []
Match = 0

for line in Scas_Mut_Extended_Ortholog_List:
    for i in Scer_Gene_List:
        if re.search(i[0], line[8]):
            Match+=1
            Scer_Orthologs-=1
            Scas_Ortholog_List.append(line[0]+"\t"+i[0]+"\t"+i[6]+"\t"+i[8]) 
            Temp_Scas_Ortholog_List.append(line[0])

for line in Temp_Scas_Ortholog_List:
    for row in Scas_Mut_No_Ortholog_List_2:
        if row == line:
            Scas_Mut_No_Ortholog_List_2.remove(row)

for line in Scas_Mut_No_Ortholog_List_2:
    for row in Homology_List:
        if re.search(line, row[4]) and re.search("---", row[11]) and re.search("---", row[21]) or re.search(line, row[28]) and re.search("---", row[11]) and re.search("---", row[21]):
            if line == "NCAS0J02110":
                Match+=1
                Scas_Ortholog_List.append(line+"\t"+"---"+"\t"+"AGO1"+"\t"+"Argonaute Protein")
                Temp_Scas_Ortholog_List.append(line)
            elif line == "NCAS0C00230":
                Match+=1
                Scas_Ortholog_List.append(line+"\t"+"---"+"\t"+"DCR1"+"\t"+"Dicer Protein")
                Temp_Scas_Ortholog_List.append(line)
            else:
                Match+=1
                Scas_Ortholog_List.append(line+"\t"+"---"+"\t"+"---"+"\t"+"Likely part of transposable element")
                Temp_Scas_Ortholog_List.append(line)
        elif re.search(line, row[4]):
            for i in Scer_Gene_List:
                if re.search(row[11], i[0]):
                    Match+=1
                    Scas_Ortholog_List.append(line+"\t"+i[0]+"\t"+i[6]+"\t"+i[8])
                    Temp_Scas_Ortholog_List.append(line)
                elif re.search(row[21], i[0]) and not re.search(row[11], i[0]):
                    Match+=1
                    Scas_Ortholog_List.append(line+"\t"+i[0]+"\t"+i[6]+"\t"+"Paralog to"+" "+i[8])
                    Temp_Scas_Ortholog_List.append(line)                    
        elif re.search(line, row[28]):
            for i in Scer_Gene_List:
                if re.search(row[21], i[0]):
                    Match+=1
                    Scas_Ortholog_List.append(line+"\t"+i[0]+"\t"+i[6]+"\t"+i[8])
                    Temp_Scas_Ortholog_List.append(line)
                elif re.search(row[11], i[0]) and not re.search(row[21], i[0]):
                    Match+=1
                    Scas_Ortholog_List.append(line+"\t"+i[0]+"\t"+i[6]+"\t"+"Paralog to"+" "+i[8])
                    Temp_Scas_Ortholog_List.append(line)

SNPS_Table_out_open.write("Genes with S. cerevisiae orthologs (Includes paralogs):"+" "+str(Match)+'\n')

for line in Temp_Scas_Ortholog_List:
    for row in Scas_Mut_No_Ortholog_List_2:
        if row == line:
            Scas_Mut_No_Ortholog_List_2.remove(row)

SNPS_Table_out_open.write("Genes without S. cerevisiae orthologs:"+" "+str(len(Scas_Mut_No_Ortholog_List_2))+'\n')
SNPS_Table_out_open.write('\n'+'\n')

for line in Scas_Ortholog_List:
    Scas_Scer_Gene_Annotation_out_open.write(line+"\n")

for line in Scas_Mut_No_Ortholog_List_2:
    Scas_Scer_No_Ortholog_out_open.write(line+"\n")

        
Scas_Scer_Gene_Annotation_out_open.close()
Scas_Scer_No_Ortholog_out_open.close()
SNPS_Table_out_open.close()
