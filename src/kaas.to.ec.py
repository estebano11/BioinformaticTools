import argparse
import csv
from pprint import pprint as pp
        ## Modified 15-12-2015: From KAAS_to_ECnumber.py
##  Use this program to identified EC number, EC abundance and K000 list in an annotation made with KASS. The program generates
##  a file (table) of EC numbers, EC abundance and K000s identified for each gene (genes can get more than one EC number).
##  All the K000 for each EC are listed separated by spaces. Also makes a file with a description of each gene annotation.
        ## 24-12-2015: Modified to include genes list for each EC number in output file


parser = argparse.ArgumentParser(description="PROGRAM TO IDENTIFY ECs and K000s IN A KAAS ANNOTATION.\n\
It required the ko_number-database.db in the same folder")
parser.add_argument("file1", help='file format: gene_id  KO-number (KAAS)')
args = parser.parse_args()

data_base = 'ko_number-database.db'
Dic_db = dict()
complete_lst = []


#readline database
## line1A format:[ K00008    SORD, gutB    L-iditol 2-dehydrogenase   [EC:1.1.1.14] ]

with open (data_base) as db:
    line1A = db.readline()[:-1]
    while line1A !='':
        
## Process lines of database

        line1B = [i for i in line1A.split('\t') if i != '']

        if len(line1B) == 3:   ## line1A = no EC number present
            line1B.append('[EC:-.-.-.-]')
        Dic_db[line1B[0]]=line1B
        line1A = db.readline()[:-1]
           
## Open file of genome information and KAAS annotation (gene_id  KOOOOO)

with open (args.file1) as annot:
    gene_id = annot.readline()[:-1].split('\t')[0:]
    while gene_id != ['']:
        if len(gene_id) == 2:
            gene_id[len(gene_id):len(gene_id)] = Dic_db.get(gene_id[1],['KO not assigned'])
            ##print (gene_id)                                                                     ## print on screen
            complete_lst.append(gene_id)
                            
        else:
            ##print(str(gene_id[0]),'NON-annotated')
            gene_id[len(gene_id):len(gene_id)] = ['non-annotated']
            complete_lst.append(gene_id)



        gene_id = annot.readline()[:-1].split('\t')[0:]
##pp (complete_lst)                                                                           ## print on screen


with open(args.file1+'.result.complete.txt','w') as fileOUT2:
    for i in complete_lst:
        fileOUT2.write(str(i)+'\n')


"""
def process_ec(i):
    if len(i) == 6:
        #for j in i[5][4:-1].split():
        for j in i[5].split():
            return((j),(i[1]),(i[0]))  ### add print K00 (i[0])
    else:
        pass #print('EC\t-.-.-.-')
"""
EC = []
NEC = []
for i in complete_lst:
    #N = [] #process_ec(i)
    if len(i) == 6:
        for j in i[5][4:-1].split():
            N = (j),(i[1]),(i[0])
            NEC.append(N)
    else:
        pass
#    if len(i) == 6:
#        NEC.append(N)
        
        

### process EC and K00 to get the list as : EC ....   3   K001 K002 K003 

SETEC = []
for i in NEC:
    EC.append(i[0])

SETEC = set(EC)
print SETEC
N = []  ## list of abundance
K = []  ## list of KOs
G = []  ## lsit of Genes
for i in SETEC:
    n=0
    k=''
    g=''
    for j in NEC:
        if i == j[0]:
            n+=1
            k+=' '+str(j[1])
            g+=' '+str(j[2])
    N.append(n)
    K.append(k)
    G.append(g)

print '\nNumber of genes with ECs:',len(NEC)    
print 'Number of ECs:',len(SETEC)
##print N
##print K
##print G

list_EC_K00 = zip(SETEC, N, K, G)


with open (args.file1+'.result.ECs.tab', 'wt') as file_result:   # Save in a file a column of values
        Tablewriter = csv.writer(file_result, delimiter='\t')
        Tablewriter.writerow(['ECnumber','Abundance','K0s','Genes'])
        for n in list_EC_K00:
            Tablewriter.writerow (n)

print 'Output '+args.file1+'.result.ECs.tab: Use for cytoscape node atributes\n'












		
