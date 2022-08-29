from common import *

x = get_file('GoldStandard_1008_final_NoJuvs.txt', '\t')

sex = {'Male': '1', "Female": '0'}
censor = {'Y': '1', 'N': '0'}
#Sample_id Subject_id Missing Gender SurvivalTime CensoringIndicator
b = open('pheno_file.txt', 'w')
b.write('Sample_id Subject_id Missing Gender SurvivalTime CensoringIndicator\n')
for i in x[1:]:
    #print(i)
    k = [i[0], i[0], '0', sex[i[1]], i[3], censor[i[4]]]
    u = ' '.join(k)+'\n'
    b.write(u)
b.close()
