import sys

def gt2int(x):
    return sum(map(int, x.split('|')))

chromD = {}

for i in sys.stdin:
    i = i.strip().split('\t')
    chrom, nil, pos, ref, alt = i[:5]
    geno = list(map(gt2int, i[5:]))
    #print(len(geno))
    #print(chrom, nil, pos, ref, alt, geno, sep=' ')
    if chrom not in chromD:
        chromD[chrom] = []
    k = ['', '.', pos, ref, alt]
    k.extend(geno)
    chromD[chrom].append(k)

allx = ''
idx_all = 0
for i in sorted(chromD):
    idx_indv = 0
    b = open('indv.chrom/'+i+'.gen', 'w')
    for j in chromD[i]:
        j[0] = 'snp_'+str(idx_indv)
        p = ' '.join(map(str, j))+'\n'
        b.write(p)
        j[0] = 'snp_'+str(idx_all)
        allx+=' '.join(map(str, j))+'\n'
        idx_indv+=1
        idx_all+=1
        if not idx_all % 1000: print(j[1:3])

c = open('all.chrom/all.chrom.gen', 'w')
c.write(allx)

