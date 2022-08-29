set -eoux pipefail

python3 transform.pheno.py #makes pheno input

bcftools query -f "%CHROM\tNA\t%POS\t%REF\t%ALT\t[%GT\t]\n" imputed.turtle.vcf.gz  |\
   python3 conv.geno.file.py  #makes geno files


#run gwas

mono SurvivalGWAS_SV\ v1.3.2/survivalgwas-sv.exe \
  -gf=indv.chrom/NW_007359877.1.gen \
  -sf=pheno_file.txt -threads=1 -m=cox \
  -t=SurvivalTime -c=CensoringIndicator \
  -chr=1 -lstart=0 -lstop=99999 \
  -o='test.surv.gwas.out.txt'
