set -eoux pipefail

python3 transform.pheno.py #makes pheno input

#run gwas

mono SurvivalGWAS_SV\ v1.3.2/survivalgwas-sv.exe \
  -gf=imputed.turtle.vcf.gz \
  -sf=pheno_file.txt -threads=1 -m=cox \
  -t=SurvivalTime -c=CensoringIndicator \
  -chr=NC_024226.1 -lstart=0 -lstop=99999 \
  -p=onlysnp -o='test.surv.gwas.out.txt'
