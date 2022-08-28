set -eoxu pipefail

bcftools view -H finalvcf_GRM_nodups_renamedsamples.recode.vcf.gz | cut -f 1 |\
   sort | uniq -c | tr -s " " | grep "^ 1 " | cut -f 3 -d " "  > contigs.with.1.site.txt

bcftools view -H finalvcf_GRM_nodups_renamedsamples.recode.vcf.gz | cut -f 1,2 |\
  grep -f contigs.with.1.site.txt > contig.and.pos.1.site.txt

bcftools view -T ^contig.and.pos.1.site.txt finalvcf_GRM_nodups_renamedsamples.recode.vcf.gz \
  -Oz -o turtle.vcf.no.single.snp.contigs.vcf.gz

java -jar beagle.22Jul22.46e.jar gt=turtle.vcf.no.single.snp.contigs.vcf.gz  out=imputed.turtle ne=40000
