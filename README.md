# turtle-survival-gwas

## To run

The code runs thru the shell scripts:

`1.rm.single.site.contig.impute.sh`

`2.run.gwas.sh`

Run in the order they are numbered. They call all the other necessary functions. Must be run with 
`finalvcf_GRM_nodups_renamedsamples.recode.vcf.gz` and `GoldStandard_1008_final_NoJuvs.txt` in the working dir. You also need mono installed for Survival-GWAS and bcftools installed.
