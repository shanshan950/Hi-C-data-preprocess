## The GM12878 in-situ Hi-C data is from GEO accession number GSE63525 PMID: 25497547
## restriction enzyme: Mbol
### Step1: Download data
#### Download data (fastq-dump)
#### Download snp bed file 
`wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_SNPs.txt.gz`
### reformat snp file to a bed file for processing bams
`cat GSE63525_GM12878_SNPs.txt.gz | gunzip | sed s/":"/\\t/g | awk '{print "chr"$1 "\t" $2 "\t" $2  "\t" $3"/"$4}' > GSE63525_GM12878_SNPs.reformat.bed`
#### mask hg19.fa by GSE63525_GM12878_SNPs.reformat.bed
`bedtools maskfasta -fi hg19.fa -bed GSE63525_GM12878_SNPs.reformat.bed -fo hg19.masked.fa`
#### rebuild hg19 bowtie2Index with hg19.masked.fa
`bowtie2-build hg19.masked.fa`
### Step2: run mapping with bowtie2 (hg19 build)
#### To cover more reads overlapping snps, we use the full length to do mapping (the mapping step is same as HiC-Pro(https://github.com/nservant/HiC-Pro))
#### Edit the path to bowtiepath and path to lib and fragment bed file
```
hg19=Your_hg19_Bowtie2IndexPath/YourIndexPrefix
hg19Ind=Your_PathTo_hg19.fa.fai
lib=Path_to_lib
bed=Path_to_fragbed
SRR=SRR1658587
```
```
bowtie2 --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder --rg-id BMG --rg SM:${SRR}_R1 -p 4 -x $hg19 -U ${SRR}_1.fastq | samtools view -bS - >  ${SRR}_R1_hg19.masked.bam &
bowtie2 --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder --rg-id BMG --rg SM:${SRR}_R2 -p 4 -x $hg19 -U ${SRR}_2.fastq | samtools view -bS - >  ${SRR}_R2_hg19.masked.bam &
wait
samtools sort -@ 2 -n -T ${SRR}_R1_hg19.masked -o ${SRR}_R1_hg19.masked.sorted.bam ${SRR}_R1_hg19.masked.bam &
samtools sort -@ 2 -n -T ${SRR}_R2_hg19.masked -o ${SRR}_R2_hg19.masked.sorted.bam${SRR}_R2_hg19.masked.bam &
wait
mv ${SRR}_R1_hg19.masked.sorted.bam ${SRR}_R1_hg19.masked.bam
mv ${SRR}_R2_hg19.masked.sorted.bam ${SRR}_R2_hg19.masked.bam
python $lib/mergeSAM.py -q 0 -t -v -f ${SRR}_R1_hg19.masked.bam -r ${SRR}_R2_hg19.masked.bam -o ${SRR}_hg19.masked.bwt2pairs.bam
```

```
for expt in SRR4271980 SRR4271981 SRR4271982 SRR4271983;do
  fq1=${expt}_1.fastq
  fq2=${expt}_2.fastq
  length=`head $fq1 | tail -1 | wc -m`
  let length=$length-1
  let trlen=$length-36
  bowtie -v 3 -m 1 --trim3 $trlen --best --strata --time -p 5  --sam $hg19  $fq1  $expt.R1.sam  &
  bowtie -v 3 -m 1 --trim3 $trlen --best --strata --time -p 5  --sam $hg19  $fq2  $expt.R2.sam  &
  wait
  echo Total reads count for $expt is `samtools view $expt.R1.sam | grep -vE ^@ | wc -l` >> $expt.summary.total.read_count &
  samtools view -u $expt.R1.sam | samtools sort -@ 5 -n -T $expt.R1 -o $expt.R1.sorted.bam  &
  samtools view -u $expt.R2.sam | samtools sort -@ 5 -n -T $expt.R2 -o $expt.R2.sorted.bam  &
  wait
  $lib/pairing_two_SAM_reads.pl <(samtools view $expt.R1.sorted.bam) <(samtools view $expt.R2.sorted.bam) | samtools view -bS -t $hg19Ind -o - - > $expt.bam
  echo Uniquely mapped read pairs for $expt is `samtools view $expt.bam | wc -l` >> $expt.summary.total.read_count
done &
```
### Step3: Process the bam files
#### check the bam files:
`ls *.bam | grep -v sorted` # You will have SRR4271980.bam SRR4271981.bam SRR4271982.bam SRR4271983.bam
#### These SRR files belong to the same biological replicate, therefore we merge the bam file first and then remove the duplicates
`samtools merge $outputname.bam SRR4271980.bam SRR4271981.bam SRR4271982.bam SRR4271983.bam`
#### remove duplicates
`samtools sort $outputname.bam | samtools view - | $lib/remove_dup_PE_SAM_sorted.pl | samtools view -bS -t $hg19Ind -o - - > $outputname.sorted.nodup.bam `
### Step4: map reads pair to fragment pairs, 36 is the read length for mapping
```
$lib/reads_2_cis_frag_loop.pl $bed 36 $outputname.loop.inward $outputname.loop.outward $outputname.loop.samestrand summary.frag_loop.read_count $outputname $outputname.temp &
$lib/reads_2_trans_frag_loop.pl $bed 36 $outputname.loop.trans $outputname.temp &
wait
for file in $outputname.loop.inward $outputname.loop.outward $outputname.loop.samestrand;do
        cat $file | $lib/summary_sorted_frag_loop.pl $bed  > temp.$file &
done
wait
for file in $outputname.loop.inward $outputname.loop.outward $outputname.loop.samestrand;do
        $lib/resort_by_frag_id.pl $bed temp.$file &
done
cat $outputname.loop.trans | $lib/summary_sorted_trans_frag_loop.pl - > temp.$outputname.loop.trans
wait
$lib/merge_sorted_frag_loop.pl temp.$outputname.loop.samestrand temp.$outputname.loop.inward temp.$outputname.loop.outward > frag_loop.$outputname.cis &
$lib/merge_sorted_frag_loop.pl temp.$outputname.loop.trans > frag_loop.$outputname.trans &
wait
```
### frag_loop.$outputname.cis and frag_loop.$outputname.trans will be used to run HiCorr

