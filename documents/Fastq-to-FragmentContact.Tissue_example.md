## The Adrenal Hi-C data is from GEO accession number [GSM2322539](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2322539)
### Step1: Download data 
```
for file in SRR4271980 SRR4271981 SRR4271982 SRR4271983;do fastq-dump --split-files $SRR &;done wait
```
### Step2: run mapping with bowtie (hg19 build)
#### check read length:
```
for file in `ls *_1.fastq`;do echo $file `cat $file | head -2 | tail -1 | wc -c`; done`` # Different read length, therefore using 36bp for mapping for a fair processing
# We chose 36bp as read length because the shortest read length for this example is 36bp.
```
#### Edit the path to bowtiepath and path to lib and fragment bed file
```
hg19=Your_hg19_BowtieIndexPath/YourIndexPrefix
hg19Ind=Your_PathTo_hg19.fa.fai
HiCorrPath=<where you put HiCorr>
lib=Path_to_lib
bed=Path_to_fragbed # <chr> <start> <end> <frag_id>
outputname=Adrenal
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
done &
```
### Step3: Process the bam files
#### check the bam files:
```
ls *.bam | grep -v sorted` # You will have SRR4271980.bam SRR4271981.bam SRR4271982.bam SRR4271983.bam
```
#### These SRR files belong to the same biological replicate, therefore we merge the bam file first and then remove the duplicates
```
samtools merge $outputname.bam SRR4271980.bam SRR4271981.bam SRR4271982.bam SRR4271983.bam
```
#### remove duplicates
```
samtools sort $outputname.bam | samtools view - | $lib/remove_dup_PE_SAM_sorted.pl | samtools view -bS -t $hg19Ind -o - - > $outputname.sorted.nodup.bam
samtools view $outputname.sorted.nodup.bam | cut -f2-8 | $lib/bam_to_temp_HiC.pl > $outputname.temp
```
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
$lib/merge_sorted_frag_loop.pl temp.$outputname.loop.samestrand <(cat temp.$outputname.loop.inward | awk '{if($4>1000)print $0}') <(cat temp.$outputname.loop.outward | awk '{if($4>25000)print $0}') > frag_loop.$outputname.cis &
$lib/merge_sorted_frag_loop.pl temp.$outputname.loop.trans > frag_loop.$outputname.trans &
wait
```
### frag_loop.$outputname.cis and frag_loop.$outputname.trans will be used to run [HiCorr](https://github.com/JinLabBioinfo/HiCorr)
### Step5: Run HiCorr
```
$HiCorrPath/HiCorr HindIII frag_loop.$outputname.cis frag_loop.$outputname.trans $outputname hg19 
```
#### When it's done, the "HiCorr_output/" will include "anchor_2_anchor.loop" files for each chromosome. "HiCorr_output/" will be the input directory for DeepLoop
#### The file format is

<table><tr><td>anchor_id_1</td><td>anchor_id_2</td> <td>observed_reads_count</td> <td>expected_reads_count</td></tr></table>

#### The ratio will be (observed_reads_count + dummy)/ (expected_reads_count + dummy), we use 5 as the default dummy.

### Step6: Check heatmaps from HiCorr
```
mkdir plots # In the directory as HiCorr_output
cd plots
$HiCorrPath/HiCorr Heatmap chr1 119565703 120357702 ../HiCorr_output/anchor_2_anchor.loop.chr1 hg19 HindIII
# This will generate 3 png plots in the current directory, "raw.matrix", "expt.matrix" and "ratio.matrix"
```
### Step7: Run DeepLoop

cat `ls HiCorr_output/* | grep -v p_val` | awk '{sum+=$3}END{print sum/2}' # check reads within 2Mb

### Step8: Check heatmaps from DeepLoop

