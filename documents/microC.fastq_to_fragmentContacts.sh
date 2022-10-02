#/bin/bash
#############################################################################
fq1=
fq2=
expt
hg19=/hg19_bowtie2Index/hg19.fa.fai
##############################################################################

cat ${fq1} | gunzip | bowtie -v 3 -m 1 --best --strata --time -p 15  --sam $hg19  -  $name.R1.sam  &
cat ${fq2} | gunzip | bowtie -v 3 -m 1 --best --strata --time -p 15  --sam $hg19  -  $name.R2.sam  &
wait
echo Total reads count for $name is `samtools view $name.R1.sam | grep -vE ^@ | wc -l` >> summary.total.read_count &
samtools view -u $name.R1.sam | samtools sort -@ 12 -n -T $name.R1 -o $name.R1.sorted.bam  &
samtools view -u $name.R2.sam | samtools sort -@ 12 -n -T $name.R2 -o $name.R2.sorted.bam  &
wait
$lib/pairing_two_SAM_reads.pl <(samtools view $name.R1.sorted.bam) <(samtools view $name.R2.sorted.bam) | samtools view -bS -t $hg19 -o - - > $name.bam &
samtools view $name.sorted.bam | $lib/remove_dup_PE_SAM_sorted.pl | samtools view -bS -t $hg19 -o - - > $name.sorted.nodup.bam
echo Total non-duplicated read pairs for $name is `samtools view $name.sorted.nodup.bam | wc -l` >> summary.total.read_count &

samtools view $name.sorted.nodup.bam | cut -f2-8 | /lab/solexa_weissman/cweng/PlayRoom/Shanshan/lib/bam_to_temp_HiC.pl > $name.temp
$lib/reads_2_cis_frag_loop.pl $bed 50 $name.loop.inward $name.loop.outward $name.loop.samestrand summary.frag_loop.read_count $expt $name.temp &
$lib/reads_2_trans_frag_loop.pl $bed 50 $name.loop.trans $name.temp &
wait
for file in $name.loop.inward $name.loop.outward $name.loop.samestrand;do
        cat $file | $lib/summary_sorted_frag_loop.pl $bed  > temp.$file &
done
cat $name.loop.trans | $lib/summary_sorted_trans_frag_loop.pl - > temp.$name.loop.trans &
wait
for file in $name.loop.inward $name.loop.outward $name.loop.samestrand;do
        $lib/resort_by_frag_id.pl $bed temp.$file &
done
wait
cat temp.$name.loop.inward | awk '{if($4>25000)print $0}' | $lib/merge_sorted_frag_loop.pl - > frag_loop.$name.inward &
cat temp.$name.loop.outward | awk '{if($4>5000)print $0}' | $lib/merge_sorted_frag_loop.pl - > frag_loop.$name.outward &
cat temp.$name.loop.samestrand | $lib/merge_sorted_frag_loop.pl - > frag_loop.$name.samestrand &
cat temp.$name.loop.trans | $lib/merge_sorted_frag_loop.pl - > frag_loop.$name.trans &
wait
$lib/merge_sorted_frag_loop.pl frag_loop.$name.inward frag_loop.$name.outward frag_loop.$name.samestrand > frag_loop.$name.cis

echo $expt "trans:" `cat frag_loop.$name.trans | awk '{sum+=$3}END{OFMT="%f";print sum/2}'` "cis:" `cat frag_loop.$name.cis | awk '{sum+=$3}END{OFMT="%f";print sum/2}'`  >> reads.summary 

