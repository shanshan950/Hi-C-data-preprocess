```
# specifiy the paths below
#################################################
lib=./lib
bed=../test_HiCorr/test_V2/HiCorr/ref/Arima/hg19
allValidPairs=<your allValidPairs.gz file>
name=<output name>
#################################################

cat $allValidPairs | gunzip | cut -f2-7 | $lib/reads_2_trans_frag_loop.pl $bed 36 $name.loop.trans - &
cat $allValidPairs | gunzip | cut -f2-7 | $lib/reads_2_cis_frag_loop.pl $bed 36 $name.loop.inward $name.loop.outward $name.loop.samestrand summary.frag_loop.read_count $name - &
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

$lib/merge_sorted_frag_loop.pl temp.$name.loop.samestrand <(cat temp.$name.loop.inward | awk '{if($4>1000)print $0}') <(cat temp.$name.loop.outward | awk '{if($4>5000)print $0}') > frag_loop.$name.cis &
$lib/merge_sorted_frag_loop.pl temp.$name.loop.trans > frag_loop.$name.trans &
wait

# Run HiCorr
./HiCorr DPNII frag_loop.$name.cis frag_loop.$name.trans $name <reference_genome> [options]
