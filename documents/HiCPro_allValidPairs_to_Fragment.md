### Map HiCPro.all.valid.pairs to fragment pairs
```
lib=lib/
validPair=
fragbed=
outputname=
HiCorrPath=

cat $validPair | gunzip | cut -f2-7 | $lib/reads_2_cis_frag_loop.pl $fragbed 36 $outputname.loop.inward $outputname.loop.outward $outputname.loop.samestrand summary.frag_loop.read_count $outputname -
cat $validPair | gunzip | cut -f2-7 | $lib/reads_2_trans_frag_loop.pl ../test_HiCorr/test_V2/HiCorr/ref/Arima/hg19.Arima.frag.bed 36 $outputname.loop.trans - &
for file in $outputname.loop.inward $outputname.loop.outward $outputname.loop.samestrand;do
        cat $file | $lib/summary_sorted_frag_loop.pl $bed  > temp.$file &
done
cat $outputname.loop.trans | $lib/summary_sorted_trans_frag_loop.pl - > temp.$outputname.loop.trans &
wait
for file in $outputname.loop.inward $outputname.loop.outward $outputname.loop.samestrand;do
        $lib/resort_by_frag_id.pl $bed temp.$file &
done
wait
$lib/merge_sorted_frag_loop.pl temp.$outputname.loop.samestrand <(cat temp.$outputname.loop.inward | awk '{if($4>1000)print $0}') <(cat temp.$outputname.loop.outward | awk '{if($4>5000)print $0}') > frag_loop.$outputname.cis &
$lib/merge_sorted_frag_loop.pl temp.$outputname.loop.trans > frag_loop.$outputname.trans &
wait
cat frag_loop.$outputname.cis <(cat frag_loop.$outputname.cis | awk '{print $2 "\t" $1 "\t" $3 "\t" $4}') | sed s/"frag_"//g | sort -k1,2n -k2,2n | awk '{print "frag_"$1 "\t" "frag_"$2 "\t" $3 "\t" $4}' > frag_loop.$outputname.cis.tmp &
cat frag_loop.$outputname.trans <(cat frag_loop.$outputname.trans | awk '{print $2 "\t" $1 "\t" $3 "\t" $4}') | sed s/"frag_"//g | sort -k1,2n -k2,2n | awk '{print "frag_"$1 "\t" "frag_"$2 "\t" $3 "\t" $4}' > frag_loop.$outputname.trans.tmp &
wait
mv frag_loop.$outputname.cis.tmp frag_loop.$outputname.cis
mv frag_loop.$outputname.trans.tmp frag_loop.$outputname.trans
```
### Run HiCorr on cis and trans loop
```
$HiCorrPath/HiCorr 
bash Arima.sh ../../test_HiCorr/test_V2/HiCorr/ref/Arima/ ../../test_HiCorr/test_V2/HiCorr/bin/Arima/ frag_loop.HK2663.cis frag_loop.HK2663.trans HK2663 hg19 &

```
### Run DeepLoop on HiCorr_output
```
#Check mid-range reads:
echo `cat HiCorr_output/anchor* | awk '{sum+=$3}END{print sum/2}'` # choose model by this depth
```
# go to the prediction/ under DeepLoop
`cd /PATH_to_DeepLoop/DeepLoop/prediction`
```
for i in {1..22} X Y;do
  python3 predict_chromosome.py --full_matrix_dir HiCorr_output/ \
                                --input_name anchor_2_anchor.loop.chr${i} \
                                --h5_file DeepLoop_models/CPGZ_trained/50M.h5 \
                                --out_dir HK2662/DeepLoop/ \
                                --anchor_dir DeepLoop/ref/hg19_Arima_anchor_bed/  \
                                --chromosome chr${i} --small_matrix_size 128  --step_size 128 --dummy 5
done

```
