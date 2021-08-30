# Hi-C data preprocessing
This repository introduces the **step-by-step** analytical pipleine for Hi-C data from raw fastq files to [HiCorr](https://github.com/JinLabBioinfo/HiCorr) bias correction and [DeepLoop](https://github.com/JinLabBioinfo/DeepLoop) noise-reduction or signal-enhance contacts at 5kb-10kb resolution. 
>Please install HiCorr and DeepLoop in advance.
>
>The example documents include mapping step, bowtie or bowtie2 index for hg19/mm10 are needed.
## In the documents folder, you will find examples for processing:
#### [Hi-C (restriction enzyme HindIII) in human Adrenal tissue](https://github.com/shanshan950/Hi-C-data-preprocess/blob/master/documents/Fastq-to-FragmentContact.Tissue_example.md)
[GSM2322539](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2322539)
#### [in-situ Hi-C (restriction enzyme Mbol) in human GM12878 cell line, allele-resolved](https://github.com/shanshan950/Hi-C-data-preprocess/blob/master/documents/Fastq-to-FragmentContact.Allele_resolved_Hi-C_example.md)
[GSE63525](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525) (PMID: 25497547)
#### [sn-m3C-seq(restriction enzyme DpNII) processed contact pairs in human Astro cell type](https://github.com/shanshan950/Hi-C-data-preprocess/blob/master/documents/Processed_readsPair-to-FragmentContact.sn-m3C-seq_example.md)
[GSE130711](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130711) (PMID: 31501549)
#### [Take HiCPro "allValidPairs" to run HiCorr and DeepLoop](https://github.com/shanshan950/Hi-C-data-preprocess/blob/master/documents/HiCPro_allValidPairs_to_Fragment.md)
