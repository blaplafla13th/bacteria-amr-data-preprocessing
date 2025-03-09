bwa-mem2 index sequence.fasta
python ../scripts/filter_data.py
python ../scripts/select_data.py
bash ../scripts/ena-downloader.sh
nextflow ../scripts/VariantCalling.nf -resume
python ../scripts/final.py