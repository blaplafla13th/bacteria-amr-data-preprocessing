bwa-mem2 index sequence.fasta &
python ../scripts/filter_data.py
python ../scripts/select_data.py
bash ../scripts/ena-downloader.sh &
watch_process=$!
sleep 100
while ps -p $watch_process > /dev/null; do nextflow ../scripts/VariantCalling.nf -resume; sleep 100; done
nextflow ../scripts/VariantCalling.nf -resume
python ../scripts/final.py