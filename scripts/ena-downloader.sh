download() {
  java -jar /home/linuxbrew/.local/bin/ena-file-downloader.jar --accessions=$1 --format=READS_FASTQ \
  --location=/home/linuxbrew/workspace --protocol=ASPERA --asperaLocation=/home/linuxbrew/.aspera/connect \
  --email=NONE # change this if you want to receive an email when the download is complete
  # change protocol to FTP if error occurs
  find reads_fastq/$1 -name "*.fastq.gz" -exec ln {} fastq \;
}
export -f download
mkdir -p fastq
cat selected_sra.txt | parallel --jobs 4 download
