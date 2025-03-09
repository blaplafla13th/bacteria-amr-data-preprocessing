FROM homebrew/brew:latest

# install essential packages
USER root
RUN apt update && apt upgrade -y && apt install -y --no-install-recommends  \
    busybox parallel python3-pandas python3-pulp python3-openpyxl python-is-python3 python3-pickleshare \
    && rm -rf /var/lib/apt/lists/*

# bio package with brew
USER linuxbrew
ENV HOMEBREW_NO_AUTO_UPDATE=1
RUN brew tap brewsci/bio && brew install trimmomatic bwa-mem2 samtools bcftools freebayes vcftools plink2  \
    && brew link openjdk@11

# nextflow
RUN curl -s https://get.nextflow.io | bash  \
    && chmod +x nextflow && mkdir -p /home/linuxbrew/.local/bin/ && mv nextflow /home/linuxbrew/.local/bin/

# aspera for ena-file-downloader
RUN curl -L https://ak-delivery04-mul.dhe.ibm.com/sar/CMA/OSA/0adrj/0/ibm-aspera-connect_4.1.3.93_linux.tar.gz | tar -xzv && \
    bash ibm-aspera-connect_4.1.3.93_linux.sh && rm -rf ibm-aspera-connect_4.1.3.93_linux.sh
    # This is why i used this version: https://www.biostars.org/p/9528910/
# ena-file-downloader
RUN curl -L https://github.com/enasequence/ena-ftp-downloader/releases/download/v1.1.8/ena-file-downloader-1.1.8.zip | \
    busybox unzip -p - ena-file-downloader.jar > /home/linuxbrew/.local/bin/ena-file-downloader.jar

# other script
COPY ./scripts /home/linuxbrew/scripts

USER root
ENV PATH=/home/linuxbrew/.local/bin/:$PATH
WORKDIR /home/linuxbrew/workspace
ENTRYPOINT ["/bin/bash"]