#!/bin/bash
# Download the 9 BioC folders though FTP, convert them to pubtator format, 
# and save them into BioCXML.pubtator folder

OUTPUT_FOLDER=data/preprocessing/BioCXML.pubtator
mkdir -p $OUTPUT_FOLDER
for n in {8..9}; do
    echo
    echo Downloading BioCXML.${n}
    # Obtain and untar BioC folder (~90GB)
    wget https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/BioCXML.${n}.tar.gz
    tar -xzf BioCXML.${n}.tar.gz
    rm BioCXML.${n}.tar.gz

    # Convert from BioC to PubTator
    echo converting BioCXML.${n} to PubTator
    ./bioc_to_pubtator.sh output/BioCXML/ $OUTPUT_FOLDER/BioCXML.${n}.pubtator
    rm -r output/BioCXML
done
