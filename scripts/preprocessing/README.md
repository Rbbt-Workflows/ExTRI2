## Summary of the main documents in preprocessing/ folder

This folder contains the scripts to download and conver BioCXML documents of PubTator3 into filtered .pubtator documents to be used in the workflow.rb script.

```
# To get the list of NCBI Gene IDs considered as TF IDs run the cells in:
# get_NCBI_TF_IDs.ipynb

# Download the BioC folders through FTP, process them & save them in Pubtator form into the BioCXML.pubtator folder
./get_all_pubtators.sh

# Filter the BioCXML.pubtator to only contain relevant abstracts to be inputted to workflow.rb
python prepare_pubtator_for_ExTRI2.py
```