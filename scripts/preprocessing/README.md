## Content of this folder
<span style="color: red;">TODO - Revise this whole file</span>

This folder is dedicated to convert the BioC format into the pubtator input for the ExTRI2 pipeline. To do so, run:

```
./get_all_pubtators.sh
python prepare_pubtator_for_ExTRI2.py
```

It contains the following scripts:

`get_all_pubtators.sh` downloads the BioC folders through FTP, processes them using the `bioc_to_pubtator.sh` script, and saves the output into the `BioCXML.pubtator` folder.
* `bioc_to_pubtator.sh` takes a folder of BioC collections as input and outputs a folder with pubtator format. For that, it:
    * Truncates the file to contain only the title and abstract, and only those PMIDs with gene mentions.
    * Uses `bioc_to_pubtator.py` to convert the truncated BioC into PubTator format. Saves those PMIDs that raise an AssertionError.
    * Uses `get_pmid_list.sh` to obtain the problematic PMIDs directly from PubTator3 through GET.

`prepare_pubtator_for_ExTRI2.py` prepares the output of the previous script for ExTRI2. 
___

Code:
* `postprocessing.sh` <span style="color: red;">TODO. Old. Must be removed. It was the old version of BioCXML.pubtator-to_ExTRi2.sh</span>:

Files / datasets:
* `PubMed_pubtator/`: <span style="color: red;">TODO. Old. Must be removed. It was what I gave to ExTRI2.</span>:
* `NTNU.pubtator`:
* `NTNU_IDs`:

Other:
* `mutation_analysis/`: <span style="color: red;">TODO. Explan what that does. Also, create a similar one but to check that all PMIDs are present :)</span>:

Usage:

```
# Code to run it
```


___

<span style="color: grey;">last updated: 5/8/24</span>