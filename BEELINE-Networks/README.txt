This folder contains the following files:

├── human-tfs.csv
├── mouse-tfs.csv
└── Networks/
    ├── human/
    │   ├── HepG2-ChIP-seq-network.csv
    │   ├── hESC-ChIP-seq-network.csv
    │   ├── Non-specific-ChIP-seq-network.csv
    │   └── STRING-network.csv
    └── mouse/
        ├── mDC-ChIP-seq-network.csv
        ├── mESC-ChIP-seq-network.csv
        ├── mESC-lofgof-network.csv
        ├── mHSC-ChIP-seq-network.csv
        ├── Non-Specific-ChIP-seq-network.csv
        └── STRING-network.csv

File/Folder description:
1) human-tfs.csv: List of human TFs collected from various sources. The file contains the TF, followed by the data source.If there are multiple sources, they are separated by a pipe '|'. You can use this file along with the GeneOrdering.csv to generate datasets of type 500+TFs, 1000+TFs used in BEELINE.

2) mouse-tfs.csv: List of mouse TFs collected from various sources. The file contains the TF, followed by the data source.If there are multiple sources, they are separated by a pipe '|'.  You can use this file along with the GeneOrdering.csv to generate datasets of type 500+TFs, 1000+TFs used in BEELINE.

3) Networks/human/: Folder containing networks corresponding to various human cell types. The naming convention is similar to the one used in the BEELINE paper.

4) Networks/mouse/: Folder containing networks corresponding to various mouse cell types. The naming convention is similar to the one used in the BEELINE paper. 



