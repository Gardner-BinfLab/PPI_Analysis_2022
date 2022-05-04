### LazyPair standalone
- To install requirements, run ```pip install -r requirements.txt```.
- Download the classifiers at [Zenodo](https://doi.org/10.5281/zenodo.6071630) and place it in ```clfs``` or other path required for the option -c.
- Run ```python rf.py --help```.
- For example, ```python rf.py -a ../test/test_file1.fasta -p m -c ../clfs/lazypair_clfs.pickle | awk '!seen[$1,$2]++' > output.txt```. Redundant predictions could be subsequently removed using awk to avoid out of memory error for a large number of sequences.
