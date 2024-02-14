# non-host-de-novo-assembly
Extract reads unmatched to the host genome and do de novo assembly.


**Python pipeline for non-host _de novo_ assembly**

Install following tools and setup in your path:

bwa: https://bio-bwa.sourceforge.net/

bedtools: https://bedtools.readthedocs.io/en/latest/content/installation.html

samtools: http://www.htslib.org/download/ 

SPADES: https://github.com/ablab/spades

**usage:**

denono.py -h 

**example:**

```
denovo.py align-pe <fqfile1> <fqfile2> REF [-t <x>] -o OUTPUT
denovo.py align-se <fqfile> REF [-t <x>] -o OUTPUT
denovo.py assembly-pe <home_dir> <fqfile1> <fqfile2> <out_dir> [REF]
denovo.py assembly-se <home_dir> <fqfile1> <out_dir> [REF]
denovo.py (-h | --help)
```
