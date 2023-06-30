pyViKO
======
A Python tool to generate viral knockouts.

### What is Pyviko?
pyviko stands for **P**ython **v**iral **k**nock**o**uts. pyviko is a tool for designing molecular cloning protocols in complex viruses or other organisms with overlapping genes. Our manuscript has details and worked examples: [Taylor LJ, Strebel K. Pyviko: an automated Python tool to design gene knockouts in complex viruses with overlapping genes. BMC Microbiol. 2017 Jan 7;17(1):12.](https://bmcmicrobiol.biomedcentral.com/articles/10.1186/s12866-016-0920-3).

### What is an “overprinted gene”?
An overprinted gene is defined as the extension of one gene's open reading frame into the reading frame of a second gene. A single DNA sequence can code for multiple proteins in different reading frames or by reading in different directions. For more information, see the [Wikipedia article on reading frames](https://en.wikipedia.org/wiki/Reading_frame) or this (open access) [paper on origins of overprinted genes](http://www.ncbi.nlm.nih.gov/pubmed/22821011).

### How do I install Pyviko?
Using the python package manager:

    pip install pyviko

Otherwise, you can install it directly using `setup.py`:

    python `setup.py` install
    
([What is setup.py?](http://stackoverflow.com/questions/1471994/what-is-setup-py))

### Web interface
The basic workflow is [available as a web-based JavaScript user interface](http://louiejtaylor.github.io/pyviko/). Also check out the [Quick-start guide](http://louiejtaylor.github.io/pyviko/doc/Pyviko_quick-start.pdf) for more information on using the web interface.

### How do I use Pyviko?
Here's a simple example:

    >>> from pyviko import mutation
    >>> m = mutation.Mutant(        "ATGCATCCCTCAAGTGACTAA")
    >>> m.set_over_gene(over_seq = "ATGTATGCATCCCTCAAGTGA")
    >>> m.find_mutants()
    [(0, 'ACG'), (3, 'TAA'), (3, 'TGA')]
    
There are more sample scripts in the `scripts` folder. Also check out the [Pyviko documentation](http://louiejtaylor.github.io/pyviko/doc).
