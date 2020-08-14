# Contributing downstream analyses

This project uses the [Snakemake](https://snakemake.readthedocs.io) build system to generate reproducible results. See the [README](../README.md) for installation instructions. [Why Snakemake](https://vincebuffalo.com/blog/2020/03/04/understanding-snakemake.html)?

1. A build system offers an explicit way to represent the dependencies between data, code, and results.
2. Compared with [Make](https://www.gnu.org/software/make), Snakemake offers the flexibility of a scripting language and, as an extension of Python, is easy to read and write. 
3. Workflow tools like [bpipe](http://docs.bpipe.org) are more structured for defining how input files map on to output files. With Snakemake, the conceptual focus is more heavily on telling the system what *output* files we want, and the workflow automatically determines the dependencies necessary  