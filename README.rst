# autoGDC

![logo](https://raw.githubusercontent.com/chasealanbrown/autoGDC/master/autoGDC_logo.png)

This package was created for automatic downloading and processing of Genomic Data Commons data.

## Why?
The GDC is a great resource; however the process of downloading the data and combining it into a simple format can be time consuming. This package makes this more streamlined.

Examples of a few questions that may be explored are provided:

* DNA Methylation Clock: What sites are chosen when training a DNA methylation clock model?  How do they differ when different techniques are applied to retreive these sites?
* RNA vs. DNA methylation: What is the correlation between transcription and DNAm for paired data?
* DNA methylation sequence modeling: How much information is contained within sequence information of DNA methylation with respect to transcription?
* Glioblastoma Differential Gene Expression: What is the DEG signature for glioblastoma patients?

These examples can help provide as a basis for how to use this package effectively.

## Installation
The package autoGDC is available on pypi:

```sh
pip install autoGDC
```

## Installation from source

```sh
git clone https://github.com/chasealanbrown/autoGDC
cd autoGDC
```

followed by either:

```sh
python setup.py install
```

or for development mode installation

```sh
python -m pip install -e .
```
