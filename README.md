*SuperGSL* is a toy reimplementation of the Genome Specification Language written in python.

SuperGSL is inspired by the original Genome Specification Language (See: [Paper](https://pubs.acs.org/doi/abs/10.1021/acssynbio.5b00194) & [Code](https://github.com/Amyris/GslCore)) written by Erin Wilson, Darren Platt and many others at Amyris. 

I'm refering to the original GSL as *fGSL* (original gsl is written in F#).

### New Features in SuperGSL

#### Explicit Imports

In *fGSL* you specify the refgenome pragmas so that the compiler knows where to find genes and what their sequences are.

```
#refgenome S288C
uHO ; pADH1 ; gERG10[1:728] ; ### ; dHO
```
In contrast, in superGSL the import statement allows you to bring in genomes of organisms:

```
from S288C import ADH1, ERG10, HO

uHO ; pADH1 ; gERG10[1:728] ; ### ; dHO
```


#### Part Slicing and Hierarchial Parts

*fGSL* had the concept of part slicing. There are two types of slicing.. part region prefix and index slice: 

##### Prefix Slice
Region prefix slicing involves prefixing a part with a single letter to identify which region of the part was desired. For example, in the case of the GAL1 gene, `pGAL1` referes to the first 500 basepairs (bp) of the gene loosely corresponding to the promoter region.

##### Index-based Slice
Index-based slicing will feel familiar to anyone who has used a programming language with arrays, but particularly python or numpy/pandas. Index-based slices allow a user to specify an exact range of basepairs to include. For example, `GAL1[200:300]` will yield a 100 bp part corresponding to the nucleotides at 200-300 of the GAL1 gene.

The SBOL definition defines the concept of hierarchial genetic parts.

<find image>
  
Proposed syntax: `GAL1.ORF` or `GAL1['ORF']`
  



### Inspiring Links

https://blog.usejournal.com/writing-your-own-programming-language-and-compiler-with-python-a468970ae6df
