# GSL vs SuperGSL

SuperGSL is inspired by the original Genome Specification Language (See: [Paper](https://pubs.acs.org/doi/abs/10.1021/acssynbio.5b00194) & [Code](https://github.com/Amyris/GslCore)) written by Erin Wilson, Darren Platt and many others at Amyris.

I refer to the original GSL as *fGSL* because the original gsl is written in F#.

### Things I don't like about *fGSL*

1. *fGSL* is written in F#

2. *fGSL*'s architecture makes it hard to extend

### Other stuff that made this seem like a good idea

3. I wanted tweak and extend the language syntax

4. Writing a compiler seemed like a fun challenge


## New Features in SuperGSL

#### Explicit Imports

In *fGSL* you specify the refgenome pragma so that the compiler knows where to find genes and what their sequences are.

```
#refgenome S288C
uHO ; pADH1 ; gERG10[1:728] ; ### ; dHO
```
In contrast, in superGSL the import statement allows you to bring in genomes of organisms:

```
from S288C import ADH1, ERG10, HO

uHO ; pADH1 ; gERG10[1:728] ; ### ; dHO
