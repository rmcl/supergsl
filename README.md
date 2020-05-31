*SuperGSL* is a toy reimplementation of the Genome Specification Language written in python.

SuperGSL is inspired by the original Genome Specification Language (See: [Paper](https://pubs.acs.org/doi/abs/10.1021/acssynbio.5b00194) & [Code](https://github.com/Amyris/GslCore)) written by Erin Wilson, Darren Platt and many others at Amyris. 

I'm refering to the original GSL as *fGSL* (original gsl is written in F#).

### New Features in SuperGSL

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
```


#### Part Providers

Use the `supergsl-config.json` file to define many part providers to import parts from local files or remote services such as [synbiohub](https://synbiohub.org/) or your companies private repository.

```
{
    "part_providers": [
        {
            "name": "genebank.S288C",
            "provider_class": "supergsl.plugins.file_part.GenbankFilePartProvider",
            "sequence_file_path": "examples/GCF_000146045.2_R64_genomic.gbff.gz"
        }, {
            "name": "tab.S288C",
            "provider_class": "supergsl.plugins.file_part.FeatureTableWithFastaPartProvider",
            "fasta_file_path": "examples/S288C/S288C.fsa",
            "feature_file_path": "examples/S288C/S288C_features.tab"
        }
    ]
}
```


#### Part Slicing and Hierarchial Parts

In SuperGSL, parts can be sliced to return subsequences of supplied genetic parts.

Classic *fGSL* Slice Syntax:

`<prefix>PARTNAME[<slice>]`
1. Use the <prefix> to determine the part type slice, for example `pGAL1`.
2. Use the <slice> to determine the sub-region of the part type for example pGAL1[1:300] will
    return the first 300 bp of the promoter region.


Hierarchical Part Syntax

```
PARTNAME[<subcomponent>]
PARTNAME[<subcomponent>]
```

Use the bracket notation to access child parts. This syntax supports infinite child components, though in practice
more than one or two is likely to rather confusing.

```
PARTNAME[<subcomponent>][<child-of-subcomponent]
```

You cannot utilize *fGSL* part prefix with hierarchical parts, but you can prepend `.promoter`, `.terminator`, etc
to access these regions based on standard GSL semantics. i.e `promoter` standards for first 500 bp by default or is
overriden by PROMOTER_LENGTH setting.

```
PARTNAME[<subcomponent>].promoter
```
  



### Inspiring Links

https://blog.usejournal.com/writing-your-own-programming-language-and-compiler-with-python-a468970ae6df
