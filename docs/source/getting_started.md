# Getting Started with SuperGSL

## Ten-minutes to genetic engineering

### Install *supergsl*

```
pip install supergsl

```

### Option 1: Fungal Engineering with SLiC/CPAC or related technique



### Option 2: I just need a sequence to synthesize



### Option 3: All I want are primers!

Say all you want from SuperGSL is to use its convenient slice syntax to retrieve primers of your favorite gene.

Create your `primers-only.gsl` source file:
```

from S288C import HO, SDH1, GAL3

primers {
    gHO
    tSDH1
    pGAL3    
}
```

Then execute the compiler with the `primers` output format:

```
sgsl -f primers primers-only.gsl
```

You should see a `primers.txt` file in your current directory with the following:

```
<enter primer.txt example here>
```


## Toy examples complete! Time to get serious...

Before you can start editing organisms you probably need a few things...

## Get a Part Provider

Sequence of DNA whether they are part of your target organisms genome or exogeneous sequences that you what to transform are provided to the compiler by a **Part Provider**. The part provider you utilize depends on the format your DNA sequences are in. If they are local files on your computer such as FASTA or GenBank files then you will likely want one of the "File Part" providers. If you want to retrieve parts from an online database such as SynBioHub then you are going to want a specific provider for that source. If you are a biotech company with a private repository you may need to write your own provider to gain access to your part repository.

## Select a Assembly Method


## Choose your Output

Once you design your assemblies with SuperGSL, you are probably going to want to perform the needed reactions in lab or order parts from your vendor of choice.

### SBOL, GeneBank or Flat Files



### SynBioHub

This is not quite supported yet, but adding an output provider to SynBioHub is on the list!
