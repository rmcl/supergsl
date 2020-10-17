# Getting Started with SuperGSL

Before you can start editing organisms you probably need a few things.

## Get a Part Provider

Sequence of DNA whether they are part of your target organisms genome or exogeneous sequences that you what to transform are provided to the compiler by a **Part Provider**. The part provider you utilize depends on the format your DNA sequences are in. If they are local files on your computer such as FASTA or GenBank files then you will likely want one of the "File Part" providers. If you want to retrieve parts from an online database such as SynBioHub then you are going to want a specific provider for that source. If you are a biotech company with a private repository you may need to write your own provider to gain access to your part repository.

## Select a Assembly Method
