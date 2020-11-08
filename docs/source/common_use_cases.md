# Common Use Cases

**Strain engineers often want to do things:**

1. Knockout a gene
2. Add an extra copy of a gene (overexpress)
5. Introduce non-native genes
4. Introduce a mutation to an gene
6. Introduce an exact nucleotide sequence into a gene


## Knockout a Gene

*Using a CRISPR gRNA and Homologous Recombination*

A user may want to remove a gene from an organism's genome. One approach is to first generate a template DNA part that connects the upstream and downstream flanking regions of target gene to be removed. Then introduce a double stranded break in the gene of interest. Allow homologous recombination to fix the cut with the provided template DNA.

Generate a DNA template that connects upstream and downstream flanking regions:
```
uHO ; dHO
```

Generate a CRISPR gRNA:
```
?????
```

## Add an extra copy of a native gene (overexpress a gene)

```
from S288C import <gene>



```


## Express non-native genes

There is nothing to say that your template need only contain DNA from the host organism. As long as you have a DNA template or the budget to synthesize some DNA you can just as easily introduce genes from other organisms.

```
from <cool organims> import <sweet gene>



```

## Introduce a mutation to an gene

## Introduce an exact nucleotide sequence into a gene
