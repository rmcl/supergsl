#############################
Part Slicing
#############################

Part slicing is the extraction of a sub-regions of genetic parts or genomes.
The idea is that you can extract regions of DNA from multiple organism and recombine them
with other parts to create novel genes.

SuperGSL's part slicing functionality is based on the syntax of the original fGSL
compiler. I good place to learn more about how fGSL works is this `Paper <https://pubs.acs.org/doi/abs/10.1021/acssynbio.5b00194>`_.


*******************************************************************************
Use positional slice notation
*******************************************************************************

.. warning::
   SuperGSL uses zero-relative indexing!! fGSL starts indexing its sequences at 1.


=============================================================================
Approximate Positions
=============================================================================

The "~" prefix in slice notation specifies that the range is approximate and can be adjusted by the compiler to optimize primer anneal efficiency (if supported by your desired :ref:`assemblies:Assembly Strategies`).

=============================================================================
Postfixes
=============================================================================


* Prepending "S" specifies that the slice coordinate should be relative to the Start of the Open Reading Frame
* Prepending "E" the End of an open reading frame

=============================================================================
An example
=============================================================================

Let s1 equal the following sequence:

3'  TACTTATGTTTGCAAGGTTACAAGAGGCCAGTCTCTAAATGGTTCCAGAAAGCTTGTTT    5'
5'  ATGAATACAAACGTTCCAATGTTCTCCGGTCAGAGATTTACCAAGGTCTTTCGAACAAA    3'

in sGSL shell you can define this sequence via:

.. code-block:: gsl
    let s1 = /ATGAATACAAACGTTCCAATGTTCTCCGGTCAGAGATTTACCAAGGTCTTTCGAACAAA/


.. code-block:: gsl
    s1[0:3]
    -> ATG
    s1[0S:3]
    -> ATG
    s1[-5E:0E]
    -> ACAAA


**********************************************
Use slice "prefixes"
**********************************************

Many Part Providers allow you to use "prefix notation" to address sub regions of
genomic features.

======== ========================== ========= ===================================================
 prefix   function                   Example   Implementation
======== ========================== ========= ===================================================
 g        gene locus                 gADH1
 p        promoter                   pERG10    Return approx. 500bp upstream of the ORF
 t        terminator                 tERG10    Return approx. 500bp downstream of the stop codon
 u        upstream                   uHO       Return approx. 500bp upstream of the ORF
 d        downstream                 dHO       Return approx. 500bp downstream of the ORF
 o        open reading frame (ORF)   oERG10    Return the ORF through stop codon.
 f        fusible ORF                fERG10    Return the ORF with the stop codon omitted
 m        mRNA                       mERG10
======== ========================== ========= ===================================================

For example, you can slice from the S. cerevisiae (S288C) to extract genes as well
as their promoters, terminators, ORFS and homology regions.

.. code-block:: gsl
    from S288C import GAL3
    from builtin import detail

Retrieving the details of GAL3 returns the ORF sequence by default.

.. code-block:: gsl
    detail(GAL3)
    -> ATGAATACAAACGTT....

You can also retrieve the GAL3 promoter region which is defined as the 500bp immediately
upstream of the primary ORF.

.. code-block:: gsl
    pGAL3
    -> CGCTTTTACTATTA...

Or the terminator region...

.. code-block:: gsl
    tGAL3
    -> CACTAAACACCTTCT...

The exact semantics of what the prefixes mean is dependent on your part provider. Clearly,
promoter, terminator and upstream regions are underspecified terms so SuperGSL leaves
it to the part provider to give specific conext specific definitions.
