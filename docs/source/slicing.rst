#############################
Part Slicing
#############################

Part slicing is the extraction of a sub-regions of genetic parts or genomes.
The idea is that you can extract regions of DNA from multiple organism and recombine them
with other parts to create novel genes.

SuperGSL's part slicing functionality is based on the syntax of the original fGSL
compiler. I good place to learn more about how fGSL works is this `Paper <https://pubs.acs.org/doi/abs/10.1021/acssynbio.5b00194>`_.


**********************************************
Use slice "prefixes"
**********************************************


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

*******************************************************************************
Use positional slice notation
*******************************************************************************

.. warning::
   SuperGSL uses zero-relative indexing!! fGSL starts indexing its sequences at 1.


=============================================================================
Approximate Positions
=============================================================================

The "~" prefix in slice notation specifies that the range is approximate and can be adjusted by the compiler to optimize primer anneal efficiency (if supported by your desired [Assembler](assemblies)).

=============================================================================
Postfixes
=============================================================================


* Prepending "S" specifies that the slice coordinate should be relative to the Start of the Open Reading Frame
* Prepending "E" the End of an open reading frame

=============================================================================
Combine positional slice and prefix notation
=============================================================================
