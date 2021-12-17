######################################
Assembly Strategies
######################################

SuperGSL supports building parts for several different assembly protocols and allows users to define new methods of assembly based on their needs.


******************************************************************************
Built in Assembly Protocols
******************************************************************************

===================================================
Simple Fusion
===================================================

.. autoclass:: supergsl.plugins.builtin.fuse.FusionAssembler


===================================================
Seamless Assembly
===================================================

.. autoclass:: supergsl.plugins.pydna.seamless_ligation.SeamlessLigationAssembler


===================================================
Synthentic Oligo Synthesis
===================================================

.. autoclass:: supergsl.plugins.builtin.oligos.SyntheticOligoAssembler


******************************************************************************
Assembly Provider Setup
******************************************************************************

To avail yourself of the fusion assembler or seamless ligation based assembly strategies, add the following to your `supergsl-config.json`.

.. code-block:: json

        {
            ...
            "assemblers": [
                ...
                {
                    "name": "seamless-ligation",
                    "assembler_class": "supergsl.plugins.seamless_ligation.SeamlessLigationAssembler",
                    "assembler_options": {
                        "target_Tm": 32
                    }
                }
            ],
            ...
        }


******************************************************************************
Other Standard Methods i'd like to implement someday
******************************************************************************


=============================================================================
BioBricks
=============================================================================

"BioBrick parts are DNA sequences which conform to a restriction-enzyme assembly standard" [1]. Some great documentation of the BioBricks assembly method can be found in the j5 documentation `BioBrick article <https://j5.jbei.org/j5manual/pages/21.html>`_.

.. code-block:: gsl

        from igem import biobrick
        from S288C import ADH1, ERG10, HO

        biobrick {
            HO_pADH1_gERG10: uHO ; pADH1 ; gERG10[1:728] ; dHO
            HO_pTDA1_gERG10: uHO ; pTDA1 ; gERG10[1:728] ; dHO
        }


=============================================================================
Golden Gate Assembly
=============================================================================

For an excellent overview of this method see the `J5 Documentation <https://j5.jbei.org/j5manual/pages/23.html>`_

.. code-block:: gsl

        from S288C import ADH1, ERG10, HO

        assemble_golden_gate {
            HO_pADH1_gERG10: uHO ; pADH1 ; gERG10[1:728] (gERG10_trunc) ; dHO
            HO_pADH1_gERG10: uHO ; pADH1 ; gERG10[1:728] (gERG10_trunc) ; dHO
        }

******************************************************************************
Custom Assembly Protocol
******************************************************************************

Many biotechs have proprietary asssembly strategies and the infrastructure for building parts. SuperGSL supports this by allowing the creation of "custom assembly methods" which can capture the specific requirements of these proprietary processes.

.. code-block:: gsl

        from mycompany import secret-assembly-method
        from S288C import ADH1, ERG10, HO

        secret-assembly-method {
            HO_pADH1_gERG10: uHO ; pADH1 ; gERG10[1:728] ; dHO
            HO_pTDA1_gERG10: uHO ; pTDA1 ; gERG10[1:728] ; dHO
            HO_pGAL3_gERG10: uHO ; pGAL3 ; gERG10[1:728] ; dHO
            HO_pGAL7_gERG10: uHO ; pGAL7 ; gERG10[1:728] ; dHO
        }

You can implement your own Assembly strategy by subclassing `AssemblerBase`

.. autoclass:: supergsl.core.assembly.AssemblerBase
    :members: assemble


===================================================
Registering Custom Assembly Types
===================================================

Once you implement your Assembler, you should register it in a SuperGSL plugin and optionally add an assembly config section to your `supergsl-config.json` file.

.. code-block:: json

        {
            "assemblers": [
                {
                    "name": "assemble_secretmethod",
                    "provider_class": "mycompany.assembly.SecretMethodAssembler",
                    "assembly_options": {
                        "Tm": 32,
                        "max_part_len": 5000
                    }
                }
            ],
            "plugins": [
                ...
                "path.to.your.plugin"
            ]
        }


******************************************************************************
References
******************************************************************************


1. Reshma P Shetty, Drew Endy, Thomas F Knight, Jr. Engineering BioBrick vectors from BioBrick parts. J Biol Eng. 2008; 2: 5.
Published online 2008 Apr 14. doi: 10.1186/1754-1611-2-5. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2373286/.
