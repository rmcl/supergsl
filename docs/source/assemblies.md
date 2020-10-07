# Assembly Strategies

SuperGSL supports building parts for several different assembly protocols and allows users to define new methods of assembly based on their needs.

## Built in Assembly Protocols

### Simple Fusion

The most basic strategy is a "fusion" strategy where each part is annealed to its neighbors. The product of this strategy is likely useful for direct synthesis methods.
```

from S288C import ADH1, ERG10, HO

fuse {
    HO_pADH1_gERG10: uHO ; pADH1 ; gERG10[1:728] ; dHO
	HO_pTDA1_gERG10: uHO ; pTDA1 ; gERG10[1:728] ; dHO
}

```

### BioBricks

"BioBrick parts are DNA sequences which conform to a restriction-enzyme assembly standard" [1]. Some great documentation of the BioBricks assembly method can be found in the j5 documentation (BioBrick article)[https://j5.jbei.org/j5manual/pages/21.html]

```
from S288C import ADH1, ERG10, HO

biobrick {
    HO_pADH1_gERG10: uHO ; pADH1 ; gERG10[1:728] ; dHO
	HO_pTDA1_gERG10: uHO ; pTDA1 ; gERG10[1:728] ; dHO
}

```

### Golden Gate Assembly
*For an excellent overview of this method see the J5 Documentation (article)[https://j5.jbei.org/j5manual/pages/23.html]*

```
from S288C import ADH1, ERG10, HO

assemble_golden_gate {
    HO_pADH1_gERG10: uHO ; pADH1 ; gERG10[1:728] (gERG10_trunc) ; dHO
    HO_pADH1_gERG10: uHO ; pADH1 ; gERG10[1:728] (gERG10_trunc) ; dHO
}
```

## Custom Assembly Protocol

Many biotechs have proprietary asssembly strategies and the infrastructure for building parts. SuperGSL supports this by allowing the creation of "custom assembly methods" which can capture the specific requirements of these proprietary processes.

```
from S288C import ADH1, ERG10, HO

assemble_company_parts {
    HO_pADH1_gERG10: uHO ; pADH1 ; gERG10[1:728] ; dHO
    HO_pTDA1_gERG10: uHO ; pTDA1 ; gERG10[1:728] ; dHO
    HO_pGAL3_gERG10: uHO ; pGAL3 ; gERG10[1:728] ; dHO
    HO_pGAL7_gERG10: uHO ; pGAL7 ; gERG10[1:728] ; dHO
}
```

### Registering Custom Assembly Types
You can use the `assembly_strategies` section of `supergsl-config.json` to register custom assembly types or to override parameters of built in assembly types.

```
{
    "assembly_strategies": [
        {
            "name": "assemble_secretmethod",
            "provider_class": "mycompany.assembly.SecretMethodAssembler",
            "assembly_options": {
                "Tm": 32,
                "max_part_len": 5000
            }
        }
    ]
}
```


# References

1. Reshma P Shetty, Drew Endy, Thomas F Knight, Jr. Engineering BioBrick vectors from BioBrick parts. J Biol Eng. 2008; 2: 5.
Published online 2008 Apr 14. doi: 10.1186/1754-1611-2-5. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2373286/.
