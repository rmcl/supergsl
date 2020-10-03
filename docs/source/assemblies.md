# Assembly Strategies

SuperGSL supports building parts for several different assembly strategies and allows users to define new assembly types based on their needs.

## Built in Assembly Strategies

### Simple Fusion

The most basic strategy is a "fusion" strategy where each part is annealed  
```

from S288C import ADH1, ERG10, HO

assemble_fusion:
    HO_pADH1_gERG10: uHO ; pADH1 ; gERG10[1:728] ; dHO

```

### BioBricks

```
from S288C import ADH1, ERG10, HO

assemble_biobricks:
    HO_pADH1_gERG10: uHO ; pADH1 ; gERG10[1:728] ; dHO

```

### Golden Gate Assembly
*For an excellent overview of this method see the J5 Documentation: https://j5.jbei.org/j5manual/pages/23.html*

```
from S288C import ADH1, ERG10, HO

assemble_golden_gate:
    HO_pADH1_gERG10: uHO ; pADH1 ; gERG10[1:728] (gERG10_trunc) ; dHO
```

## Custom Assembly Types

Many biotechs have proprietary asssembly strategies and the infrastructure for building parts. SuperGSL supports this by allowing the creation of "custom assembly types" which can capture the specific requiremetns of these proprietary processes.

```
from S288C import ADH1, ERG10, HO

assemble_companyparts:
    HO_pADH1_gERG10: uHO ; pADH1 ; gERG10[1:728] ; dHO
    HO_pTDA1_gERG10: uHO ; pTDA1 ; gERG10[1:728] ; dHO
    HO_pGAL3_gERG10: uHO ; pGAL3 ; gERG10[1:728] ; dHO
    HO_pGAL7_gERG10: uHO ; pGAL7 ; gERG10[1:728] ; dHO
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
