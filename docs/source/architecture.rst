Architecture
=============

SuperGSL is a python library which defines an "embedded language" for genetic engineering. It also defines a domain specific language (DSL) which can be used to encode genetic constructs.


Sequence Store and Sequence Entries
-----------------------------------

SuperGSL provides an in-memory "Sequence Store" which allows manipulation of sequences of nucleic and amino acids in a hierarchical format similar to that described by SBOL. Sequences can be imported from "references" such as a sequenced genomes or existing parts. These sequences can then be sliced into sub-regions or concatenated together to form novel sequences. The sequence store data structure does its best to preserve sequence annotations accross these annotations and recapitulate them in the final assembled constructs.


Type System
-----------



Domain Specific Language
------------------------


The DSL aspect is processed using concepts from [multi-pass compiler](https://en.wikipedia.org/wiki/Multi-pass_compiler) and the apply/eval magic described in Structure and Interpretation of Computer Programs. The author fondly remembers his undergraduate compiler course at UC Santa Barbara and refreshed his memories on the finer points of compiler design by referring to [1].

The superGSL compiler converts the superGSL source code into an Abstract Syntax Tree (AST) and then executes a series of passes over the AST to incrementally resolve parts, functions, and other primitives. It then translates to AST into an IR which is used to build assemblies, generate primers and a number of other tasks to ultimately generate final DNA sequences.

The parts of the compiler that convert raw source code into a AST is called the "front-end". The front-end of the compiler first performs [lexical analysis](https://en.wikipedia.org/wiki/Lexical_analysis) or "lexing" to chunk the source code into a series of tokens. These tokens are then passed to the [Parser](https://en.wikipedia.org/wiki/Parsing) which enforces the syntax of the SuperGSL language and generates the AST. From there, eval() methods in the AST classes generate the IR classes from the AST graph nodes.

.. mermaid::
    :align: center

    graph TB
        subgraph "Front End"
            A[Source] -->|Lexing & Parsing| B(AST)
        end
        subgraph "Backend"
        B --> C(AST)
        C -->|Eval| D(IR)
        D -->|Assemble| E[IRv2]
        end

Providers
----------------

Part providers
Function Providers
Assembler Providers


Inspiring Links
---------------

1. "Writing your own programming language and compiler with Python" by Marcelo Andrade. June 2018. https://blog.usejournal.com/writing-your-own-programming-language-and-compiler-with-python-a468970ae6df


References
-----------

1. Keith D. Cooper and Linda Torczon. Engineering a compiler, second edition. 2003.
