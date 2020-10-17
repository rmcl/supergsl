# Architecture

SuperGSL is a fairly standard [multi-pass compiler](https://en.wikipedia.org/wiki/Multi-pass_compiler). The author fondly remembers his undergraduate compiler course at UC Santa Barbara and refreshed his memories on the finer points of compiler design by refering to [1].

The superGSL compiler converts the superGSL source code into an Abstract Syntax Tree (AST) and then executes a series of passes over the AST which incrementally resolve parts, perform assemblies, generate primers and a number of other tasks before ultimately generating the final DNA sequences.

The parts of the compiler that convert raw source code into a AST is called the "front-end". The front-end of the compiler first performs [lexical analysis](https://en.wikipedia.org/wiki/Lexical_analysis) or "lexing" to chunk the source code into a series of tokens. These tokens are then passed to the [Parser](https://en.wikipedia.org/wiki/Parsing) which enforces the syntax of the SuperGSL language and generates the AST.

## References

1. Keith D. Cooper and Linda Torczon. Engineering a compiler, second edition. 2003.
