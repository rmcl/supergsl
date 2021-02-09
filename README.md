# *SuperGSL*

*SuperGSL* is a python implementation of the Genome Specification Language (GSL) for genetic engineering. *SuperGSL* is very **alpha** so it's probably a bad idea to rely on it for, well, anything.

SuperGSL is inspired by the original Genome Specification Language (See: [Paper](https://pubs.acs.org/doi/abs/10.1021/acssynbio.5b00194) & [Code](https://github.com/Amyris/GslCore)) written by Erin Wilson, Darren Platt and many others at Amyris.

[![Documentation Status](https://readthedocs.org/projects/supergsl/badge/?version=latest)](https://supergsl.readthedocs.io/en/latest/?badge=latest)
![SuperGSL Tests](https://github.com/rmcl/supergsl/workflows/SuperGSL%20Tests/badge.svg)

### New Features in SuperGSL

See how SuperGSL differs from the original gsl [here](/docs/build/gsl_vs_supergsl)

### Install *SuperGSL*

SuperGSL can be installed using the python package manager.

```
pip install supergsl
```

This will add two commands into your environment- `sgsl` and `sgsl-util`. These two commands can used to invoke the compiler and utility commands respectively.


### Install & Running *SuperGSL* with Docker

*SuperGSL* comes with a Dockerfile and docker-compose file such that you can run superGSL without polluting your local environment.

Using docker you can do the following:

#### Run the tests

To run tests and mypy type validations:

```
docker-compose run supergsl test
```

#### Run the Compiler

To run the compiler (Note: This sort of doesn't work right now because you would somehow need to get the input file into the docker container, but is useful to me for testing):
```
docker-compose run supergsl sgsl <input-file>
```
