*SuperGSL* is a reimplementation of the Genome Specification Language written in python.

SuperGSL is inspired by the original Genome Specification Language (See: [Paper](https://pubs.acs.org/doi/abs/10.1021/acssynbio.5b00194) & [Code](https://github.com/Amyris/GslCore)) written by Erin Wilson, Darren Platt and many others at Amyris.

### New Features in SuperGSL

See how SuperGSL differs from the original gsl [here](/docs/build/gsl_vs_supergsl)


### Install & Running *SuperGSL*

*SuperGSL* comes with a Dockerfile and docker-compose file such that you can run superGSL without polluting your local environment. If for some reason you don't want to use Docker feel free to install the python dependencies specified in requirements.txt and then execute the ./sqsl python script.

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

### Inspiring Links

https://blog.usejournal.com/writing-your-own-programming-language-and-compiler-with-python-a468970ae6df
