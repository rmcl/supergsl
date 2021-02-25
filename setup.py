import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="supergsl",
    version="0.0.1",
    author="Russell McLoughlin",
    author_email="russ.mcl@gmail.com",
    description="A python implementation of the Genome Specification Language (GSL) for genetic engineering.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rmcl/supergsl",
    packages=setuptools.find_packages(),
    entry_points={
        'console_scripts': [
            'sgsl = supergsl.sgsl:main',
            'sgsl-util = supergsl.sgsl_util:main',
        ],
    },
    install_requires=[
        'invoke',
        'llvmlite',
        'rply',
        'nose',
        'coverage',
        'biopython',
        'pydna',
        'sbol2',
        'mock',
        'mypy',
        'mypy-extensions',

    ],
    extras_require={
        "plugins": [
            'docker',
            'dnachisel',
            'graphviz',
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='>=3.7',
)
