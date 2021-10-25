# Faster Modular Composition

Date of version: October 2021.

This code has been written by 
Vincent Neiger, Bruno Salvy, Éric Schost, Gilles Villard,
to accompany the article
[A] Faster Modular Composition.

Links to preprint: [arXiv](https://arxiv.org/abs/2110.08354) and [hal](https://hal.archives-ouvertes.fr/hal-03380258) .

# Software requirements

This code is based on SageMath; the version 9.4 (or later) is mandatory for the code to run. SageMath 9.4 can be downloaded from [SageMath's website](https://www.sagemath.org/).

# Quick first usage:

- launch SageMath on your machine from the folder containing this README;
- from SageMath's command line interface, type "%runfile examples.sage".

# Short file descriptions:

- `constants.sage` defines a few constants and flags used in other files;
- `simultaneous_modular_operations.sage` contains all algorithms in [A, Section 3];
- `matrix_of_relations.sage` contains the algorithm of [A, Sections 4 and 5];
- `change_of_basis.sage` contains the algorithm of [A, Section 6];
- `modular_composition_basecase.sage` contains the first algorithm of [A, Section 8];
- `tests.sage` is mainly for testing purpose;
- `examples.sage` shows basic usage of the provided algorithms and runs the main algorithms in verbose mode to show the main steps and objects they compute on four different examples;
- `output_example.txt` gives an example of output when running `examples.sage` .

# Trying your own examples:

The code in `examples.sage` can easily be used and modified to observe the run of Algorithm ModularCompositionBaseCase [A, Algorithm 8.1] on any given input polynomials f(x), a(x), g(y).

# Remarks

- This code is not meant to be an efficient, optimized implementation. In particular, many operations on matrices over `K[x]` (such as approximant bases, kernel bases, ..) are currently rather slow in SageMath, because at the moment this software lacks a good implementation of K[x]-matrix multiplication.
- If you find a bug in this code or have questions about it, please notify the authors (preferably via GitHub's "Issues" tracking system).
