# Test Environment

This directory is used for creating test cases and using them as benchmarks for implementations.
The important scripts are `benchmark.py` and `gen_cases.py`.

## `gen_cases.py`
This generates test cases for the TSP problem.
It has the following syntax: `python3 gen_cases.py DIRECTORY NUM_OF_CASES [MAX_N] [SEED]`.

The `directory` specifies to which directory the test cases should be stored to.
`NUM_OF_CASES` defines how many test scenarios should be generated.
`MAX_N` limits the number of maximum points for each test case in order to make smaller or bigger test sets.
`SEED` is the seed used by the PRNG to make it reproducible.

Example usage: `python3 gen_cases.py env 10 40 101`

## `benchmark.py`
This benchmarks a given executable on some test cases.
Optionally, it compares the results to the naive implementation.
Warning, this does not account for the runtime.
Usage syntax: `python3 benchmark.py FILE TESTS [IGNORE_NAIVE=0]`.

`FILE` is the executable to be run on the tests found in `TESTS`.
When `IGNORE_NAIVE` is set to some integer except for 0, the naive implementation will not be run.

Example usage: `python3 benchmark.py /path/to/tsp env`, `python3 benchmark.py /path/to/tsp env2 1`
