# Finite program search

`python -m discovery search grammar.json --out-dir runs/experiment` enumerates
construction programs using the same registry, runner, coordinate verifiers,
compressed verifiers, and operator budgets as `run`. Every registered operator
is available to this scheduler in every dimension where its mathematical
preconditions and explicit resource budgets permit it. There is no scheduler
rule choosing a family by dimension, cardinality, or a small/large threshold.

The input is an explicit **finite grammar**. It is a construction language and
a uniform scheduling mechanism, not a claim that the implemented operators or
the chosen alternatives cover all mathematical constructions.

The runnable `discovery/recipes/program_search.json` example declares six Weil
parameter choices plus a supplied three-dimensional signed-permutation action.
Every branch feeds the same orbit-search operator. With default budgets the
largest Weil matrix branch is rejected by an explicit matrix-entry ceiling;
the search still records and attempts the other branches. Its summary is
therefore incomplete even when all seven indices were attempted:

```bash
python -m discovery search discovery/recipes/program_search.json --out-dir runs/unitary-grammar
```

```json
{
  "schema_version": 1,
  "name": "generated unitary alternatives",
  "slots": [
    {
      "id": "representation",
      "choices": [
        {
          "op": "groups.weil",
          "params": {"parity": 1},
          "parameter_grid": {"p": [3, 5], "rank": [1, 2]}
        },
        {
          "op": "groups.from_generators",
          "params": {"generators_real": [[[0, 1], [-1, 0]]]}
        }
      ]
    },
    {
      "id": "candidate",
      "choices": [
        {
          "op": "groups.search",
          "inputs": ["representation"],
          "params": {"trials": 4, "linecap": 100},
          "parameter_grid": {"phase_orders": [[1, 2], [1, 3, 6]]}
        }
      ]
    }
  ],
  "outputs": ["candidate"]
}
```

Each slot chooses one node. Its alternatives can change the operator and input
edges. Inputs refer to earlier slots, giving every generated program an
explicit topological order. `params` are literal values. `parameter_grid` maps
parameter names to nonempty lists of alternative values; its keys must be
disjoint from `params`. An array-valued parameter is one value within the outer
list, as `phase_orders` illustrates. Different choices may use different input
edges or operators. Type compatibility and mathematical preconditions are
checked by the ordinary runner for each generated program.

Enumeration preserves slot order and choice-list order. Parameter names are
sorted; the last parameter axis and last slot vary fastest. Duplicate declared
alternatives remain separate candidates. Each candidate has a stable index,
slot-branch indices, and generated-program digest. The scheduler does not
remove inconvenient branches or infer infeasibility from numerical failure.

The fresh output directory contains `search.json` and an ordinary
`candidate-00000000.json` run receipt for each completed candidate. Existing
directories and files are never overwritten. A receipt can be checked and
compared with the existing commands:

```bash
python -m discovery verify runs/experiment/candidate-00000000.json
python -m discovery compare runs/experiment/candidate-00000000.json
```

`compare` still requires an independently exact spherical configuration and a
dated catalog entry for that dimension. Numerical orbit search stays numerical.
`best_exact_by_dimension` ranks only independently verified exact spherical
outputs by cardinality; it is a result within this run, not a worldwide record.

The scheduler limits are `--max-programs` (default 100),
`--time-limit-seconds` (60), `--per-program-seconds` (30),
`--max-grammar-choices` (100000 expanded node alternatives, counted before
enumerating complete programs), and `--max-receipt-bytes` (16777216).
Each candidate runs in an isolated interpreter and is killed on its wall-time
timeout. Total time is checked before each candidate, and the remaining time
caps that candidate's timeout. Process startup, result serialization, and
summary writes have overhead beyond these execution time limits. Existing
`--budget KEY=VALUE` options govern operator allocations and proof work.

An invalid grammar is rejected before creating the output directory. Individual
parameter errors, incompatible input kinds, operator budget failures, worker
failures, and timeouts are recorded with their candidate index; later candidates
are still attempted while resources remain. Such outcomes are **not** proofs
that a mathematical construction is impossible.

`enumeration_complete` means every declared candidate index was attempted from
zero. `status: complete` additionally means every candidate finished normally.
Any failure, timeout, skipped prefix, count cutoff, or time cutoff yields
`status: incomplete`. CLI exit code is 0 for complete, 2 for incomplete, and 1
for an invalid request or output error. `next_index` identifies the next unseen
index; it does not include failed candidates that should be retried.
`--start-index N` runs a suffix into a new directory with the same ordering.
Suffix runs remain incomplete with respect to the full declared grammar;
combining coverage from multiple receipts is a separate review task.

The summary stores the original grammar and digest, all attempted branch
indices, budgets, cutoffs, failed reasons, output statuses, and receipt hashes.
Receipts replay mathematical evidence. Their hashes detect changed contents;
they do not authenticate authorship, prove novelty, or establish unrestricted
optimality.
