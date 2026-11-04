# Search-strength diagnostic, 4 November 2026

Both tests start the existing augmented local-search engine from an empty
independent set, with eight workers and a 30-second limit. The candidate pools
already contain the known configurations. This tests recovery within supplied
pools, not discovery of their geometry.

| Model | Known selected set | Recovered | Time |
| --- | --- | --- | --- |
| 19D sampled exact polar pool | 1792 auxiliary points | 1792 | 0.0025 seconds; stopped at target |
| 21D full line-label component | 960 word/label pairs | 922 | 30 seconds; limit reached |

Both returned sets were independently checked against their saved conflict
graphs. The 21D test shows that this heuristic struggles even to rediscover
the structured incumbent without a warm start. It does not measure what a
longer or different solver could achieve.

The original exact 19D pool solve ruled out 1793 auxiliary points. The full
21D CP-SAT solve timed out with a feasible component of size 960 and a reported
upper bound of 975, so a 961-point improvement remains unresolved. Neither
statement is an unrestricted upper bound on the kissing number.

The earlier 83-million move count measures correlated local moves, and should
not be interpreted as coverage of the configuration space. Future discovery
work should benchmark recovery, strengthen the unresolved 21D formulation,
and broaden candidate generation or allow coordinated base changes.
