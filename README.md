# PG-Pipeline

## IN ACTIVE DEVELOPMENT

Nothing here is particularly polished or complete — if you're looking for a working demo, this probably isn't the branch you're looking for!

## Commenting Convention

These "codetags" help to quickly grep for different categories of information throughout the codebase that's useful to keep track of.

- `FIXME` — Issues that need addressing before the code can be considered stable or complete
- `TODO` — Features that can be added, but that aren't essential for completeness or stability
- `SAFETY` — Lays out the reasoning for why code that could theoretically panic, never will (e.g. justifying the use of `.unwrap()`)
- `PERF` — Flags up where future performance optimizations could be made
- `DESIGN` — Explains why things were done the way they were: which alternatives were discarded and why, and which downsides were explicitly accepted
- `MISSING` — Explains an intentional omission, or why something expected to be there isn't
- `NOTE` — Yes, I know this code is a bit unusual... Let me explain myself...
