# Rules

Hard constraints. These project rules take precedence over convenience and deadlines. Follow system and
developer instructions and explicit user instructions when they take precedence.

## 1. Never bend the code toward a result

Never write code so that a chosen answer comes out. Not a compound ranked first, not a metric
reaching a threshold, not a number matching a previous run.

**Allowed, and encouraged:** correcting a fault, optimising, adding a method, changing a parameter
for a stated methodological reason. Making Vina deterministic with a fixed seed is allowed — it
fixes *which* structure is sampled, not *what score it earns*. Seeding the hydrogen placement is
allowed for the same reason: it removes noise, it does not choose an outcome.

The line: **a change may decide what the method is, never what the answer is.** If a change would
be embarrassing to describe in the Methods section, it is the wrong change.

A test that asserts a specific score is a trap — it invites making the score match. Assert
*reproducibility* (two runs agree) and *properties*, not values.

## 2. Never leak anything personal

No absolute paths from a personal machine, no usernames, no machine names, no home directories,
no credentials, in code, comments, tests, docs, commit messages, or these notes.

**The one exception:** the author's initials and surname, for citation and authorship.

Discover paths at runtime. Where a path must be written down, use `<project root>`,
`<projects root>`, `<install prefix>`.

## 3. Never hardcode

No hardcoded paths, no hardcoded lists of files, receptors, compounds, or tools. Discover from
disk. Honour the environment variables. The software must work on a machine nobody anticipated.

Discover tool paths at runtime and honour documented environment variables. Runtime options
currently come from CLI arguments and input JSON. Tests honour TOPICS_TEST_DATA.

A path that only exists on the machine the code was written on is the specific failure this rule
exists to prevent: it does not error, it silently reports that everything is absent.

## 4. English, always

All code, comments, docstrings, identifiers, on-disk names, tests, documentation, commit messages
and these notes are written in English.

Preserve existing translations, legacy identifiers needed for compatibility, and the
author's name. Conversation with the user may remain in Spanish.

## 5. Keep the workspace clean

- Announce every new directory: where it is, what it is for, and whether it can be deleted.
- Do not leave scratch files in the project tree.
- Temporary work goes outside the repository.

## 6. Backward compatibility is mandatory

Older projects must keep opening. Add to the legacy maps and the layout resolution rather than
renaming in place.

## 7. Verify what only shows up at runtime

Changes that cannot fail at import time — interface widget keys, packaging, installer layout,
Windows paths, external tool invocation — get a test. Use focused runtime checks for these boundaries. PoliScreen's packaging, Windows-path,
process-lifecycle and reproducibility tests are reference patterns, not existing TOPICS tests.

## Scientific interpretation

- Record inputs, content hashes, parameters, seeds, software versions and failures per run.
- Distinguish repeatability within a recorded environment from cross-platform agreement.
- Invalidate caches when inputs, parameters or relevant tool versions change.
- Preserve molecular stereochemistry and complete atom/residue/chain identities.
- Report structural heuristics as heuristics. Do not label docking scores as measured affinity,
  cooperativity or degradation, or geometric proximity as proof of ubiquitination.
- State hypotheses and uncertainty explicitly. Check contradictory evidence before claiming a cause.

## Applying skills

Use relevant skills on demand. Ponytail may guide implementation simplicity, but it must not
remove explicit requirements, scientific validation, provenance or the rules above.
