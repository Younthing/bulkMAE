# Domain Docs

How engineering skills should consume this repository's domain documentation.

## Layout

This is a single-context repository:

```
/
├── CONTEXT.md
├── docs/adr/
└── R/
```

## Before exploring, read these

- `CONTEXT.md` at the repository root.
- Relevant ADRs under `docs/adr/`.

If these files do not exist, proceed silently. Do not create them speculatively. The `/domain-modeling` skill creates them lazily when terminology or durable decisions are resolved.

## Use the glossary's vocabulary

When output names a domain concept—in an issue title, proposal, hypothesis, or test—use the term defined in `CONTEXT.md`. Do not drift to synonyms the glossary explicitly avoids.

If a needed concept is absent, reconsider whether the term belongs to the project or note the gap for `/domain-modeling`.

## Flag ADR conflicts

If proposed work contradicts an existing ADR, surface the conflict explicitly instead of silently overriding it.
