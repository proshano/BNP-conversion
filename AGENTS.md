# AGENTS.md

## Project priorities
- Our goal is publication‑ready work for the medical literature.
- Favor readable, clean, and robust code over clever or “fancy” solutions.
- Prefer explicit, traceable workflows: commands should be visible to the user rather than running silently in the background.
- Maintain high methodological (statistical) correctness in all analyses and reporting.

## Execution and transparency
- Run commands in the foreground and surface their outputs.
- Use `assets/` as the sole analysis output directory.
- Avoid hidden or background jobs unless explicitly requested.
- When changing analysis logic, explain the statistical rationale clearly.

## Code quality
- Write self-explanatory code with minimal but useful comments.
- Prefer clear naming and simple control flow.
- Guard against edge cases and missing data.

## Protocol
- The study protocol and methodological plan for this substudy is in `BNP NTproBNP protocol.md` and should be treated as the primary reference for analysis decisions.
