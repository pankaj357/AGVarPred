# Frontend plan: React + Tailwind

Build the frontend as a separate app in `webserver/frontend/` in the next milestone. Suggested stack: React, TypeScript, Vite, Tailwind CSS, and a maintained table/chart library.

Initial screens:
1. Home: concise purpose, intended use, limitations and link to documentation.
2. Predict: GRCh38 coordinate form and VCF/TSV upload (upload route to be added only with safe size limits).
3. Results: variant ID, gene, calibrated AGVarPred score, operating-threshold class, allele frequency and consequence.
4. Variant detail: available functional-track summaries and annotation provenance.
5. Help: example input, input conventions, output interpretation, supported variant types and known limitations.

Scientific safeguards:
- Never label the output as a clinical diagnosis.
- Call 0.42 an empirically chosen operating threshold, not a clinical cutoff.
- Clearly disclose that calibration was fitted on the ClinVar calibration distribution and is not guaranteed to transport to different prevalence settings.
- Show per-variant errors/unsupported cases instead of silently discarding variants.
- Do not imply that a score exists until the validated live feature-provider adapter returns all required data.
