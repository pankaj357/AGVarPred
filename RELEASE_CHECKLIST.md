# AGVarPred v1.0.0 Release Checklist

This file lists the remaining steps to complete the public GitHub repository
and Zenodo release. Items marked with `[x]` are complete; `[ ]` still require
action.

## Repository ownership and URLs

- [x] Replace `<GITHUB_ORG>` with the GitHub user or organization name
      (`pankaj357`).
- [x] Decide the final repository name (`AGVarPred`).
- [ ] Create the public GitHub repository at:
      `https://github.com/pankaj357/AGVarPred`
- [x] Update all URLs in:
  - `README.md`
  - `CITATION.cff`
  - `AGVarPred-training/CITATION.cff`
  - `AGVarPred-zenodo/CITATION.cff`
  - `pyproject.toml`
  - `AGVarPred-zenodo/code_reference.txt`

## Authors and citation metadata

- [x] Replace `[GIVEN_NAME]`, `[FAMILY_NAME]`, `[INSTITUTIONAL_AFFILIATION]`,
      and the ORCID placeholder in all `CITATION.cff` files with the actual
      author list.
- [x] Add ORCID iDs for every author.
- [x] Add institutional affiliations.
- [x] Add contact emails to:
  - `README.md`
  - `SECURITY.md`
  - `CODE_OF_CONDUCT.md`
  - `pyproject.toml`

## Zenodo and DOI

- [x] Initialize git and tag `v1.0.0` locally.
- [ ] Create a GitHub Release tagged `v1.0.0` on the public repository.
- [ ] Link the GitHub repository to Zenodo and create the v1.0.0 Zenodo record.
- [ ] Replace `10.5281/zenodo.TODO` in all `CITATION.cff` files, `README.md`,
      and `MODEL_CARD.md` with the real DOI.
- [x] Update `AGVarPred-zenodo/code_reference.txt` with the real GitHub Release URL.

## Code and model metadata

- [x] Make the first git commit and tag it `v1.0.0`.
- [x] Replace the `git_commit` placeholder in:
  - `model/model_full/manifest.yaml`
  - `model/model_no_af/manifest.yaml`
  - `AGVarPred-zenodo/model/model_full/manifest.yaml`
  - `AGVarPred-zenodo/model/model_no_af/manifest.yaml`
  with the actual commit SHA (`52a8ede5e1c7b437d7bdf760016325a35f9848f7`).
  Note: this SHA refers to the initial release commit that contains the model
  artifacts. The `v1.0.0` tag points to the final metadata commit.
- [x] Regenerate `AGVarPred-zenodo/checksums.sha256` after all edits.

## AlphaGenome access

- [ ] Confirm how users will obtain the AlphaGenome SDK and API key.
- [ ] Add the official AlphaGenome SDK installation instructions to
      `docs/installation.md` when they become publicly available.

## Final checks

- [x] Run `pytest` in a clean environment and confirm all tests pass.
- [x] Run `python -m build` and confirm the wheel/sdist build succeeds.
- [x] Run the CLI on `examples/sample.vcf` both with and without `GNOMAD_VCF`.
- [x] Verify the Zenodo archive checksums with `sha256sum -c checksums.sha256`.
- [x] Choose an open-source license (Apache-2.0).
- [ ] Add repository topics/badges on GitHub once the repo is public.

## Acknowledgements

- [x] Add acknowledgements to `README.md` and `MODEL_CARD.md`.
