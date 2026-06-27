# Contributing to AGVarPred

Thank you for your interest in contributing!

## Getting started

1. Fork the repository and clone your fork.
2. Create a virtual environment and install in editable mode:
   ```bash
   python -m venv .venv
   source .venv/bin/activate
   pip install -e ".[test]"
   ```
3. Run the tests:
   ```bash
   pytest
   ```

## Reporting bugs

Please open a GitHub issue with:

- A minimal reproducible example.
- The exact command or code you ran.
- Your Python version and operating system.
- The output of `AGVarPred --version` or `pip show AGVarPred`.

## Reporting security issues

Please see [`SECURITY.md`](SECURITY.md) for responsible disclosure instructions.

## Pull requests

- Keep changes focused and well-scoped.
- Add or update tests for new functionality.
- Ensure all tests pass before requesting review.
- Follow the existing code style (PEP 8).
- For scientific changes (feature engineering, model training, evaluation), open
  an issue first so the impact on reproducibility can be discussed.

## Scientific reproducibility

Do not change feature engineering, model training logic, or evaluation without
extensive validation and discussion in an issue. The primary goal of this
repository is reproducible research.

## Code of conduct

All contributors are expected to follow the [Code of Conduct](CODE_OF_CONDUCT.md).
