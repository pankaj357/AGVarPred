"""Command-line interface for AGVarPred."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from agvarpred_core.feature_generator import FeatureGenerator

from .pipeline import list_models
from .predictor import AGVarPredAutoPredictor, AGVarPredPredictor


def _add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--model-dir",
        type=str,
        default=None,
        help="Root directory containing model versions (default: env AGVARPRED_MODEL_DIR or repo model/)",
    )


def cmd_predict(args: argparse.Namespace) -> int:
    """Run prediction on a VCF."""
    if not Path(args.input).exists():
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        return 1

    # Determine explicit model request if --model was provided
    requested_model = args.model
    if requested_model == "auto":
        requested_model = None

    print("Loading model (auto-selecting best available)...")
    auto_predictor = AGVarPredAutoPredictor(
        model_dir=args.model_dir,
        requested_model=requested_model,
        gnomad_vcf=args.gnomad_vcf,
    )
    print(
        f"Using {auto_predictor.predictor.model_name} "
        f"({auto_predictor.n_features} features, type={auto_predictor.model_type}, "
        f"af_source={auto_predictor.af_source_name})"
    )

    print(f"Generating features for {args.input} (alpha_mode={args.alpha_mode})...")
    generator = FeatureGenerator(
        af_source=auto_predictor.af_source,
        alpha_mode=args.alpha_mode,
        alpha_api_key=args.alpha_api_key,
        alpha_dir=args.alpha_dir,
    )
    features = generator.from_vcf(args.input)
    print(f"Generated full feature matrix: {features.shape}")

    predictions = auto_predictor.predict(features)
    print(f"Predictions: {len(predictions)} variants")

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(output_path, index=False)
    print(f"Saved predictions to {output_path}")
    return 0


def cmd_list_models(args: argparse.Namespace) -> int:
    """List available model versions."""
    models = list_models(args.model_dir)
    if not models:
        print("No models found.", file=sys.stderr)
        return 1
    print("Available models:")
    for m in models:
        print(f"  - {m}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="AGVarPred",
        description=(
            "AGVarPred: germline variant pathogenicity prediction. "
            "Automatically uses a local gnomAD VCF if available, otherwise falls back "
            "to the bundled no-AF model."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # predict
    predict_parser = subparsers.add_parser(
        "predict",
        help="Predict pathogenicity for variants in a VCF",
    )
    predict_parser.add_argument("input", help="Input VCF file")
    predict_parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output CSV path for predictions",
    )
    _add_common_args(predict_parser)
    predict_parser.add_argument(
        "--model",
        type=str,
        default="auto",
        help=(
            "Model to use: 'auto' (default), 'full', 'no_af', or a model directory name. "
            "In auto mode the best model is selected based on AF source availability."
        ),
    )
    predict_parser.add_argument(
        "--gnomad-vcf",
        type=str,
        default=None,
        help="Path to gnomAD exomes VCF (default: GNOMAD_VCF env var)",
    )
    predict_parser.add_argument(
        "--alpha-mode",
        choices=["auto", "sdk", "precomputed"],
        default="auto",
        help="How to obtain AlphaGenome features (default: auto)",
    )
    predict_parser.add_argument(
        "--alpha-api-key",
        type=str,
        default=None,
        help="AlphaGenome API key (default: ALPHAGENOME_API_KEY env var)",
    )
    predict_parser.add_argument(
        "--alpha-dir",
        type=str,
        default=None,
        help="Directory with precomputed AlphaGenome feature matrices",
    )
    predict_parser.set_defaults(func=cmd_predict)

    # list-models
    list_parser = subparsers.add_parser(
        "list-models",
        help="List installed model versions",
    )
    _add_common_args(list_parser)
    list_parser.set_defaults(func=cmd_list_models)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
