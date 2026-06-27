"""Allele-frequency / annotation source providers.

This module defines a small provider interface so that local gnomAD, future
online annotation services, and the no-AF fallback can be handled uniformly.
"""

from __future__ import annotations

import os
import warnings
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

from .gnomad import query_gnomad_by_variants


class AFSource(ABC):
    """Abstract base class for AF + VEP annotation providers."""

    name: str = "abstract"

    @abstractmethod
    def is_available(self) -> bool:
        """Return True if this provider can be used."""
        raise NotImplementedError

    @abstractmethod
    def query(
        self, variant_ids: list[str]
    ) -> tuple[dict[str, float], dict[str, dict[str, Any]]]:
        """Return (af_map, vep_map) for the requested variant IDs.

        Missing annotations should return empty dictionaries or zero AF.
        """
        raise NotImplementedError


class LocalGnomADSource(AFSource):
    """Query a local gnomAD VCF (bgzipped + tabixed)."""

    name = "local_gnomad"

    def __init__(self, vcf_path: str | Path | None = None):
        self.vcf_path = vcf_path or os.environ.get("GNOMAD_VCF")
        self._available = self.vcf_path is not None and Path(self.vcf_path).exists()

    def is_available(self) -> bool:
        return self._available

    def query(
        self, variant_ids: list[str]
    ) -> tuple[dict[str, float], dict[str, dict[str, Any]]]:
        if not self.is_available():
            raise RuntimeError("Local gnomAD VCF is not available")
        return query_gnomad_by_variants(variant_ids, self.vcf_path)


class OnlineAFSource(AFSource):
    """Placeholder for a future online AF/VEP annotation provider.

    When implemented, this provider will query a remote service (e.g., a
    gnomAD/VEP REST endpoint) and return AF + VEP annotations. For now it is
    always unavailable.
    """

    name = "online"

    def __init__(self, config: dict[str, Any] | None = None):
        self.config = config or {}
        self._available = self.config.get("enabled", False)

    def is_available(self) -> bool:
        return self._available

    def query(
        self, variant_ids: list[str]
    ) -> tuple[dict[str, float], dict[str, dict[str, Any]]]:
        if not self.is_available():
            raise RuntimeError("Online AF provider is not available")
        raise NotImplementedError("Online AF provider not yet implemented")


class NoAFSource(AFSource):
    """No annotation source available.

    Returns empty AF/VEP maps so the caller can fall back to a no-AF model.
    """

    name = "none"

    def is_available(self) -> bool:
        return True

    def query(
        self, variant_ids: list[str]
    ) -> tuple[dict[str, float], dict[str, dict[str, Any]]]:
        return {}, {}


def resolve_af_source(
    gnomad_vcf: str | Path | None = None,
    online_config: dict[str, Any] | None = None,
) -> AFSource:
    """Choose the best available AF/VEP provider.

    Priority:
        1. Local gnomAD VCF (explicit path or GNOMAD_VCF env var)
        2. Online provider (if enabled)
        3. No-AF fallback
    """
    local = LocalGnomADSource(gnomad_vcf)
    if local.is_available():
        return local

    online = OnlineAFSource(online_config)
    if online.is_available():
        warnings.warn(
            "Online AF/VEP provider selected. Accuracy depends on the remote service."
        )
        return online

    return NoAFSource()
