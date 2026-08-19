"""Validation result recorder — controls results_<name>.csv schema.

Identical schema/behavior to the R->Python skill's ``_recorder.py``; see
``references/differential-testing-and-validation.md``. Kept in Python
because the differential-test harness itself is a Python script that
imports both the original Python implementation and the
``maturin develop``-installed Rust extension in the same process --
which requires the ``crates/pyo3-bindings/python/<pkg>/__init__.py``
wiring described in that file's own docstring to already be done (real
package, not the placeholder shim); see
``differential-testing-and-validation.md`` for the two-environment
fallback if that wiring isn't practical.

Usage::

    from _recorder import ResultRecorder

    rec = ResultRecorder("scoring_differential")
    rec.record("score_sequence", "abs_diff", value=1e-12, threshold=1e-9, tier=1)
    rec.log("score_sequence", "n_inputs=128")
    rec.save()
"""

import csv
from pathlib import Path

__all__ = ["ResultRecorder", "FIELDNAMES"]

# Canonical CSV columns — shared with update_port_status.py parser.
FIELDNAMES = ["tutorial", "function", "tier", "metric", "value", "threshold", "pass"]


class ResultRecorder:
    """Accumulate validation results and write them to a CSV file.

    Parameters
    ----------
    tutorial : str
        Comparison run name (e.g. ``"scoring_differential"``). Used as the
        ``tutorial`` column value and to derive the output filename
        ``results_<tutorial>.csv``.
    output_dir : str or Path or None
        Directory for the output CSV. Defaults to the directory containing
        this file (``validation/``).
    """

    def __init__(self, tutorial: str, output_dir: str | Path | None = None):
        if not tutorial:
            raise ValueError("tutorial name is required (e.g. 'scoring_differential')")
        self.tutorial = tutorial
        self._rows: list[dict[str, str]] = []
        if output_dir is None:
            output_dir = Path(__file__).resolve().parent
        self._path = Path(output_dir) / f"results_{tutorial}.csv"

    def record(
        self,
        function: str,
        metric: str,
        *,
        value: float,
        threshold: float,
        tier: int = 1,
        passed: bool | None = None,
    ) -> None:
        """Record a single validation check.

        Parameters
        ----------
        function : str
            What was tested (e.g. ``"score_sequence"``).
        metric : str
            Measurement type (e.g. ``"abs_diff"``, ``"pearson_r"``).
        value : float
            Measured value.
        threshold : float
            Pass/fail threshold (see the tier tables in
            ``differential-testing-and-validation.md``).
        tier : int
            Validation tier (1, 2, or 3).
        passed : bool or None
            Explicit pass/fail. If ``None``, auto-determined as
            ``value >= threshold`` for similarity metrics — invert the
            comparison yourself for "smaller is better" metrics like
            ``abs_diff`` before calling, or pass ``passed`` explicitly.
        """
        if passed is None:
            passed = float(value) >= float(threshold)
        self._rows.append({
            "tutorial": self.tutorial,
            "function": function,
            "tier": str(tier),
            "metric": metric,
            "value": f"{value:.6f}" if isinstance(value, float) else str(value),
            "threshold": str(threshold),
            "pass": str(passed).lower(),
        })
        mark = "PASS" if passed else "FAIL"
        print(f"  [{mark}] {function}: {metric}={value} (threshold={threshold})")

    def log(self, function: str, message: str) -> None:
        """Print informational output (not written to CSV)."""
        print(f"  [INFO] {function}: {message}")

    def save(self) -> None:
        """Write accumulated results to the CSV file."""
        self._path.parent.mkdir(parents=True, exist_ok=True)
        with open(self._path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
            writer.writeheader()
            writer.writerows(self._rows)
        n_pass = sum(1 for r in self._rows if r["pass"] == "true")
        n_total = len(self._rows)
        print(f"\nResults: {n_pass}/{n_total} passed -> {self._path}")

    @property
    def path(self) -> Path:
        """Output file path."""
        return self._path
