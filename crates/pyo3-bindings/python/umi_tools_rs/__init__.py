"""``umi_tools_rs`` — the Rust port of ``umi_tools``, as a Python extension.

Deliberately a **separate top-level package** from the original
``umi_tools`` rather than an ``umi_tools._native`` submodule.

``references/differential-testing-and-validation.md`` documents that
``maturin develop``-installing a mixed-layout ``umi_tools`` package into a
venv that already has the original ``pip install -e .``-ed silently shadows
it, so ``import umi_tools`` resolves to the shim and the not-yet-ported pure
Python modules disappear — destroying the oracle the whole validation
strategy depends on. Its preferred fix (point ``python-source`` at the real
package's source root) would install a non-editable *copy* of the original
package alongside the editable install, which is a different flavour of the
same hazard.

Using a distinct name avoids both: ``import umi_tools`` still resolves to the
untouched editable oracle, ``import umi_tools_rs`` resolves to the Rust port,
and both live in one interpreter — so the differential harness stays
same-process with no serialization round-trip, which is the actual
requirement. Recorded in ``05_design.md`` § Intentional deviations.
"""

from ._native import *  # noqa: F401,F403
from . import _native  # noqa: F401

__all__ = ["_native"]
