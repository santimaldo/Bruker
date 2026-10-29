"""
Shared pytest configuration for the Bruker analysis test-suite.

* Puts the repository root on ``sys.path`` so tests can ``import Datos``,
  ``import kpopt_utils`` and (later) ``import nmrtools`` exactly as the
  analysis scripts do.
* Forces the non-interactive matplotlib backend, so no figure window can
  block a test run.
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)
