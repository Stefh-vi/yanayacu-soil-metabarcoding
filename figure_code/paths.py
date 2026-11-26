"""Portable paths for the current public figure component."""
from pathlib import Path
import os
ROOT=Path(__file__).resolve().parents[1]
INP=ROOT/'data/figure_inputs'
OUT=Path(os.environ['YANAYACU_FIGURE_OUTPUT']).resolve()
QA=Path(os.environ['YANAYACU_FIGURE_WORK']).resolve()
TAB=INP/'tables'
