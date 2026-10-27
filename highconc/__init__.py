"""highconc — figure code for the high-concentration odor-mapping manuscript.

Panels live in `highconc.panels.*` and each one draws into an `ax` you pass in;
composition stays visible in the notebook. Style comes from `highconc.style`,
cluster identity and the 1-7 numbering from `highconc.clusters`, and the data
from `highconc.data.load()`.

    import matplotlib.pyplot as plt
    from highconc import style as hs, clusters as cl

    hs.apply()
    fig, ax = plt.subplots(figsize=(hs.WIDTH_1COL, 2.6))
"""

from . import (clusters, data, dimsweep, flow, pls, recruitment, schematic,
               style, transfer)
from .data import Dataset, load

__all__ = ["clusters", "data", "dimsweep", "flow", "pls", "recruitment",
           "schematic", "style", "transfer", "Dataset", "load"]
__version__ = "0.1.0"
