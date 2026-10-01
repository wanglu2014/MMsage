# MMSage matrix trajectory core

The R functions in upstream.R are extracted/adapted from monocle3 1.4.27,
Copyright (c) 2019 Cole Trapnell, under the MIT license (see LICENSE).
settings.R contains the upstream runtime defaults. runtime.R defines the
SingleCellExperiment subclass/accessors locally and replaces the two compiled
helpers used by the production path with R implementations.

Supported contract: in-memory matrix input, PCA with norm_method="none", UMAP,
unweighted Leiden clustering, learned principal graph, and explicit root cells.
The principal graph optimization, projection, loop closure and pseudotime
algorithms are preserved. BPCells storage, interactive root picking and weighted
clustering are outside this module's contract and are not supported.

Runtime does not load monocle3, BPCells, ggrastr, ragg or lme4. It still needs
SingleCellExperiment and numeric packages (see packages.txt). In particular uwot
may still require RcppEigen/RSpectra when installed from source. Removing the
monocle3 package does not remove all compilation requirements.

M1-M3 statistical processing, M4 parameter grid and output columns are unchanged.
Existing input quoting and no-root behavior are not changed by this extraction.
