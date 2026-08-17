"""Domain knowledge ported from the original R Shiny app (legacy/app.R).

Contains the AI heuristic cell-type marker dictionary used for local,
API-free cluster annotation, and the reference-resource catalog surfaced in
the deconvolution workflow. Keeping this data in one module makes it easy to
extend without touching pipeline logic.
"""

from __future__ import annotations

# AI heuristic dictionary for automated spatial-domain annotation. Maps
# canonical marker gene sets to cell-type labels; cluster annotation is done
# by Jaccard similarity between each cluster's top markers and each entry.
# This runs entirely locally -- no external API calls, no data leaves the app.
AI_CELL_TYPE_DICT: dict[str, list[str]] = {
    "T-cells (CD8+)": ["CD8A", "CD8B", "GZMB", "PRF1", "NKG7"],
    "T-cells (CD4+)": ["CD4", "IL7R", "CD40LG", "LEF1"],
    "B-cells / Plasma": ["CD19", "MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN"],
    "Macrophages / Microglia": ["CD68", "CD163", "C1QA", "C1QB", "AIF1", "TREM2", "APOE"],
    "Astrocytes": ["GFAP", "AQP4", "SLC1A3", "ALDH1L1", "FGFR3"],
    "Oligodendrocytes": ["MBP", "MOG", "MAG", "PLP1", "OPALIN"],
    "Neurons (Excitatory)": ["SLC17A7", "CAMK2A", "NEUROD6", "SNAP25"],
    "Neurons (Inhibitory)": ["GAD1", "GAD2", "SLC32A1", "PVALB", "SST"],
    "Epithelial / Tumor": ["EPCAM", "KRT18", "KRT19", "CDH1", "MUC1"],
    "Fibroblasts / CAFs": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "FAP"],
    "Endothelial": ["PECAM1", "CDH5", "VWF", "CLDN5", "FLT1"],
}

# Reference resource catalog. These are pointers to public datasets a user can
# download and upload; the platform never silently downloads them.
REFERENCE_CATALOG: list[dict[str, str]] = [
    {"resource": "10x adult mouse brain FFPE", "species": "mouse", "modality": "Visium FFPE", "source": "10x Genomics"},
    {"resource": "10x adult mouse kidney FFPE", "species": "mouse", "modality": "Visium FFPE", "source": "10x Genomics"},
    {"resource": "10x human breast cancer FFPE", "species": "human", "modality": "Visium FFPE", "source": "10x Genomics/figshare"},
    {"resource": "Human DLPFC HumanPilot10x", "species": "human", "modality": "Visium", "source": "Globus jhpce#HumanPilot10x / spatialLIBD"},
    {"resource": "10x CytAssist mouse brain sections 1/2", "species": "mouse", "modality": "Visium CytAssist", "source": "10x Genomics"},
    {"resource": "SpatialResearch breast cancer", "species": "human", "modality": "legacy spatial", "source": "spatialresearch.org"},
    {"resource": "HDST breast cancer H&E", "species": "human", "modality": "H&E image", "source": "HDST publication files"},
    {"resource": "Mouse hippocampus snRNA-seq SCP1", "species": "mouse", "modality": "snRNA-seq", "source": "Broad SCP1"},
    {"resource": "Allen adult mouse cortex taxonomy", "species": "mouse", "modality": "scRNA-seq", "source": "Allen Institute"},
    {"resource": "Breast cancer scRNA-seq GSE176078", "species": "human", "modality": "scRNA-seq", "source": "NCBI GEO GSE176078"},
    {"resource": "Lieber spatialLIBD DLPFC markers", "species": "human", "modality": "marker genes", "source": "LieberInstitute/spatialLIBD"},
]
