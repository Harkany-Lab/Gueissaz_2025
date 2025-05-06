---
title: "Metabolic genes expression analysis of Young Mice Hypothalamus with focus on PVN (intersection sets)"
authors:
  - name: Evgenii O. Tretiakov, PhD
    affiliation:
      - Department of Molecular Neurosciences, Center for Brain Research, Medical University of Vienna, Vienna A-1090, Austria
    roles:
      - Conceptualization
      - Formal analysis
      - Data curation
      - Investigation
      - Methodology
      - Software
      - Visualization
      - Resources
      - Writing
    corresponding: true
    orcid: 0000-0001-5920-2190
    email: evgenii.tretiakov@meduniwien.ac.at
  - name: Tibor Harkany, PhD
    affiliation:
      - Department of Molecular Neurosciences, Center for Brain Research, Medical University of Vienna, Vienna A-1090, Austria
      - Division of Molecular and Cellular Neuroendocrinology, Department of Neuroscience, Biomedicum 7D, Karolinska Institutet, Solna SE-17165, Sweden
    roles:
      - Conceptualization
      - Funding acquisition
      - Resources
      - Supervision
      - Project administration
      - Writing
    corresponding: true
    orcid: 0000-0002-6637-5900
    email: tibor.harkany@meduniwien.ac.at
date: "2025-05-06"
format:
  html:
    default: true
    toc: true
    df-print: paged
    code-fold: true
    fig-width: 9
    fig-height: 12
    fig-format: retina
    fig-responsive: true
    fig-dpi: 300
    embed-resources: true
  pdf:
    default: true
    colorlinks: true
    fontsize: 12pt
execute:
  keep-md: true
  echo: true
  error: false
  message: false
  warning: false
  debug: false
knitr:
  opts_chunk:
    autodep: true
    fig.align: center
    fig.retina: 2
    fig.width: 14
    fig.height: 12
bibliography: references.bib
---







## Setup parameters




::: {.cell layout-align="center"}

```{.r .cell-code}
# Load tidyverse infrastructure packages
suppressPackageStartupMessages({
  library(future)
  library(here)
  library(tidyverse)
  library(magrittr)
  library(stringr)
  library(skimr)
  library(RColorBrewer)
  library(viridis)
})


# Load packages for scRNA-seq analysis and visualisation
suppressPackageStartupMessages({
  library(UpSetR)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(ggstatsplot)
  library(Seurat)
  library(SeuratWrappers)
  library(scCustomize)
})
```
:::




### Set paths




::: {.cell layout-align="center"}

```{.r .cell-code}
src_dir <- here("code")
data_dir <- here("data")
output_dir <- here("output")
```
:::




### Load gene-sets




::: {.cell layout-align="center"}

```{.r .cell-code}
source(here(src_dir, "genes.R"))
```
:::




### Set fixed variables




::: {.cell layout-align="center"}

```{.r .cell-code}
# set seed
reseed <- 42
set.seed(seed = reseed)

# Parameters for parallel execution
n_cores <- min(parallelly::availableCores() / 2 - 1, 16)
plan("multisession", workers = n_cores)
options(
  future.globals.maxSize = 100000 * 1024^2,
  future.rng.onMisuse = "ignore"
)
plan()
```

::: {.cell-output .cell-output-stdout}

```
multisession:
- args: function (..., workers = 16, envir = parent.frame())
- tweaked: TRUE
- call: plan("multisession", workers = n_cores)
```


:::

```{.r .cell-code}
# ggplot2 theme
theme_set(ggmin::theme_powerpoint())
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
bioproject <- "PRJNA547712"
project <- "kim2020_Hypoth-dev"
cb_fpr <- 0.001
low_cutoff_gene <- 500
high_cutoff_gene <- NULL
high_cutoff_gene <- 5000
low_cutoff_umis <- NULL
low_cutoff_umis <- -Inf
high_cutoff_umis <- 25000
high_cutoff_pc_mt <- 15
high_cutoff_pc_ribo <- 20
high_cutoff_pc_hb <- 0.1
high_cutoff_doublet_score <- 0.33
high_cutoff_complexity <- 0.85
connectivity_model <- "min_tree"
k <- 10
metric <- "euclidean"
signature <- 100
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
opioid_system_genes <- c(
  # Classical opioid receptors
  "Oprd1", # Delta opioid receptor
  "Oprk1", # Kappa opioid receptor
  "Oprl1", # Nociceptin/orphanin FQ receptor
  "Oprm1", # Mu opioid receptor

  # Processing enzymes
  "Pcsk1", # Proprotein convertase 1
  "Pcsk2", # Proprotein convertase 2

  # Endogenous opioid precursors
  "Pdyn", # Prodynorphin
  "Penk", # Proenkephalin
  # "Pomc",   # Proopiomelanocortin
  "Pnoc" # Prepronociceptin
)
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
metabolic_signaling_genes <- c(
  # Receptor tyrosine kinases and ligands
  "Alk", # Anaplastic lymphoma kinase - neural development, metabolism
  "Fam150a", # ALK ligand 1/Augmentor-β - ALK receptor activator
  "Fam150b", # ALK ligand 2/Augmentor-α - ALK receptor activator

  # Melanocortin system
  "Mc3r", # Melanocortin 3 receptor - energy homeostasis, inflammation
  "Mc4r", # Melanocortin 4 receptor - appetite control, energy balance

  # Metabolic hormone receptors
  "Lepr", # Leptin receptor - energy balance, satiety
  "Insr", # Insulin receptor - glucose homeostasis
  # "Igf1r",    # Insulin-like growth factor 1 receptor - growth, development

  # Signaling adaptors/regulators
  "Lmo4", # LIM domain only 4 - transcriptional regulation, metabolism
  "Irs1", # Insulin receptor substrate 1 - insulin signaling
  "Irs4" # Insulin receptor substrate 4 - insulin/leptin signaling
)
```
:::




## Load Kim DW et al., 2020[@kim2020]




::: {.cell layout-align="center"}

```{.r .cell-code}
srt.kim <- schard::h5ad2seurat(
  here(
    data_dir,
    "kim2020_combined.h5ad"
  ),
  use.raw = TRUE
)

X_umap <- srt.kim@meta.data |>
  dplyr::select(X, Y, Z) |>
  as.matrix()
colnames(X_umap) <- c("UMAP_1", "UMAP_2", "UMAP_3")
rownames(X_umap) <- colnames(srt.kim)
srt.kim[["umap"]] <- CreateDimReducObject(
  embeddings = X_umap,
  key = "umap_",
  assay = DefaultAssay(srt.kim)
)
srt.kim$Age %<>%
  forcats::fct(
    levels = c(
      "E10",
      "E11",
      "E12",
      "E13",
      "E14",
      "E15",
      "E16",
      "E17",
      "E18",
      "P0",
      "P2",
      "P4",
      "P8",
      "P10",
      "P14",
      "P23",
      "P45"
    )
  )

Idents(srt.kim) <- "Age"
srt.kim <- Store_Palette_Seurat(
  seurat_object = srt.kim,
  palette = rev(brewer.pal(n = 11, name = "Spectral")),
  palette_name = "expr_Colour_Pal"
)
```
:::




## Load Romanov et al., 2020[@romanov2020]




::: {.cell layout-align="center"}

```{.r .cell-code}
print(srt.kim)
```

::: {.cell-output .cell-output-stdout}

```
An object of class Seurat 
27998 features across 128006 samples within 1 assay 
Active assay: RNA (27998 features, 0 variable features)
 2 layers present: counts, data
 1 dimensional reduction calculated: umap
```


:::

```{.r .cell-code}
srt.romanov.pub <- readRDS("/data/PRJNA548917/old/oldCCA_nae_srt.rds")
srt.romanov.pub <- UpdateSeuratObject(srt.romanov.pub)
Idents(srt.romanov.pub) <-
  factor(srt.romanov.pub$wtree, ordered = TRUE)

# Consistent colours and clusters names
colours_wtree <- setNames(read_lines(here(data_dir, "colours_wtree.tsv")), 1:45)

srt.romanov.pub$age <-
  Cells(srt.romanov.pub) |>
  str_split(pattern = ":", simplify = T) %>%
  .[, 1] %>%
  str_split_fixed(pattern = "_", n = 3) %>%
  .[, 3]
print(srt.romanov.pub)
```

::: {.cell-output .cell-output-stdout}

```
An object of class Seurat 
24340 features across 51199 samples within 1 assay 
Active assay: RNA (24340 features, 3500 variable features)
 3 layers present: counts, data, scale.data
 3 dimensional reductions calculated: pca, tsne, umap
```


:::

```{.r .cell-code}
glimpse(srt.romanov.pub@meta.data)
```

::: {.cell-output .cell-output-stdout}

```
Rows: 51,199
Columns: 20
$ nGene            <int> 1652, 782, 447, 1706, 1106, 894, 727, 734, 669, 617, …
$ nUMI             <dbl> 2787, 1090, 544, 2709, 1817, 1220, 995, 1036, 920, 86…
$ orig.ident       <fct> Hypothalamus, Hypothalamus, Hypothalamus, Hypothalamu…
$ res.0.2          <chr> "23", "23", "23", "23", "23", "23", "23", "23", "23",…
$ res.0.4          <chr> "34", "34", "34", "34", "34", "34", "34", "34", "34",…
$ res.0.8          <chr> "42", "42", "42", "42", "42", "42", "42", "42", "42",…
$ res.1.2          <chr> "47", "47", "47", "47", "47", "47", "47", "47", "47",…
$ res.2            <chr> "54", "54", "54", "54", "54", "54", "54", "54", "54",…
$ tree.ident       <int> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
$ pro_Inter        <chr> "41", "41", "41", "41", "41", "41", "41", "41", "41",…
$ pro_Enter        <chr> "41", "41", "41", "41", "41", "41", "41", "41", "41",…
$ tree_final       <fct> 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 1…
$ subtree          <fct> 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 4…
$ prim_walktrap    <fct> 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 3…
$ umi_per_gene     <dbl> 1.687046, 1.393862, 1.217002, 1.587925, 1.642857, 1.3…
$ log_umi_per_gene <dbl> 0.22712693, 0.14421974, 0.08529138, 0.20082998, 0.215…
$ nCount_RNA       <dbl> 2787, 1090, 544, 2709, 1817, 1220, 995, 1036, 920, 86…
$ nFeature_RNA     <int> 1652, 782, 447, 1706, 1106, 894, 727, 734, 669, 617, …
$ wtree            <fct> 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 3…
$ age              <chr> "P23", "3P2", "3P2", "P2", "P2", "P2", "P2", "P2", "P…
```


:::

```{.r .cell-code}
table(Idents(srt.romanov.pub))
```

::: {.cell-output .cell-output-stdout}

```

    1     2     3     4     5     6     7     8     9    10    11    12    13 
 2344  8146   395   402  3234   712   552   374   259   952 13727  1615   765 
   14    15    16    17    18    19    20    21    22    23    24    25    26 
  832  1244   792   590   808  2486  1683   628  1039  1750   292   394   547 
   27    28    29    30    31    32    33    34    35    36    37    38    39 
  391   407   507    93    81   402   143   701   222   353   324    73    78 
   40    41    42    43    44    45 
  328   190    73    37   179    55 
```


:::

```{.r .cell-code}
srt.romanov.pub %<>% RenameIdents(object = ., `43` = "mneOXY")
srt.romanov.pub %<>% RenameIdents(object = ., `26` = "mneVAS")
srt.romanov.pub %<>% RenameIdents(object = ., `31` = "pneSS")
srt.romanov.pub %<>% RenameIdents(object = ., `24` = "pneCRH")
srt.romanov.pub %<>% RenameIdents(object = ., `15` = "pneTRH")
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
srt.romanov.pub$stage <-
  srt.romanov.pub$age %>%
  forcats::fct_collapse(
    Embryonic = c("E15", "E17"),
    Neonatal = c("P0", "P2", "3P2"),
    Pubertal = c("1P10", "P10"),
    Adult = c("P23")
  )
srt.romanov.pub$stage %<>%
  factor(
    levels = c("Embryonic", "Neonatal", "Pubertal", "Adult"),
    ordered = TRUE
  )
srt.romanov.pub$stage %>% forcats::fct_count()
```

::: {.cell-output-display}

`````{=html}
<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["f"],"name":[1],"type":["fct"],"align":["left"]},{"label":["n"],"name":[2],"type":["int"],"align":["right"]}],"data":[{"1":"Embryonic","2":"19503"},{"1":"Neonatal","2":"20316"},{"1":"Pubertal","2":"8965"},{"1":"Adult","2":"2415"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>
`````

:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
srt.romanov.pub$age <-
  plyr::mapvalues(
    x = srt.romanov.pub$age,
    from = c("E15", "E17", "P0", "P2", "3P2", "1P10", "P10", "P23"),
    to = c("E15", "E17", "P00", "P02", "P02", "P10", "P10", "P23")
  )
srt.romanov.pub$age %>% forcats::fct_count()
```

::: {.cell-output-display}

`````{=html}
<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["f"],"name":[1],"type":["fct"],"align":["left"]},{"label":["n"],"name":[2],"type":["int"],"align":["right"]}],"data":[{"1":"E15","2":"8290"},{"1":"E17","2":"11213"},{"1":"P00","2":"7492"},{"1":"P02","2":"12824"},{"1":"P10","2":"8965"},{"1":"P23","2":"2415"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>
`````

:::
:::




## Prepare query mapping between datasets




::: {.cell layout-align="center"}

```{.r .cell-code}
srt.kim <- NormalizeData(srt.kim)
srt.kim <- FindVariableFeatures(
  srt.kim,
  selection.method = "vst",
  nfeatures = 3000
)
# all.genes <- rownames(srt.kim)
# srt.kim <- ScaleData(srt.kim, features = all.genes)
srt.kim <- ScaleData(srt.kim)
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
hypoth.anchors <- FindTransferAnchors(
  reference = srt.romanov.pub,
  query = srt.kim,
  dims = 1:30,
  reference.reduction = "pca"
)
predictions <- TransferData(
  anchorset = hypoth.anchors,
  refdata = srt.romanov.pub$wtree,
  dims = 1:30
)
srt.kim <- AddMetaData(srt.kim, metadata = predictions)
table(srt.kim$predicted.id)
```

::: {.cell-output .cell-output-stdout}

```

     1     10     11     12     13     15     16     17     19      2     20 
  3700      1   2198    225    545     58    218     19 100308   9886    924 
    22     23     26     27     28     29      3     38      4      5      6 
   163     11    115     24     76    481    951     49    156   2895   4756 
     7      8      9 
   199     47      1 
```


:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
srt.romanov.pub <- RunUMAP(
  srt.romanov.pub,
  dims = 1:30,
  reduction = "pca",
  return.model = TRUE
)
srt.kim <- IntegrateEmbeddings(
  anchorset = hypoth.anchors,
  reference = srt.romanov.pub,
  query = srt.kim,
  new.reduction.name = "ref.pca"
)
srt.kim <- ProjectUMAP(
  query = srt.kim,
  query.reduction = "ref.pca",
  reference = srt.romanov.pub,
  reference.reduction = "pca",
  reduction.model = "umap"
)
Idents(srt.kim) <- srt.kim$Cluster
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
all.genes <- rownames(srt.kim)
gene.scale <- c(
  cnbn,
  opioid_system_genes,
  metabolic_signaling_genes,
  np,
  npr,
  nmr,
  neurotrans
) |>
  unique() %>%
  .[. %in% all.genes]
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
srt.kim$stage <-
  srt.kim$Age %>%
  forcats::fct_collapse(
    Embryonic = c(
      "E10",
      "E11",
      "E12",
      "E13",
      "E14",
      "E15",
      "E16",
      "E18"
    ),
    Neonatal = c("P4", "P8"),
    Pubertal = c("P14"),
    Adult = c("P45")
  )
srt.kim$stage %<>%
  factor(
    levels = c("Embryonic", "Neonatal", "Pubertal", "Adult"),
    ordered = TRUE
  )
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
selected_cells <- read_lines(here(data_dir, "kim2020_pvn_neurons.txt"))
srt.kim <- subset(
  srt.kim,
  cells = c(
    selected_cells,
    WhichCells(
      srt.kim,
      expression = (Crh > 0 & (Scgn > 0 | Alk > 0 | Fam150b > 0 | Fam150a > 0))
    )
  )
)
# srt.kim <- subset(srt.kim, subset = refUMAP_1 > 4 & refUMAP_2 > -1)

srt.kim@meta.data <- srt.kim@meta.data |>
  dplyr::rename(wtree = predicted.id, age = Age)

srt.kim <- subset(srt.kim, subset = stage %in% c("Pubertal", "Adult"))

srt.kim
```

::: {.cell-output .cell-output-stdout}

```
An object of class Seurat 
27998 features across 954 samples within 1 assay 
Active assay: RNA (27998 features, 3000 variable features)
 3 layers present: counts, data, scale.data
 3 dimensional reductions calculated: umap, ref.pca, ref.umap
```


:::
:::




We subset Kim et al., 2020 dataset to only Pubertal and Adult stages.

## Intersection sets analysis

### PVN Neurons from Kim et al. 2020, Nature Communications




::: {.cell layout-align="center"}

```{.r .cell-code}
sbs_mtx <-
  srt.kim@assays$RNA@data %>%
  as.data.frame() %>%
  t()

rownames(sbs_mtx) <- colnames(srt.kim)
colnames(sbs_mtx) <- rownames(srt.kim)

# Filter features
filt_low_genes <-
  colSums(sbs_mtx) %>%
  .[. > quantile(., 0.4)] %>%
  names()
sbs_mtx %<>% .[, filt_low_genes]

min_filt_vector2 <-
  sbs_mtx %>%
  as_tibble() %>%
  dplyr::select(all_of(filt_low_genes)) %>%
  summarise(across(.fns = ~ quantile(.x, .005))) %>%
  as.list() %>%
  map(as.double) %>%
  simplify() %>%
  .[filt_low_genes]

# Prepare table of intersection sets analysis
content_sbs_mtx_kim <-
  (sbs_mtx > min_filt_vector2) %>%
  as_tibble() %>%
  mutate_all(as.numeric) %>%
  bind_cols(
    srt.kim@meta.data |> dplyr::select(wtree, age, stage)
  )
```
:::




### PVN Neurons from Romanov et al. 2020, Nature




::: {.cell layout-align="center"}

```{.r .cell-code}
srt.romanov.pvn <-
  subset(
    x = srt.romanov.pub,
    cells = unique(c(
      WhichCells(
        srt.romanov.pub,
        idents = c(
          "mneOXY",
          "mneVAS",
          "pneSS",
          "pneCRH",
          "pneTRH"
        )
      ),
      WhichCells(
        srt.romanov.pub,
        expression = (Crh > 0 &
          (Scgn > 0 | Alk > 0 | Fam150b > 0 | Fam150a > 0))
      )
    )),
    invert = FALSE
  )

table(srt.romanov.pvn$age)
```

::: {.cell-output .cell-output-stdout}

```

E15 E17 P00 P02 P10 P23 
332 615 374 389 537  67 
```


:::

```{.r .cell-code}
table(srt.romanov.pvn$stage)
```

::: {.cell-output .cell-output-stdout}

```

Embryonic  Neonatal  Pubertal     Adult 
      947       763       537        67 
```


:::

```{.r .cell-code}
srt.romanov.pvn <- subset(
  srt.romanov.pvn,
  subset = stage %in% c("Pubertal", "Adult")
)
```
:::




We subset Romanov et al., 2020 dataset to only Pubertal and Adult
stages.




::: {.cell layout-align="center"}

```{.r .cell-code}
table(srt.romanov.pvn$age)
```

::: {.cell-output .cell-output-stdout}

```

P10 P23 
537  67 
```


:::

```{.r .cell-code}
table(srt.romanov.pvn$stage)
```

::: {.cell-output .cell-output-stdout}

```

Embryonic  Neonatal  Pubertal     Adult 
        0         0       537        67 
```


:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
sbs_mtx <-
  srt.romanov.pvn@assays$RNA@data %>%
  as.data.frame() %>%
  t()

# Filter features
filt_low_genes <-
  colSums(sbs_mtx) %>%
  .[. > quantile(., 0.4)] %>%
  names()
sbs_mtx %<>% .[, c("Fam150b", filt_low_genes)]

min_filt_vector2 <-
  sbs_mtx %>%
  as_tibble() %>%
  dplyr::select(all_of(c("Fam150b", filt_low_genes))) %>%
  summarise(across(.fns = ~ quantile(.x, .005))) %>%
  as.list() %>%
  map(as.double) %>%
  simplify() %>%
  .[c("Fam150b", filt_low_genes)]

# Prepare table of intersection sets analysis
content_sbs_mtx_romanov <-
  (sbs_mtx > min_filt_vector2) %>%
  as_tibble() %>%
  mutate_all(as.numeric) %>%
  bind_cols(
    srt.romanov.pvn@meta.data |> dplyr::select(wtree, age, stage)
  )
```
:::




### PVN Neurons from both datasets joined




::: {.cell layout-align="center"}

```{.r .cell-code}
# Prepare table of intersection sets analysis
to_select <-
  c(gene.scale, "wtree", "age", "stage") %>%
  .[. %in% colnames(content_sbs_mtx_kim)] %>%
  .[. %in% colnames(content_sbs_mtx_romanov)]

content_sbs_mtx <-
  bind_rows(
    content_sbs_mtx_kim |> dplyr::select(all_of(to_select)),
    content_sbs_mtx_romanov |> dplyr::select(all_of(to_select))
  )
```
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
metabolic_signaling_genes <- c(
  "Alk",
  "Fam150a",
  "Fam150b",
  "Mc3r",
  "Mc4r",
  "Lepr",
  "Insr",
  "Lmo4",
  "Irs1",
  "Irs4"
)

# 1. Prepare the data for the UpSet plot (specifically for Pubertal stage)
data_for_upset_pubertal <- content_sbs_mtx |>
  dplyr::filter(stage == "Pubertal") |>
  # Select the genes of interest, including Fam150b via metabolic_signaling_genes
  dplyr::select(dplyr::any_of(c(
    metabolic_signaling_genes,
    "Crh",
    "Trh",
    "Oxt"
  ))) |>
  # IMPORTANT: This step removes any gene column that has no expressing cells in this subset
  dplyr::select(dplyr::where(~ is.numeric(.x) && sum(.x, na.rm = TRUE) > 0))

# Convert to a standard data frame for UpSetR
data_for_upset_pubertal_df <- as.data.frame(data_for_upset_pubertal)

# 2. Check if Fam150b is present after filtering
if ("Fam150b" %in% colnames(data_for_upset_pubertal_df)) {
  cat("Fam150b is present in the data for the UpSet plot.\n")
  cat(
    "Number of cells expressing Fam150b in Young mice hypothalamic neurons (after binarization):",
    sum(data_for_upset_pubertal_df$Fam150b, na.rm = TRUE),
    "\n"
  )
} else {
  cat("Fam150b is NOT present in the data for the UpSet plot.\n")
  cat(
    "This means no cells in the Young stage expressed Fam150b. \n"
  )

  # Optional: Check Fam150b in the original Pubertal subset before the sum filter
  pubertal_subset_before_sum_filter <- content_sbs_mtx |>
    dplyr::filter(stage == "Pubertal") |>
    dplyr::select(dplyr::any_of("Fam150b"))

  if ("Fam150b" %in% colnames(pubertal_subset_before_sum_filter)) {
    cat(
      "Original expression of Fam150b in Young mice hypothalamic neurons (sum after binarization):",
      sum(pubertal_subset_before_sum_filter$Fam150b, na.rm = TRUE),
      "\n"
    )
  } else {
    cat(
      "Fam150b was not even found in the selection from content_sbs_mtx for Young mice hypothalamic neurons.\n"
    )
  }
}
```

::: {.cell-output .cell-output-stdout}

```
Fam150b is present in the data for the UpSet plot.
Number of cells expressing Fam150b in Young mice hypothalamic neurons (after binarization): 5 
```


:::

```{.r .cell-code}
# 3. Generate the UpSet plot using the prepared and filtered data
# The 'sets' argument should dynamically get column names from the prepared data.
if (
  ncol(data_for_upset_pubertal_df) > 0 && nrow(data_for_upset_pubertal_df) > 0
) {
  upset(
    data_for_upset_pubertal_df,
    # --- Ordering and Filtering ---
    order.by = "freq", # Order intersections by frequency (size)
    # keep.order = TRUE, # If you want to maintain the order of sets as in 'sets' argument
    nsets = 30, # Max number of sets to display in the matrix
    nintersects = 100, # Max number of intersections to display
    # cutoff = 1,        # Minimum size of an intersection to be shown

    # --- Aesthetics: Colors and Points/Lines ---
    main.bar.color = "gray20", # Color of the main intersection size bars
    sets.bar.color = "black", # Color of the set size bars on the left
    matrix.color = "red", # Color of the lines connecting the dots in the matrix
    matrix.dot.alpha = 0.5, # Transparency of the dots (0=transparent, 1=opaque)
    # The example shows solid dots, so alpha close to 1.
    # Note: matrix.color might also affect dot color if specific dot color param is not available or used.
    # UpsetR's parameters can sometimes be a bit tricky.
    shade.color = "gray80", # Color of the shading for row highlighting (if enabled)
    shade.alpha = 0.25, # Transparency of the shading

    # --- Aesthetics: Points and Lines ---
    point.size = 10.0, # Size of the dots in the matrix (adjust as needed)
    line.size = 4.0, # Thickness of the lines in the matrix (adjust as needed)

    # --- Text and Labels ---
    sets.x.label = "Number of cells", # Label for the set size bars
    mainbar.y.label = "Intersection size", # Label for the main intersection size bars
    text.scale = c(
      intersection_size_title = 4.0, # Scales various text elements
      intersection_size_tick_labels = 3.5,
      set_size_title = 4.0,
      set_size_tick_labels = 3.5,
      set_names = 3.8, # Size of the set names (gene names)
      numbers_above_bars = 3.0 # Size of numbers above intersection bars
    ),
    number.angles = 0, # Angle of the numbers above intersection bars

    # --- Matrix Configuration ---
    show.numbers = "yes", # Show numbers above intersection bars ("no" to hide)

    # --- Set Configuration ---
    sets = colnames(data_for_upset_pubertal_df), # Explicitly define the sets from your data
    empty.intersections = NULL,

    # --- Layout Ratio ---
    mb.ratio = c(0.70, 0.30) # Ratio height_of_main_bar_plot / height_of_matrix (manual default)
  )
} else {
  cat(
    "The data frame for the UpSet plot is empty (either no rows or no columns after filtering).\n"
  )
  cat("Number of rows:", nrow(data_for_upset_pubertal_df), "\n")
  cat("Number of columns:", ncol(data_for_upset_pubertal_df), "\n")
}


# 4. Display skim output for the data being plotted
print(skim(data_for_upset_pubertal_df))
```

::: {.cell-output .cell-output-stdout}

```
── Data Summary ────────────────────────
                           Values                      
Name                       data_for_upset_pubertal_d...
Number of rows             753                         
Number of columns          11                          
_______________________                                
Column type frequency:                                 
  numeric                  11                          
________________________                               
Group variables            None                        

── Variable type: numeric ──────────────────────────────────────────────────────
   skim_variable n_missing complete_rate    mean     sd p0 p25 p50 p75 p100
 1 Alk                   0             1 0.369   0.483   0   0   0   1    1
 2 Fam150b               0             1 0.00664 0.0813  0   0   0   0    1
 3 Mc3r                  0             1 0.0146  0.120   0   0   0   0    1
 4 Mc4r                  0             1 0.0332  0.179   0   0   0   0    1
 5 Lepr                  0             1 0.0279  0.165   0   0   0   0    1
 6 Lmo4                  0             1 0.284   0.451   0   0   0   1    1
 7 Irs1                  0             1 0.0916  0.289   0   0   0   0    1
 8 Irs4                  0             1 0.284   0.451   0   0   0   1    1
 9 Crh                   0             1 0.161   0.367   0   0   0   0    1
10 Trh                   0             1 0.325   0.469   0   0   0   1    1
11 Oxt                   0             1 0.390   0.488   0   0   0   1    1
   hist 
 1 ▇▁▁▁▅
 2 ▇▁▁▁▁
 3 ▇▁▁▁▁
 4 ▇▁▁▁▁
 5 ▇▁▁▁▁
 6 ▇▁▁▁▃
 7 ▇▁▁▁▁
 8 ▇▁▁▁▃
 9 ▇▁▁▁▂
10 ▇▁▁▁▃
11 ▇▁▁▁▅
```


:::

::: {.cell-output-display}
![](03-upset_files/figure-html/fig-upset-group-metabolic-pvn-f1d-1.png){#fig-upset-group-metabolic-pvn-f1d fig-align='center' width=10800}
:::
:::

::: {.cell layout-align="center"}

```{.r .cell-code}
sessioninfo::session_info()
```

::: {.cell-output .cell-output-stdout}

```
─ Session info ───────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.4.2 (2024-10-31)
 os       Ubuntu 22.04.5 LTS
 system   x86_64, linux-gnu
 ui       X11
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Europe/Vienna
 date     2025-05-06
 pandoc   3.4 @ /data/1_heteroAstrocytes/.pixi/envs/default/bin/ (via rmarkdown)
 quarto   1.6.41 @ /data/1_heteroAstrocytes/.pixi/envs/default/bin/quarto

─ Packages ───────────────────────────────────────────────────────────────────
 package          * version    date (UTC) lib source
 abind              1.4-5      2016-07-21 [1] CRAN (R 4.4.1)
 base64enc          0.1-3      2015-07-28 [1] CRAN (R 4.4.1)
 bayestestR         0.15.3     2025-04-28 [1] CRAN (R 4.4.3)
 beeswarm           0.4.0      2021-06-01 [1] RSPM
 BiocManager        1.30.25    2024-08-28 [1] CRAN (R 4.4.1)
 bit                4.5.0.1    2024-12-03 [1] CRAN (R 4.4.2)
 bit64              4.5.2      2024-09-22 [1] CRAN (R 4.4.1)
 circlize           0.4.16     2024-02-20 [1] RSPM
 cli                3.6.4      2025-02-13 [1] CRAN (R 4.4.2)
 cluster            2.1.8      2024-12-11 [1] CRAN (R 4.4.2)
 codetools          0.2-20     2024-03-31 [1] CRAN (R 4.4.1)
 colorspace         2.1-1      2024-07-26 [1] CRAN (R 4.4.1)
 correlation        0.8.7      2025-03-03 [1] CRAN (R 4.4.3)
 cowplot          * 1.1.3      2024-01-22 [1] CRAN (R 4.4.1)
 crayon             1.5.3      2024-06-20 [1] CRAN (R 4.4.1)
 data.table         1.17.0     2025-02-22 [1] CRAN (R 4.4.2)
 datawizard         1.0.2      2025-03-24 [1] CRAN (R 4.4.3)
 deldir             2.0-4      2024-02-28 [1] CRAN (R 4.4.1)
 digest             0.6.37     2024-08-19 [1] CRAN (R 4.4.1)
 dotCall64          1.2        2024-10-04 [1] CRAN (R 4.4.1)
 dplyr            * 1.1.4      2023-11-17 [1] CRAN (R 4.4.1)
 effectsize         1.0.0      2024-12-10 [1] CRAN (R 4.4.2)
 evaluate           1.0.3      2025-01-10 [1] CRAN (R 4.4.2)
 farver             2.1.2      2024-05-13 [1] CRAN (R 4.4.1)
 fastDummies        1.7.5      2025-01-20 [1] CRAN (R 4.4.2)
 fastmap            1.2.0      2024-05-15 [1] CRAN (R 4.4.1)
 fitdistrplus       1.2-2      2025-01-07 [1] CRAN (R 4.4.2)
 forcats          * 1.0.0      2023-01-29 [1] CRAN (R 4.4.1)
 future           * 1.34.0     2024-07-29 [1] CRAN (R 4.4.1)
 future.apply       1.11.3     2024-10-27 [1] CRAN (R 4.4.2)
 generics           0.1.3      2022-07-05 [1] CRAN (R 4.4.1)
 ggbeeswarm         0.7.2      2023-04-29 [1] RSPM
 ggmin              0.0.0.9000 2025-03-07 [1] Github (sjessa/ggmin@8ada274)
 ggplot2          * 3.5.1      2024-04-23 [1] CRAN (R 4.4.1)
 ggprism            1.0.5      2024-03-21 [1] RSPM
 ggrastr            1.0.2      2023-06-01 [1] RSPM
 ggrepel            0.9.6      2024-09-07 [1] CRAN (R 4.4.1)
 ggridges           0.5.6      2024-01-23 [1] CRAN (R 4.4.1)
 ggstatsplot      * 0.13.0     2024-12-04 [1] CRAN (R 4.4.2)
 GlobalOptions      0.1.2      2020-06-10 [1] RSPM
 globals            0.16.3     2024-03-08 [1] CRAN (R 4.4.1)
 glue               1.8.0      2024-09-30 [1] CRAN (R 4.4.1)
 goftest            1.2-3      2021-10-07 [1] CRAN (R 4.4.1)
 gridExtra          2.3        2017-09-09 [1] CRAN (R 4.4.1)
 gtable             0.3.6      2024-10-25 [1] CRAN (R 4.4.1)
 here             * 1.0.1      2020-12-13 [1] CRAN (R 4.4.1)
 hms                1.1.3      2023-03-21 [1] CRAN (R 4.4.1)
 htmltools          0.5.8.1    2024-04-04 [1] CRAN (R 4.4.1)
 htmlwidgets        1.6.4      2023-12-06 [1] CRAN (R 4.4.1)
 httpuv             1.6.15     2024-03-26 [1] CRAN (R 4.4.1)
 httr               1.4.7      2023-08-15 [1] CRAN (R 4.4.1)
 ica                1.0-3      2022-07-08 [1] CRAN (R 4.4.1)
 igraph             2.0.3      2024-03-13 [1] CRAN (R 4.4.1)
 insight            1.2.0      2025-04-22 [1] CRAN (R 4.4.3)
 irlba              2.3.5.1    2022-10-03 [1] CRAN (R 4.4.1)
 janitor            2.2.1      2024-12-22 [1] RSPM
 jsonlite           1.9.1      2025-03-03 [1] CRAN (R 4.4.3)
 KernSmooth         2.23-26    2025-01-01 [1] CRAN (R 4.4.2)
 knitr              1.49       2024-11-08 [1] CRAN (R 4.4.1)
 labeling           0.4.3      2023-08-29 [1] CRAN (R 4.4.1)
 later              1.4.1      2024-11-27 [1] CRAN (R 4.4.2)
 lattice            0.22-6     2024-03-20 [1] CRAN (R 4.4.1)
 lazyeval           0.2.2      2019-03-15 [1] CRAN (R 4.4.1)
 lifecycle          1.0.4      2023-11-07 [1] CRAN (R 4.4.1)
 listenv            0.9.1      2024-01-29 [1] CRAN (R 4.4.1)
 lmtest             0.9-40     2022-03-21 [1] CRAN (R 4.4.1)
 lubridate        * 1.9.4      2024-12-08 [1] CRAN (R 4.4.2)
 magrittr         * 2.0.3      2022-03-30 [1] CRAN (R 4.4.1)
 MASS               7.3-64     2025-01-04 [1] CRAN (R 4.4.2)
 Matrix             1.7-2      2025-01-23 [1] CRAN (R 4.4.2)
 matrixStats        1.5.0      2025-01-07 [1] CRAN (R 4.4.2)
 mime               0.12       2021-09-28 [1] CRAN (R 4.4.1)
 miniUI             0.1.1.1    2018-05-18 [1] CRAN (R 4.4.1)
 munsell            0.5.1      2024-04-01 [1] CRAN (R 4.4.1)
 nlme               3.1-167    2025-01-27 [1] CRAN (R 4.4.2)
 paletteer          1.6.0      2024-01-21 [1] CRAN (R 4.4.1)
 parallelly         1.42.0     2025-01-30 [1] CRAN (R 4.4.2)
 parameters         0.25.0     2025-04-30 [1] CRAN (R 4.4.3)
 patchwork        * 1.3.0      2024-09-16 [1] CRAN (R 4.4.1)
 pbapply            1.7-2      2023-06-27 [1] CRAN (R 4.4.1)
 pillar             1.10.1     2025-01-07 [1] CRAN (R 4.4.2)
 pkgconfig          2.0.3      2019-09-22 [1] CRAN (R 4.4.1)
 plotly             4.10.4     2024-01-13 [1] CRAN (R 4.4.1)
 plyr               1.8.9      2023-10-02 [1] CRAN (R 4.4.1)
 png                0.1-8      2022-11-29 [1] CRAN (R 4.4.1)
 polyclip           1.10-7     2024-07-23 [1] CRAN (R 4.4.1)
 progressr          0.15.1     2024-11-22 [1] CRAN (R 4.4.2)
 promises           1.3.2      2024-11-28 [1] CRAN (R 4.4.2)
 purrr            * 1.0.4      2025-02-05 [1] CRAN (R 4.4.2)
 R.methodsS3        1.8.2      2022-06-13 [1] RSPM
 R.oo               1.27.0     2024-11-01 [1] RSPM
 R.utils            2.13.0     2025-02-24 [1] RSPM
 R6                 2.6.1      2025-02-15 [1] CRAN (R 4.4.2)
 RANN               2.6.2      2024-08-25 [1] CRAN (R 4.4.1)
 RColorBrewer     * 1.1-3      2022-04-03 [1] CRAN (R 4.4.1)
 Rcpp               1.0.14     2025-01-12 [1] CRAN (R 4.4.2)
 RcppAnnoy          0.0.22     2024-01-23 [1] CRAN (R 4.4.1)
 RcppHNSW           0.6.0      2024-02-04 [1] CRAN (R 4.4.1)
 RcppParallel       5.1.9      2024-08-19 [1] CRAN (R 4.4.1)
 readr            * 2.1.5      2024-01-10 [1] CRAN (R 4.4.1)
 rematch2           2.1.2      2020-05-01 [1] CRAN (R 4.4.1)
 remotes            2.5.0      2024-03-17 [1] CRAN (R 4.4.1)
 repr               1.1.7      2024-03-22 [1] CRAN (R 4.4.1)
 reshape2           1.4.4      2020-04-09 [1] CRAN (R 4.4.1)
 reticulate         1.40.0     2024-11-15 [1] CRAN (R 4.4.1)
 rhdf5              2.46.1     2023-11-29 [1] Bioconductor 3.18 (R 4.4.2)
 rhdf5filters       1.14.1     2023-11-06 [1] Bioconductor
 Rhdf5lib           1.24.2     2024-02-07 [1] Bioconductor 3.18 (R 4.4.2)
 rlang              1.1.5      2025-01-17 [1] CRAN (R 4.4.2)
 rmarkdown          2.29       2024-11-04 [1] CRAN (R 4.4.1)
 ROCR               1.0-11     2020-05-02 [1] CRAN (R 4.4.1)
 rprojroot          2.0.4      2023-11-05 [1] CRAN (R 4.4.1)
 RSpectra           0.16-2     2024-07-18 [1] CRAN (R 4.4.1)
 rstantools         2.4.0      2024-01-31 [1] CRAN (R 4.4.1)
 rsvd               1.0.5      2021-04-16 [1] RSPM
 Rtsne              0.17       2023-12-07 [1] CRAN (R 4.4.1)
 scales             1.3.0      2023-11-28 [1] CRAN (R 4.4.1)
 scattermore        1.2        2023-06-12 [1] CRAN (R 4.4.1)
 scCustomize      * 3.0.1      2025-03-07 [1] Github (samuel-marsh/scCustomize@3299b95)
 schard             0.0.1      2025-05-06 [1] Github (cellgeni/schard@1b62a7c)
 sctransform        0.4.1      2023-10-19 [1] CRAN (R 4.4.1)
 sessioninfo        1.2.3      2025-02-05 [1] CRAN (R 4.4.2)
 Seurat           * 5.2.1      2025-01-24 [1] CRAN (R 4.4.2)
 SeuratObject     * 5.0.2      2024-05-08 [1] CRAN (R 4.4.3)
 SeuratWrappers   * 0.4.0      2025-03-07 [1] Github (satijalab/seurat-wrappers@a1eb0d8)
 shape              1.4.6.1    2024-02-23 [1] CRAN (R 4.4.1)
 shiny              1.10.0     2024-12-14 [1] CRAN (R 4.4.2)
 skimr            * 2.1.5      2022-12-23 [1] CRAN (R 4.4.1)
 snakecase          0.11.1     2023-08-27 [1] RSPM
 sp               * 2.2-0      2025-02-01 [1] CRAN (R 4.4.2)
 spam               2.11-1     2025-01-20 [1] CRAN (R 4.4.2)
 spatstat.data      3.1-4      2024-11-15 [1] CRAN (R 4.4.2)
 spatstat.explore   3.3-4      2025-01-08 [1] CRAN (R 4.4.2)
 spatstat.geom      3.3-5      2025-01-18 [1] CRAN (R 4.4.2)
 spatstat.random    3.3-2      2024-09-18 [1] CRAN (R 4.4.2)
 spatstat.sparse    3.1-0      2024-06-21 [1] CRAN (R 4.4.1)
 spatstat.univar    3.1-1      2024-11-05 [1] CRAN (R 4.4.2)
 spatstat.utils     3.1-2      2025-01-08 [1] CRAN (R 4.4.2)
 statsExpressions   1.6.2      2024-12-02 [1] CRAN (R 4.4.2)
 stringi            1.8.4      2024-05-06 [1] CRAN (R 4.4.1)
 stringr          * 1.5.1      2023-11-14 [1] CRAN (R 4.4.1)
 survival           3.8-3      2024-12-17 [1] CRAN (R 4.4.2)
 tensor             1.5        2012-05-05 [1] CRAN (R 4.4.1)
 tibble           * 3.2.1      2023-03-20 [1] CRAN (R 4.4.1)
 tidyr            * 1.3.1      2024-01-24 [1] CRAN (R 4.4.1)
 tidyselect         1.2.1      2024-03-11 [1] CRAN (R 4.4.1)
 tidyverse        * 2.0.0      2023-02-22 [1] CRAN (R 4.4.1)
 timechange         0.3.0      2024-01-18 [1] CRAN (R 4.4.1)
 tzdb               0.4.0      2023-05-12 [1] CRAN (R 4.4.1)
 UpSetR           * 1.4.0      2019-05-22 [1] CRAN (R 4.4.2)
 utf8               1.2.4      2023-10-22 [1] CRAN (R 4.4.1)
 uwot               0.2.3      2025-02-24 [1] CRAN (R 4.4.2)
 vctrs              0.6.5      2023-12-01 [1] CRAN (R 4.4.1)
 vipor              0.4.7      2023-12-18 [1] RSPM
 viridis          * 0.6.5      2024-01-29 [1] RSPM
 viridisLite      * 0.4.2      2023-05-02 [1] CRAN (R 4.4.1)
 vroom              1.6.5      2023-12-05 [1] CRAN (R 4.4.1)
 withr              3.0.2      2024-10-28 [1] CRAN (R 4.4.1)
 xfun               0.51       2025-02-19 [1] CRAN (R 4.4.2)
 xtable             1.8-4      2019-04-21 [1] CRAN (R 4.4.1)
 yaml               2.3.10     2024-07-26 [1] CRAN (R 4.4.1)
 zeallot            0.1.0      2018-01-28 [1] CRAN (R 4.4.1)
 zoo                1.8-13     2025-02-22 [1] CRAN (R 4.4.2)

 [1] /data/1_heteroAstrocytes/.pixi/envs/default/lib/R/library
 * ── Packages attached to the search path.

──────────────────────────────────────────────────────────────────────────────
```


:::
:::
