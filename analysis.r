suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(tidyverse)
  library(r2r)
})
options("future.globals.maxSize" = 1e9)

SEURAT_VERBOSE <- TRUE
NON_IMMUNE_BRAIN_CELLS <- c(
  "Endothelial Cells",
  "Choroid Plexus Cells",
  "Pericytes",
  "Fibroblasts",
  "Neuroblasts",
  "Neurons"
)

# function to generate Seurat object from GSE199460 data
# - integrates data together from multiple runs, controlling for batch effects
# - adds provided metadata (after conducting a correction step)
mk_GSE199460 <- function() {
  so <-
    list(
      list(c("CTL1_total", "EAE1_total", "EAE2_total"), "GSE199460.1"),
      list(c("CTL2_total", "CTL3_total", "EAE3_total"), "GSE199460.2")
    ) |>
    map(
      partial(
        do.call,
        function(samples, batch) {
          map(
            samples,
            function(sample) {
              sname <- str_replace(sample, regex("_total$"), "")

              paste(
                "data",
                "GSE199460",
                sample,
                "filtered_feature_bc_matrix",
                sep = "/"
              ) |>
                Read10X() |>
                CreateSeuratObject() |>
                AddMetaData(sname, "sample") |>
                RenameCells(
                  add.cell.id = sname
                )
            }
          ) |>
            reduce(merge) |>
            JoinLayers() |>
            AddMetaData(batch, "batch")
        }
      )
    ) |>
    reduce(merge) |>
    JoinLayers()

  # fix broken cell identifiers in provided metadata

  meta <- read.csv("data/GSE199460/meta.csv.gz")[c(
    # "orig.ident",
    "cell_type",
    "condition"
  )]

  # initialize cells hashset
  cells <- hashset()
  walk(Cells(so), partial(insert, cells))

  # id prefixes in cell identifiers were improperly set
  # specifically, CTL2 and CTL3 prefixes were
  # swapped for a small portion of cells.
  # we can detect these incorrect identifiers
  # and swap the sample prefix accordingly,
  # then verify that the new identifier is still valid.
  rownames(meta) <- map(rownames(meta), function(ident) {
    if (!has_key(cells, ident)) {
      split <- str_split_1(ident, fixed("_"))
      new <- if (split[[1]] == "CTL2") {
        paste("CTL3", split[[2]], sep = "_")
      } else if (split[[1]] == "CTL3") {
        paste("CTL2", split[[2]], sep = "_")
      } else {
        # broken ident is not prefixed with CTL2/CTL3, should never happen
        stop("error: unrecoverable ident")
      }
      if (!has_key(cells, new)) {
        # sanity check to make sure new ident is still valid
        stop("error: unrecoverable ident")
      }
      new
    } else {
      ident
    }
  })

  # drop cells not given annotations
  so[, rownames(meta)] |>
    AddMetaData(meta) |>
    SetIdent(value = "cell_type") |>
    AddMetaData("GSE199460", "source")
}

mk_GSE254863 <- function() {
  c("EAE4p", "EAE4n") |>
    map(
      function(sample) {
        paste(
          "data",
          "GSE254863",
          sample,
          sep = "/"
        ) |>
          Read10X() |>
          CreateSeuratObject() |>
          AddMetaData(sample, "sample") |>
          RenameCells(add.cell.id = sample)
      }
    ) |>
    reduce(merge) |>
    JoinLayers() |>
    AddMetaData("GSE254863.1", "batch") |>
    AddMetaData("GSE254863", "source") |>
    AddMetaData("E", "condition") |>
    AddMetaData(
      read.csv("data/GSE254863/meta.csv.gz") |>
        mutate(
          X = str_replace(X, "Neg_", "EAE4n_"),
          X = str_replace(X, "Pos_", "EAE4p_")
        ) |>
        data.frame(row.names = "X")
    )
}

process_so <- function(so,
                       cell_types = NULL,
                       split_on = NULL,
                       integration_method = HarmonyIntegration,
                       clustering_alg = "leiden") {
  if (!is.null(cell_types)) {
    so <- subset(
      so,
      cell_type %in% cell_types
    )
  }
  if (!is.null(split_on)) {
    # split layers based on variable provided
    so[["RNA"]] <- split(so[["RNA"]], so[[split_on]][, ])

    # get batch-effect corrected embeddings
    so <- SCTransform(
      so,
      verbose = SEURAT_VERBOSE
    )
    so <- RunPCA(
      so,
      verbose = SEURAT_VERBOSE
    )
    so <- IntegrateLayers(
      so,
      integration_method,
      normalization.method = "SCT"
    )

    # rejoin layers now that batch-effect
    # corrected embeddings are calculated
    so[["RNA"]] <- JoinLayers(so[["RNA"]])
    reduction <- formals(integration_method)[["new.reduction"]]
  } else {
    # no batch-effect correction
    so <- SCTransform(
      so,
      verbose = SEURAT_VERBOSE
    )
    so <- RunPCA(
      so,
      verbose = SEURAT_VERBOSE
    )

    reduction <- "pca"
  }

  so <- PrepSCTFindMarkers(so, verbose = SEURAT_VERBOSE)

  if (!is.null(clustering_alg)) {
    so <- FindNeighbors(
      so,
      reduction = reduction,
      dims = 1:30,
      verbose = SEURAT_VERBOSE
    )
    so <- FindClusters(
      so,
      algorithm = clustering_alg,
      verbose = SEURAT_VERBOSE
    )
  }

  so <- RunUMAP(
    so,
    reduction = reduction,
    dims = 1:30,
    verbose = SEURAT_VERBOSE
  )

  so
}

mk_markers_list <- function(so) {
  so <- AddMetaData(
    so,
    paste(
      so[["cell_type"]][, ],
      so[["condition"]][, ],
      sep = " "
    ),
    "type.condition"
  )

  # partial application to reduce repetition
  get_markers <- compose(
    partial(subset, ... = , p_val_adj < 0.05),
    partial(
      FindMarkers,
      object = so,
      only.pos = TRUE,
      group.by = "type.condition",
      verbose = SEURAT_VERBOSE,
      # test.use = "MAST",
      recorrect_umi = FALSE
    )
  )

  add_exp_lvl <- function(df, ident) {
    sub_c <- subset(so, cell_type == ident & condition == "C")
    sub_e <- subset(so, cell_type == ident & condition == "E")
    df[["exp.ctl"]] <- rowSums(
      GetAssayData(sub_c, "RNA", "counts")[rownames(df), ] > 0
    ) / ncol(sub_c)
    df[["exp.eae"]] <- rowSums(
      GetAssayData(sub_e, "RNA", "counts")[rownames(df), ] > 0
    ) / ncol(sub_e)
    df
  }

  # make marker lists
  e_gliosis <- add_exp_lvl(
    get_markers(
      ident.1 = "Ependymal Cells E",
      ident.2 = "Ependymal Cells C"
    ), "Ependymal Cells"
  )
  e_gliosis_ctx <- add_exp_lvl(
    get_markers(
      ident.1 = "Ependymal Cells E",
      ident.2 = paste(NON_IMMUNE_BRAIN_CELLS, "E", sep = " ")
    ), "Ependymal Cells"
  )
  a_gliosis <- add_exp_lvl(
    get_markers(
      ident.1 = "Astrocytes E",
      ident.2 = "Astrocytes C"
    ), "Astrocytes"
  )
  a_gliosis_ctx <- add_exp_lvl(
    get_markers(
      ident.1 = "Astrocytes E",
      ident.2 = paste(NON_IMMUNE_BRAIN_CELLS, "E", sep = " ")
    ), "Astrocytes"
  )

  list(
    "a_gliosis" = a_gliosis,
    "a_gliosis_ctx" = a_gliosis_ctx,
    "e_gliosis" = e_gliosis,
    "e_gliosis_ctx" = e_gliosis_ctx
  )
}

if (FALSE) {
  brni <- process_so(
    JoinLayers(merge(mk_GSE199460(), mk_GSE254863())),
    cell_types = c("Astrocytes", "Ependymal Cells", NON_IMMUNE_BRAIN_CELLS),
    split_on = "batch",
    integration_method = HarmonyIntegration
  )

  markers <- mk_markers_list(brni)

  saveRDS(brni, "out/brni.rds")
  saveRDS(markers, "out/mlists.rds")
}
