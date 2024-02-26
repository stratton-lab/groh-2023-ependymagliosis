options("future.globals.maxSize" = 4 * 10 ^ 9)
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(enrichplot)
  library(patchwork)
  library(clusterProfiler)
})

MARKERS <- readRDS("out/mlists.rds")

BRNI <- readRDS("out/brni.rds") %>%
  SetIdent(value = "cell_type")

COL_MAP <- list(
  "Astrocytes" = "#aec7e8",
  "Endothelial Cells" = "#ff7f0e",
  "Neuroblasts" = "#8c564b",
  "Ependymal Cells" = "#2ca02c",
  "Choroid Plexus Cells" = "#d62728",
  "Pericytes" = "#17becf",
  "Neurons" = "#9467bd",
  "Fibroblasts" = "#dbdb8d"
)
FONT_SIZE <- 14
DESIGN <- "
  AAAA####
  AAAA####
  AAAA####
  CCCCDDDD
  CCCCDDDD
  CCCCDDDD
  CCCCDDDD
  CCCCDDDD
  ##GGGGGG
  FFGGGGGG
  FFGGGGGG
  FFGGGGGG
  FFGGGGGG
  HHHHHHHH
  HHHHHHHH
  HHHHHHHH
  HHHHHHHH
  HHHHHHHH
  HHHHHHHH
  HHHHHHHH
"

DROP_DESCRIPTIONS <- c(
  "Response to biotic stimulus",
  "Defense response to virus",
  "Cellular response to interferon beta",
  "Regulation of immune system process",
  "Innate immune response",
  "Defense response to other organism",
  "Immune response",
  "Regulation of response to external stimulus",
  "Plasma membrane bounded cell projection organization",
  "Regulation of mitochondrial fission",
  "Cellular response to molecule of bacterial origin"
)

# table 1
CO_REGULATED <- data.frame(arrange(
  merge(
    MARKERS[["e_gliosis"]],
    MARKERS[["a_gliosis"]],
    by = "row.names",
    suffix = c(".ependyma", ".astrocyte")
  ),
  desc(avg_log2FC.ependyma)
),
row.names = "Row.names")
write.csv(CO_REGULATED, "out/table1.csv")

# table 2
CO_REGULATED_CORE <- subset(
  CO_REGULATED,
  rownames(CO_REGULATED) %in% rownames(MARKERS[["e_gliosis_ctx"]]) &
    rownames(CO_REGULATED) %in% rownames(MARKERS[["a_gliosis_ctx"]])
)
write.csv(CO_REGULATED_CORE, "out/table2.csv")

mk_A <- function(so, cols) {
  DimPlot(
    so,
    cols = cols,
    group.by = "cell_type",
    label = TRUE,
    reduction = "umap"
  )
}

mk_B <- function() {
  # TODO: figure out venn diagram
}

mk_C <- function(so, cols) {
  wrap_plots(
    VlnPlot(so, "Aldh1l1", cols = cols) & labs(x = NULL),
    VlnPlot(so, "Foxj1", cols = cols) & labs(x = NULL),
    design = "AB",
    guides = "collect"
  )
}

compute_ratio <- function(gr) {
  str_split(gr, "/") %>%
    map_vec(compose(partial(do.call, `/`),
                    as.list,
                    as.numeric))
}

mk_goterm_dotplot <- function(markers, add_args) {
  do.call(enrichGO,
          append(
            list(markers, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL"),
            add_args
          )) %>%
    simplify() %>%
    data.frame() %>%
    subset(!(Description %in% DROP_DESCRIPTIONS)) %>%
    slice_min(qvalue, n = 35) %>%
    ggplot() +
    aes(
      x = -log(qvalue),
      y = fct_reorder(Description, qvalue, .desc = TRUE),
      size = compute_ratio(GeneRatio)
    ) +
    geom_point() +
    guides(size = guide_legend(title = "gene ratio")) +
    theme_classic() +
    labs(y = NULL)
}

mk_D <- function(marker_list) {
  mk_goterm_dotplot(marker_list,
                    add_args = list(
                      ont = "BP",
                      minGSSize = 2,
                      maxGSSize = 5000
                    )) + labs(title = "upregulated DEG over-represented GO terms")
}

mk_log2FC_plot <- function(markers) {
  ggplot(markers) +
    aes(
      y = fct_reorder(rownames(markers), avg_log2FC, .desc = TRUE),
      x = avg_log2FC,
      fill = pct.1 - pct.2
    ) +
    geom_bar(stat = "identity") +
    scale_fill_viridis_c(name = "pct diff",
                         option = "turbo",
                         limits = c(-0.1, 0.9)) +
    lims(x = c(0, 8)) +
    labs(y = NULL) +
    theme_classic()
}

mk_F <- function(markers, genes_table, vsplit = FALSE) {
  (if (vsplit) {
    inner_join(genes_table, mutate(markers, gene = rownames(markers)), by = "gene")
  } else {
    inner_join(
      full_join(
        genes_table,
        aggregate(category ~ gene, genes_table, length),
        by = "gene",
        suffix = c("", ".count")
      ),
      
      mutate(markers, gene = rownames(markers))[c("gene", "avg_log2FC.ependyma")],
      by = "gene"
    ) %>%
      mutate(avg_log2FC.ependyma = avg_log2FC.ependyma / category.count)
  }) %>%
    ggplot() +
    aes(fill = category, x = avg_log2FC.ependyma, y = gene) +
    geom_bar(position = if (vsplit) {
      "dodge"
    } else {
      "stack"
    },
    stat = "identity") +
    theme_classic() +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 2))
}

mk_G <- function(e_gliosis, genes_table) {
  e_gliosis <- mutate(e_gliosis, gene = rownames(e_gliosis))
  list(
    list(
      "unique reactive\nastrocyte genes",
      anti_join(genes_table,
                e_gliosis,
                by = "gene")
    ),
    list(
      "shared reactive\ngenes",
      semi_join(genes_table,
                e_gliosis,
                by = "gene")
    ),
    list(
      "unique reactive\nependymal genes",
      anti_join(e_gliosis,
                genes_table,
                by = "gene")
    )
  ) %>%
    imap(partial(do.call, function(t, g) {
      mk_goterm_dotplot(g[["gene"]],
                        add_args = list(
                          minGSSize = 2,
                          maxGSSize = 5000,
                          ont = "BP"
                        )) +
        labs(title = t) +
        lims(size = c(0, 1),
             x = c(0, 80))
    })) %>%
    wrap_plots(... = .,
               nrow = 1,
               guides = "collect")
}

mk_H <- function(so, marker_list) {
  c("Ependymal Cells", "Astrocytes", "Choroid Plexus Cells") %>%
    map(function(ident) {
      subset(so, cell_type == ident) %>%
        FoldChange(
          features = marker_list,
          ident.1 = "E",
          ident.2 = "C",
          group.by = "condition",
          verbose = FALSE
        ) %>%
        mutate(., gene = rownames(.)) %>%
        ggplot() +
        aes(
          y = fct_reorder(gene, avg_log2FC, .desc = TRUE),
          x = pct.1 - pct.2,
          fill = avg_log2FC,
        ) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_c(name = "average log2FC",
                             option = "turbo",
                             limits = c(-2, 15)) +
        labs(y = NULL) +
        lims(x = c(-0.1, 1)) +
        labs(title = sprintf("expression differences in %s", str_to_lower(ident)),
             x = "delta pct expression") +
        theme_classic()
    }) %>%
    wrap_plots(... = .,
               nrow = 1,
               guides = "collect")
}

if (TRUE) {
  gs <- partial(function(design, n, px, ext = "tiff") {
    clean <- function(s) {
      str_split_1(s, "\n *") %>%
        .[nzchar(.)]
    }
    count <- function(s) {
      str_extract(s, sprintf("%s+", n)) %>%
        map_vec(nchar) %>%
        max(na.rm = TRUE)
    }
    ggsave(
      sprintf("out/fig1%s.%s", str_to_lower(n), ext),
      w = px *
        8.5 /
        nchar(str_extract(design, regex(" *(.*)\n$"), group = 1)) *
        count(clean(design)),
      h = px *
        11 /
        (str_count(design, fixed("\n")) - 1) *
        count(reduce(clean(design),
                     function(a, c) {
                       paste(a, str_split_1(c, ""), sep = "")
                     },
                     .init = "")),
      units = "px",
      limitsize = FALSE
    )
  }, design = DESIGN)
  
  f.a <- mk_A(BRNI, COL_MAP) &
    theme(text = element_text(size = FONT_SIZE))
  gs("A", 1000)
  
  f.c <- mk_C(BRNI, COL_MAP) &
    theme(text = element_text(size = FONT_SIZE))
  gs("C", 1000)
  
  f.d <- mk_D(rownames(CO_REGULATED_CORE)) &
    theme(text = element_text(size = FONT_SIZE))
  gs("D", 1000)
  
  f.f <- mk_F(CO_REGULATED, read.csv("data/astro-review.csv")) &
    theme(text = element_text(size = FONT_SIZE))
  gs("F", 1000)
  
  f.g <- mk_G(MARKERS[["e_gliosis"]],
              read.csv("data/astro-review.csv")) &
    theme(text = element_text(size = FONT_SIZE))
  gs("G", 1000)
  
  f.h <- mk_H(BRNI, rownames(CO_REGULATED_CORE)) &
    theme(text = element_text(size = FONT_SIZE))
  gs("H", 1000)
  
  wrap_plots(
    A = f.a,
    # B = mk_B(),
    C = f.c,
    D = f.d,
    # E = mk_E(),
    F = f.f,
    G = f.g,
    H = f.h,
    design = DESIGN
  )
  ggsave("out/fig1.tiff",
         w = 1000 * 8.5,
         h = 1000 * 11,
         units = "px")
}
