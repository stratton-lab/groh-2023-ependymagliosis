library(purrr)
library(stringr)

if (!dir.exists("data")) {
  dir.create("data")
  dir.create("data/GSE199460")
  dir.create("data/GSE254863")
}
options(timeout = 600)

# fetch GSE199460 data
download.file(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE199460&format=file",
  "data/GSE199460.tar",
  quiet = TRUE
)
untar("data/GSE199460.tar", exdir = "data/GSE199460")
file.remove("data/GSE199460.tar")
walk(
  list.files(
    "data/GSE199460",
    pattern = "*.tar.gz",
    full.names = TRUE
  ),
  function(f) {
    untar(f, exdir = "data/GSE199460")
    file.remove(f)
  }
)
download.file(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE199460&format=file&file=GSE199460%5Fcell%5Fannotation%2Emeta%5Fdata%2Etotal%5Fcells%2Ecsv%2Egz",
  "data/GSE199460/meta.csv.gz",
  quiet = TRUE
)

# fetch GSE254863 data
download.file(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE254863&format=file",
  "data/GSE254863.tar",
  quiet = TRUE
)
untar("data/GSE254863.tar", exdir = "data/GSE254863")
file.remove("data/GSE254863.tar")
dir.create("data/GSE254863/EAE4p")
dir.create("data/GSE254863/EAE4n")
iwalk(list(".+_Pos_" = "EAE4p", ".+_Neg_" = "EAE4n"), function(d, p) {
  walk(
    list.files(path = "data/GSE254863", pattern = p),
    function(f) {
      file.rename(
        sprintf("data/GSE254863/%s", f),
        sprintf("data/GSE254863/%s/%s", d, str_replace(f, p, ""))
      )
    }
  )
})
download.file(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE254863&format=file&file=GSE254863%5Fannotations%2Ecsv%2Egz",
  "data/GSE254863/meta.csv.gz",
  quiet = TRUE
)
