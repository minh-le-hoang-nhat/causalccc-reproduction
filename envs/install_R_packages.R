# Installs the R packages of this project at fixed versions.
# Run from the repository root:  Rscript envs/install_R_packages.R
# CRAN versions are fixed by the snapshot date in .Rprofile; GitHub packages by commit.

pak::pak(c("Seurat", "tidyverse", "rjson"))   # rjson + tidyverse: needed by miic's causalccc branch but not declared
pak::pak("visNetwork")                         # interactive network plots (notebook 06)
pak::pak("miicTeam/CausalCCC@0aa8586e961153f244dcc1001aba9cb1442102c8")
pak::pak("saezlab/liana@6cab46c54234f861ea176c3de77c4b8aa45ecb3d")  # LIANA: L-R inference used by causalCCC.links()
pak::pak("miicTeam/miic_R_package@5ca0b222cc36bee902aa87e5ecfc480428d98927")  # last, so this exact commit is kept

# Self-check: stops with an error if the wrong versions are installed
stopifnot(
  packageDescription("miic")$RemoteSha == "5ca0b222cc36bee902aa87e5ecfc480428d98927",
  packageDescription("causalCCC")$RemoteSha == "0aa8586e961153f244dcc1001aba9cb1442102c8",
  packageDescription("liana")$RemoteSha == "6cab46c54234f861ea176c3de77c4b8aa45ecb3d"
)
message("R packages installed and checked.")
