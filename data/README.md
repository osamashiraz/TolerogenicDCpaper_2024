# data/

`CellChatDB.human_AAT.rda` goes here. It's the ligand-receptor database used
for the CellChat analysis (step 3): the CellChat v2 human database plus 42 pairs
relevant to tolerogenic DCs, such as IDO1/IDO2–AHR, AHR–TXN, TXN–IDO1, IFNG–IDO1
and LILRB3/LILRB4–CD80/CD86. Loading it gives one object, `db.AAT`.

It's archived on Zenodo with the other analysis objects
([10.5281/zenodo.24099265](https://doi.org/10.5281/zenodo.24099265)) and is
downloaded here, with its checksum checked, by:

```sh
Rscript scripts/00_download_data.R
```

To fetch only this file:

```sh
curl -L -o data/CellChatDB.human_AAT.rda \
  "https://zenodo.org/records/24099265/files/CellChatDB.human_AAT.rda?download=1"
```

The figures don't need it: they use the saved CellChat object, which
contains a copy of the database. It's needed to rerun CellChat
(`Rscript scripts/03_figure3_cellchat.R --recompute`), or to use the extended
database with your own data:

```r
load("data/CellChatDB.human_AAT.rda")
cellchat <- CellChat::createCellChat(seurat_obj, group.by = "cell_type")
cellchat@DB <- db.AAT
```
