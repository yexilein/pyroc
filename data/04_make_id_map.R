library(tidyverse)
library(org.Mm.eg.db)
library(rhdf5)


main = function() {
    make_mouse_id_map("mouse_ids_210421.csv")
}

make_mouse_id_map = function(filename) {
    ensembl_ids = h5read("mouse_coconet_210420.hdf5", "/col")
    # only keep 1:1 mappings
    symbols = mapIds(org.Mm.eg.db, ensembl_ids, "SYMBOL", "ENSEMBL", multiVals="asNA")
    symbols = enframe(symbols, name = "ensembl", value = "symbol")
    write_csv(symbols, filename)
}

if (sys.nframe() == 0) {
    main()
}
