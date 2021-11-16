# CoCoCoNET (co-expression)
wget http://labshare.cshl.edu/shares/gillislab/resource/CoCoCoNet/networks/mouse_HC_AggNet.hdf5 -O mouse_cococonet_210420.hdf5

# BIOGRID (PPI)
wget https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.197/BIOGRID-ALL-4.4.197.tab2.zip
unzip BIOGRID-ALL-4.4.197.tab2.zip
rm BIOGRID-ALL-4.4.197.tab2.zip

# GO
wget http://current.geneontology.org/ontology/go-basic.obo -O go-basic_181219.obo

# GO annotations
wget http://www.informatics.jax.org/downloads/reports/gene_association.mgi.gz