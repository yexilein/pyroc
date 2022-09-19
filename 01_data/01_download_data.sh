# CoCoCoNET (co-expression)
wget http://labshare.cshl.edu/shares/gillislab/resource/CoCoCoNet/networks/mouse_HC_AggNet.hdf5 -O mouse_cococonet_210420.hdf5

# BIOGRID (PPI)
wget https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.197/BIOGRID-ALL-4.4.197.tab2.zip
unzip BIOGRID-ALL-4.4.197.tab2.zip
rm BIOGRID-ALL-4.4.197.tab2.zip

# DrugBank (DTI) -> replace with own login and password 
#curl -Lfv -o drugbank_5.1.9.zip -u LOGIN:PASSWORD https://go.drugbank.com/releases/5-1-9/downloads/all-full-database
# password -> UBemejPQ8t2!kEZ
curl -Lfv -o drugbank_5.1.9.zip -u fischer@cshl.edu https://go.drugbank.com/releases/5-1-9/downloads/all-full-database

# STITCH
#wget http://stitch.embl.de/download/chemical_chemical.links.detailed.v5.0.tsv.gz
wget http://stitch.embl.de/download/protein_chemical.links.detailed.v5.0/10090.protein_chemical.links.detailed.v5.0.tsv.gz
#wget http://stitch.embl.de/download/actions.v5.0/10090.actions.v5.0.tsv.gz
wget https://stringdb-static.org/download/protein.aliases.v11.5/10090.protein.aliases.v11.5.txt.gz

# GO
wget http://current.geneontology.org/ontology/go-basic.obo -O go-basic_181219.obo

# GO annotations
wget http://www.informatics.jax.org/downloads/reports/gene_association.mgi.gz

# UniprotKB annotatinos
#python download_interpro.py > interpro_220812.tsv
#wget https://ftp.ebi.ac.uk/pub/databases/interpro/protein2ipr.dat.gz
wget "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Cft_domain&format=tsv&query=%28taxonomy_id%3A10090%29%20AND%20%28reviewed%3Atrue%29%20AND%20%28proteins_with%3A24%29" -O uniprot_220812.tsv.gz