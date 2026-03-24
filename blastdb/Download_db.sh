
# download swissprot
# mkdir -p ./blastdb/swissprot
# cd ./blastdb/swissprot

# update_blastdb.pl --showall | grep -i swiss
# update_blastdb.pl --decompress swissprot


# download pdbaa
mkdir -p ./blastdb/pdbaa

cd ./blastdb/pdbaa
update_blastdb.pl --showall | grep -i pdbaa
update_blastdb.pl --decompress pdbaa

# 
mkdir -p ./blastdb/pdb_cluster_src
cd ./blastdb/pdb_cluster_src

wget https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
wget https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-95.txt
wget https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-90.txt
wget https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-70.txt
wget https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-50.txt
gunzip -f pdb_seqres.txt.gz