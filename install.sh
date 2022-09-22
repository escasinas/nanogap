# # # # # # # # # # # # # # # # # # # # # # #
#                                           #
# To run this script, execute the command:  #
# source install.sh                         #
#                                           #
# # # # # # # # # # # # # # # # # # # # # # #

# create a conda environment called "nanogap" using the environment.yml file
conda env create -f environment.yml

# activate the conda environment
conda activate nanogap

# install the BLAST database
mkdir blastdb && cd blastdb && update_blastdb.pl 16S_ribosomal_RNA --decompress && cd ..

# deactivate the conda environment
conda deactivate

