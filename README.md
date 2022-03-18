# NcycFunGen
Perl script for construction of NcycFenGen database.

#Read HMMER Search out file and save to database, meanwhile clasify sequences as high, midle, low and homologous level with e-value and gene information
perl read_hmmer_out_2022.pl

#download CDS seuqences according to protein accession number from GenBank
perl download_CDS_Genbank_2022.pl

#import taxon information to database according to gi and taxid from GenBank
perl import_taxinfo_to_database_from_taxid_2022.pl

#sequence similarity calculation between every two aligned sequences
perl sequence_similarity_calculation_2022.pl
