#!/usr/bin/env bash

########### State of Databases at July 18th 2011 ###############################
#
# Legend:
# X => database exists for the given organism and is normalized
# _ => database doesn't exist for the given organism
# ? => database exists for the given organism but either it can't be normalized due to bad formatting or it has to be downloaded manually (mirwalk predictive)
#
#                                              |_hsa_|_mmu_|_rno_|_gga_|_dre_|_dme_|_cel_|
# Experimental Tarbase_________________________|  X  |  X  |  _  |  _  |  _  |  _  |  _  |
# Experimental MirTarBase______________________|  X  |  X  |  X  |  X  |  X  |  X  |  X  |
# Experimental miRecords_______________________|  X  |  X  |  X  |  _  |  X  |  X  |  X  |
# Experimental mirWalk_________________________|  X  |  X  |  X  |  _  |  _  |  _  |  _  |
#                                              |_hsa_|_mmu_|_rno_|_gga_|_dre_|_dme_|_cel_|
# Predictive microrna.org conserved high score_|  X  |  X  |  X  |  _  |  _  |  X  |  X  |
# Predictive microcosm_________________________|  X  |  X  |  X  |  X  |  X  |  _  |  _  |
# Predictive TargetScan conserved______________|  X  |  X  |  _  |  _  |  _  |  ?  |  ?  |
# Predictive PITA______________________________|  X  |  X  |  _  |  _  |  _  |  X  |  X  | 
# Predictive DIANA-microT______________________|  X  |  X  |  _  |  _  |  _  |  X  |  ?  |
# Predictive EIMMo_____________________________|  X  |  X  |  X  |  _  |  X  |  X  |  X  |
# Predictive MirTarget (miRDB)_________________|  X  |  X  |  X  |  X  |  _  |  _  |  _  |
# Predictive Targetspy_________________________|  X  |  X  |  X  |  X  |  _  |  ?  |  _  |
# Predictive mirWalk___________________________|  X  |  ?  |  ?  |  _  |  _  |  _  |  _  |
################################################################################


########### INITIALIZATION #####################################################

# Attempt to use an undefined variable outputs error message and exits. 
set -o nounset 

sep=';'
database_file_extension='.csv'

# TODO: update this list
all_databases="all_calculations | all_auxiliary | all_experimental | all_predictive | combined_experimental | combined_precision | ensembl_gene_id_to_genbank_id | ensembl_transcript_id_to_genbank_id | experimental_mirtarbase | experimental_mirwalk | genbank_id_to_gene_name | mirecords | mirnas | predictive_eimmo | predictive_microcosm | predictive_microrna_org_conserved_high_score | predictive_microt | predictive_mirtarget | predictive_mirwalk | predictive_pita | predictive_targetscan_conserved | predictive_targetspy | refseq_id_to_genbank_id | experimental_tarbase | precision | plot_score_distribution | plot_precision | plot_ROC"

usage ()
{
     echo    
     echo "USAGE: $(basename $0) ORGANISM DATABASE"
     echo
     echo "ORGANISM can be: hsa | mmu | rno | gga | dre | dme | cel"
     echo
     echo "DATABASE can be:"
     echo $all_databases
     echo
     echo "Generates the appropriate database for the given organism. This script needs Bash 4.x"
     echo
}

if [[ $# < 2 ]]; then
    echo "Wrong arguments: Please specify a valid action and organism"
    usage
    exit
fi

organism=$1
query_database=$2
if [[ $# > 2 ]]; then
    specific_database_name=$3 # optional, only for some CALCULATIONS
else
    specific_database_name=
fi

################################################################################


########### DIRECTORY STRUCTURE ################################################
# +-/
#   |-mirna_scroll            
#   +-plots/                  <--plots_path 
#   | +-hsa/                  <--plots_organism_path
#   |   +-plot.png            
#   +-lib/                    <--lib_path 
#   | +-translator.awk        <--translator_script
#   +-raw/                    <--raw_path
#   | +-primary/              <--raw_primary_path
#   +-normalized/             <--normalized_path
#     +-hsa/                  <--normalized_organism_path  
#     | +-primary/            <--normalized_organism_primary_path  
#     | | +-experimental_...  <--experimental_database_paths
#     | | +-predictive_...    <--predictive_database_paths
#     | +-calculations/       <--normalized_organism_calculations_path
#     +-mmu/                  
#     | +-primary/
#     | | +-experimental_... 
#     | | +-predictive_...   
#     | +-calculations/
#     +-...

declare -A raw_paths
declare -A auxiliary_paths
declare -A experimental_paths
declare -A predictive_paths

declare -A precision_paths
declare -A roc_curve_paths

raw_path="raw/"
normalized_path="normalized/"
organism_path="${organism}/"
primary_path="primary/"
calculations_path="calculations/"

raw_primary_path="${raw_path}${primary_path}"
mkdir -pv $raw_primary_path

normalized_organism_path="${normalized_path}${organism_path}"
normalized_organism_primary_path="${normalized_organism_path}${primary_path}"
normalized_organism_calculations_path="${normalized_organism_path}${calculations_path}"
mkdir -pv $normalized_organism_primary_path $normalized_organism_calculations_path

cutoffs_file="${normalized_organism_path}cutoffs.csv"
touch $cutoffs_file

lib_path="lib/"
plots_path="plots/"
plots_organism_path="${plots_path}${organism_path}"
mkdir -pv $lib_path $plots_organism_path

combine_predictive_databases_script="${lib_path}combine_predictive_databases.awk"
cutoff_script="${lib_path}cutoff.awk"
precision_script="${lib_path}precision.awk"
roc_curve_script="${lib_path}roc_curve.awk"
score_normalization_script="${lib_path}score_normalization.awk"
translator_script="${lib_path}translator.awk"
unique_interactions_script="${lib_path}unique_interactions.awk"

################################################################################


########### AUXILIARY DATABASES ################################################

hsa_normalize_genbank_id_to_gene_name(){
    # filter human genes (the file also has Neanderthal genes)
    awk -v OFS=$sep '$1 == 9606 {print $2, $3}'
}
mmu_normalize_genbank_id_to_gene_name(){
    # filter mouse genes (the file also has 10091, 10092, 39442, 57486)
    awk -v OFS=$sep '$1 == 10090 {print $2, $3}'
}
rno_normalize_genbank_id_to_gene_name(){
    # filter rat genes
    awk -v OFS=$sep '$1 == 10116 {print $2, $3}'
}
gga_normalize_genbank_id_to_gene_name(){
    # filter chicken genes (the file also has 208524, 208525, 208526)
    awk -v OFS=$sep '$1 == 9031 {print $2, $3}'
}
dre_normalize_genbank_id_to_gene_name(){
    # filter fish genes 
    awk -v OFS=$sep '$1 == 7955 {print $2, $3}'
}
dme_normalize_genbank_id_to_gene_name(){
    # filter fly genes 
    awk -v OFS=$sep '$1 == 7227 {print $2, $3}'
}
cel_normalize_genbank_id_to_gene_name(){
    # filter worm genes 
    awk -v OFS=$sep '$1 == 6239 {print $2, $3}'
}

current_database_name="genbank_id_to_gene_name"; auxiliary_paths["$current_database_name"]="${normalized_organism_path}${current_database_name}${database_file_extension}"
normalize_genbank_id_to_gene_name(){
    local current_database_name=$(echo ${FUNCNAME#*normalize_})
    case $organism in
        hsa|mmu|rno|gga|dre|dme|cel)

            raw_paths["hsa:${current_database_name}"]="${raw_path}Homo_sapiens.gene_info" # downloaded from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
            raw_paths["mmu:${current_database_name}"]="${raw_path}Mus_musculus.gene_info" # downloaded from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz
            raw_paths["rno:${current_database_name}"]="${raw_path}Rattus_norvegicus.gene_info" # downloaded from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/
            raw_paths["gga:${current_database_name}"]="${raw_path}Gallus_gallus.gene_info" # downloaded from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/
            raw_paths["dre:${current_database_name}"]="${raw_path}Danio_rerio.gene_info" # downloaded from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/
            raw_paths["dme:${current_database_name}"]="${raw_path}Drosophila_melanogaster.gene_info" # downloaded from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/
            raw_paths["cel:${current_database_name}"]="${raw_path}Caenorhabditis_elegans.gene_info" # downloaded from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INF
            
            local input_database_path=${raw_paths["${organism}:${current_database_name}"]}
            local output_database_path=${auxiliary_paths["$current_database_name"]}
            
            if [ -f $input_database_path ] ; then    
                echo "Normalizing $organism $current_database_name database"
                cat $input_database_path | 
                ${organism}_normalize_${current_database_name} |
                sort -u > $output_database_path

                echo_summary $output_database_path
            else
                echo "Path $input_database_path doesn't exist"
            fi
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}


hsa_normalize_refseq_id_to_genbank_id(){
    # skip header, filter organism, refseq_id is not empty, delete version numbers, print refseq_id genbank_id, 
    awk -v OFS=$sep 'NR > 1 && $1 == 9606 && $4 != "-" {split($4,refseq_id,".");print refseq_id[1], $2}' 
}
mmu_normalize_refseq_id_to_genbank_id(){
    # skip header, filter organism, refseq_id is not empty, delete version numbers, print refseq_id genbank_id, filter uniques
    awk -v OFS=$sep 'NR > 1 && $1 == 10090 && $4 != "-" {split($4,refseq_id,".");print refseq_id[1], $2}'
}
rno_normalize_refseq_id_to_genbank_id(){
    # skip header, filter organism, refseq_id is not empty, delete version numbers, print refseq_id genbank_id, filter uniques
    awk -v OFS=$sep 'NR > 1 && $1 == 10116 && $4 != "-" {split($4,refseq_id,".");print refseq_id[1], $2}'
}
gga_normalize_refseq_id_to_genbank_id(){
    # skip header, filter organism, refseq_id is not empty, delete version numbers, print refseq_id genbank_id, filter uniques
    awk -v OFS=$sep 'NR > 1 && $1 == 9031 && $4 != "-" {split($4,refseq_id,".");print refseq_id[1], $2}'
}
dre_normalize_refseq_id_to_genbank_id(){
    # skip header, filter organism, refseq_id is not empty, delete version numbers, print refseq_id genbank_id, filter uniques
    awk -v OFS=$sep 'NR > 1 && $1 == 7955 && $4 != "-" {split($4,refseq_id,".");print refseq_id[1], $2}'
}
dme_normalize_refseq_id_to_genbank_id(){
    # skip header, filter organism, refseq_id is not empty, delete version numbers, print refseq_id genbank_id, filter uniques
    awk -v OFS=$sep 'NR > 1 && $1 == 7227 && $4 != "-" {split($4,refseq_id,".");print refseq_id[1], $2}'
}
cel_normalize_refseq_id_to_genbank_id(){
    # skip header, filter organism, refseq_id is not empty, delete version numbers, print refseq_id genbank_id, filter uniques
    awk -v OFS=$sep 'NR > 1 && $1 == 6239 && $4 != "-" {split($4,refseq_id,".");print refseq_id[1], $2}'
}

current_database_name="refseq_id_to_genbank_id"; auxiliary_paths["$current_database_name"]="${normalized_organism_path}${current_database_name}${database_file_extension}"
normalize_refseq_id_to_genbank_id(){
    local current_database_name=$(echo ${FUNCNAME#*normalize_})
    case $organism in
        hsa|mmu|rno|gga|dre|dme|cel)
            
            raw_paths["${current_database_name}"]="${raw_path}gene2refseq" # downloaded from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq

            local input_database_path=${raw_paths["${current_database_name}"]}
            local output_database_path=${auxiliary_paths["$current_database_name"]}
            
            if [ -f $input_database_path ] ; then    
                echo "Normalizing $organism $current_database_name database"
                cat $input_database_path | 
                ${organism}_normalize_${current_database_name} |
                sort -u > $output_database_path

                echo_summary $output_database_path
            else
                echo "Path $input_database_path doesn't exist"
            fi
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}


hsa_normalize_ensembl_transcript_id_to_genbank_id(){
    # skip header, filter organism, choose ensembl_transcript_id and genbank_id columns, 
    awk -v OFS=$sep 'NR > 1 && $1 == 9606 && $5 != "-" {print $5, $2}'
}
mmu_normalize_ensembl_transcript_id_to_genbank_id(){
    # skip header, filter organism, choose ensembl_transcript_id and genbank_id columns, 
    awk -v OFS=$sep 'NR > 1 && $1 == 10090 && $5 != "-" {print $5, $2}'
}
rno_normalize_ensembl_transcript_id_to_genbank_id(){
    # skip header, filter organism, choose ensembl_transcript_id and genbank_id columns, 
    awk -v OFS=$sep 'NR > 1 && $1 == 10116 && $5 != "-" {print $5, $2}'
}
gga_normalize_ensembl_transcript_id_to_genbank_id(){
    # skip header, filter organism, choose ensembl_transcript_id and genbank_id columns, 
    awk -v OFS=$sep 'NR > 1 && $1 == 9031 && $5 != "-" {print $5, $2}'
}
dre_normalize_ensembl_transcript_id_to_genbank_id(){
    # skip header, filter organism, choose ensembl_transcript_id and genbank_id columns, 
    awk -v OFS=$sep 'NR > 1 && $1 == 7955 && $5 != "-" {print $5, $2}'
}
dme_normalize_ensembl_transcript_id_to_genbank_id(){
    # skip header, filter organism, choose ensembl_transcript_id and genbank_id columns, 
    awk -v OFS=$sep 'NR > 1 && $1 == 7227 && $5 != "-" {print $5, $2}'
}

current_database_name="ensembl_transcript_id_to_genbank_id"; auxiliary_paths["$current_database_name"]="${normalized_organism_path}${current_database_name}${database_file_extension}"
normalize_ensembl_transcript_id_to_genbank_id(){
    local current_database_name=$(echo ${FUNCNAME#*normalize_})
    case $organism in
        hsa|mmu|rno|gga|dre|dme)

            raw_paths["${current_database_name}"]="${raw_path}gene2ensembl" # downloaded from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl

            local input_database_path=${raw_paths["${current_database_name}"]}
            local output_database_path=${auxiliary_paths["$current_database_name"]}

            if [ -f $input_database_path ] ; then    
                echo "Normalizing $organism $current_database_name database"
                cat $input_database_path | 
                ${organism}_normalize_${current_database_name} |
                sort -u > $output_database_path

                echo_summary $output_database_path
            else
                echo "Path $input_database_path doesn't exist"
            fi
            
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}


hsa_normalize_ensembl_gene_id_to_genbank_id(){
    # skip header, filter organism, choose ensembl_gene_id and genbank_id columns, 
    awk -v OFS=$sep 'NR > 1 && $1 == 9606 {print $3, $2}'
}
mmu_normalize_ensembl_gene_id_to_genbank_id(){
    # skip header, filter organism, choose ensembl_gene_id and genbank_id columns, 
    awk -v OFS=$sep 'NR > 1 && $1 == 10090 {print $3, $2}'
}
rno_normalize_ensembl_gene_id_to_genbank_id(){
    # skip header, filter organism, choose ensembl_gene_id and genbank_id columns, 
    awk -v OFS=$sep 'NR > 1 && $1 == 10116 {print $3, $2}'
}
gga_normalize_ensembl_gene_id_to_genbank_id(){
    # skip header, filter organism, choose ensembl_gene_id and genbank_id columns, 
    awk -v OFS=$sep 'NR > 1 && $1 == 9031 {print $3, $2}'
}
dre_normalize_ensembl_gene_id_to_genbank_id(){
    # skip header, filter organism, choose ensembl_gene_id and genbank_id columns, 
    awk -v OFS=$sep 'NR > 1 && $1 == 7955 {print $3, $2}'
}
dme_normalize_ensembl_gene_id_to_genbank_id(){
    # skip header, filter organism, choose ensembl_gene_id and genbank_id columns, 
    awk -v OFS=$sep 'NR > 1 && $1 == 7227 {print $3, $2}'
}

current_database_name="ensembl_gene_id_to_genbank_id"; auxiliary_paths["$current_database_name"]="${normalized_organism_path}${current_database_name}${database_file_extension}"
normalize_ensembl_gene_id_to_genbank_id(){
    local current_database_name=$(echo ${FUNCNAME#*normalize_})
    case $organism in
        hsa|mmu|rno|gga|dre|dme)

            raw_paths["${current_database_name}"]="${raw_path}gene2ensembl" # downloaded from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl

            local input_database_path=${raw_paths["${current_database_name}"]}
            local output_database_path=${auxiliary_paths["$current_database_name"]}
            
            if [ -f $input_database_path ] ; then    
                echo "Normalizing $organism $current_database_name database"
                cat $input_database_path | 
                ${organism}_normalize_${current_database_name} |
                sort -u > $output_database_path
                
                echo_summary $output_database_path
            else
                echo "Path $input_database_path doesn't exist"
            fi
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}

# This file is used to download experimental_mirwalk
hsa_normalize_mirnas(){
    # filter organism, remove ">" from miRNA fasta headers
    awk '$1 ~ /hsa/{sub(/>/,"",$1); print $1}'
}
mmu_normalize_mirnas(){
    # filter organism, remove ">" from miRNA fasta headers
    awk '$1 ~ /mmu/{sub(/>/,"",$1); print $1}'
}
rno_normalize_mirnas(){
    # filter organism, remove ">" from miRNA fasta headers
    awk '$1 ~ /rno/{sub(/>/,"",$1); print $1}'
}
gga_normalize_mirnas(){
    # filter organism, remove ">" from miRNA fasta headers
    awk '$1 ~ /gga/{sub(/>/,"",$1); print $1}'
}
dre_normalize_mirnas(){
    # filter organism, remove ">" from miRNA fasta headers
    awk '$1 ~ /dre/{sub(/>/,"",$1); print $1}'
}
dme_normalize_mirnas(){
    # filter organism, remove ">" from miRNA fasta headers
    awk '$1 ~ /dme/{sub(/>/,"",$1); print $1}'
}
cel_normalize_mirnas(){
    # filter organism, remove ">" from miRNA fasta headers
    awk '$1 ~ /cel/{sub(/>/,"",$1); print $1}'
}

current_database_name="mirnas"; auxiliary_paths["$current_database_name"]="${normalized_organism_path}${current_database_name}${database_file_extension}"
normalize_mirnas(){
    local current_database_name=$(echo ${FUNCNAME#*normalize_})
    case $organism in
        hsa|mmu|rno|gga|dre|dme|cel)

            raw_paths["${current_database_name}"]="${raw_path}mirnas.fa" # downloaded from ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz

            local input_database_path=${raw_paths["${current_database_name}"]}
            local output_database_path=${auxiliary_paths["$current_database_name"]}
            
            if [ -f $input_database_path ] ; then    
                echo "Normalizing $organism $current_database_name database"
                cat $input_database_path | 
                ${organism}_normalize_${current_database_name} |
                sort -u > $output_database_path
                
                echo_summary $output_database_path
            else
                echo "Path $input_database_path doesn't exist"
            fi
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}


################################################################################


########### PRIMARY DATABASES (EXPERIMENTAL) ###################################
# FORMAT: mirna_name;genbank_id;gene_name;experimental_evidence;database_name

hsa_normalize_experimental_tarbase(){
    # remove header, filter organism rows, normalize mirna and choose ensembl_gene and experiment columns, 
    awk -F";" -v OFS=$sep 'NR > 1 && $5 == "Human" {print "hsa-"$6, $10, $3"_"$4}'
}
mmu_normalize_experimental_tarbase(){
    # remove header, filter organism rows, normalize mirna and choose ensembl_gene and experiment columns, 
    awk -F";" -v OFS=$sep 'NR > 1 && $5 == "Mouse" {print "mmu-"$6, $10, $3"_"$4}'
}

current_database_name="experimental_tarbase"; experimental_paths["$current_database_name"]="${normalized_organism_primary_path}${current_database_name}${database_file_extension}"
normalize_experimental_tarbase(){
    local current_database_name=$(echo ${FUNCNAME#*normalize_})
    case $organism in
        hsa|mmu)
            raw_paths["${current_database_name}"]="${raw_primary_path}TarBase_V5.0.csv" # downloaded from http://diana.cslab.ece.ntua.gr/data/public/TarBase_V5.0.rar

            local input_database_path=${raw_paths["${current_database_name}"]}
            local output_database_path=${experimental_paths["$current_database_name"]}
            
            if [ -f $input_database_path ] ; then    
                echo "Normalizing $organism $current_database_name database"
                cat $input_database_path | 
                ${organism}_normalize_${current_database_name} |
                ensembl_gene_id_to_genbank_id |
                filter_unique_interactions | 
                add_gene_names | 
                order_columns_add_database_name $current_database_name | 
                add_experimental_score |
                sort > $output_database_path
                
                echo_summary $output_database_path
            else
                echo "Path $input_database_path doesn't exist"
            fi
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}


hsa_normalize_experimental_mirtarbase(){
    # remove header, filter organism, choose mirna, genbank and experiment columns, 
    awk -F"\t" -v OFS=$sep 'NR > 1 && $2 ~ /hsa/ {print $2, $5, $7}'
}
mmu_normalize_experimental_mirtarbase(){
    # remove header, filter organism, choose mirna, genbank and experiment columns, 
    awk -F"\t" -v OFS=$sep 'NR > 1 && $2 ~ /mmu/ {print $2, $5, $7}'
}
rno_normalize_experimental_mirtarbase(){
    # remove header, filter organism, choose mirna, genbank and experiment columns, 
    awk -F"\t" -v OFS=$sep 'NR > 1 && $2 ~ /rno/ {print $2, $5, $7}'
}
gga_normalize_experimental_mirtarbase(){
    # remove header, filter organism, choose mirna, genbank and experiment columns, 
    awk -F"\t" -v OFS=$sep 'NR > 1 && $2 ~ /gga/ {print $2, $5, $7}'
}
dre_normalize_experimental_mirtarbase(){
    # remove header, filter organism, choose mirna, genbank and experiment columns, 
    awk -F"\t" -v OFS=$sep 'NR > 1 && $2 ~ /dre/ {print $2, $5, $7}'
}
dme_normalize_experimental_mirtarbase(){
    # remove header, filter organism, choose mirna, genbank and experiment columns, 
    awk -F"\t" -v OFS=$sep 'NR > 1 && $2 ~ /dme/ {print $2, $5, $7}'
}
cel_normalize_experimental_mirtarbase(){
    # remove header, filter organism, choose mirna, genbank and experiment columns, 
    awk -F"\t" -v OFS=$sep 'NR > 1 && $2 ~ /cel/ {print $2, $5, $7}'
}

current_database_name="experimental_mirtarbase"; experimental_paths["$current_database_name"]="${normalized_organism_primary_path}${current_database_name}${database_file_extension}"
normalize_experimental_mirtarbase(){
    local current_database_name=$(echo ${FUNCNAME#*normalize_})
    case $organism in
        hsa|mmu|rno|gga|dre|dme|cel)
            
            raw_paths["hsa:${current_database_name}"]="${raw_primary_path}mirtarbase_hsa.csv" # downloaded from http://mirtarbase.mbc.nctu.edu.tw/cache/download/hsa.xls
            raw_paths["mmu:${current_database_name}"]="${raw_primary_path}mirtarbase_mmu.csv" # downloaded from http://mirtarbase.mbc.nctu.edu.tw/cache/download/mmu.xls
            raw_paths["rno:${current_database_name}"]="${raw_primary_path}mirtarbase_rno.csv" # downloaded from http://mirtarbase.mbc.nctu.edu.tw/cache/download/rno.xls
            raw_paths["gga:${current_database_name}"]="${raw_primary_path}mirtarbase_gga.csv" # downloaded from http://mirtarbase.mbc.nctu.edu.tw/cache/download/gga.xls
            raw_paths["dre:${current_database_name}"]="${raw_primary_path}mirtarbase_dre.csv" # downloaded from http://mirtarbase.mbc.nctu.edu.tw/cache/download/dre.xls
            raw_paths["dme:${current_database_name}"]="${raw_primary_path}mirtarbase_dme.csv" # downloaded from http://mirtarbase.mbc.nctu.edu.tw/cache/download/dme.xls
            raw_paths["cel:${current_database_name}"]="${raw_primary_path}mirtarbase_cel.csv" # downloaded from http://mirtarbase.mbc.nctu.edu.tw/cache/download/cel.xls

            local input_database_path=${raw_paths["${organism}:${current_database_name}"]}
            local output_database_path=${experimental_paths["$current_database_name"]}
            
            if [ -f $input_database_path ] ; then    
                echo "Normalizing $organism $current_database_name database"
                cat $input_database_path | 
                ${organism}_normalize_${current_database_name} |
                filter_unique_interactions | 
                add_gene_names | 
                order_columns_add_database_name $current_database_name | 
                add_experimental_score |
                sort > $output_database_path
                
                echo_summary $output_database_path
            else
                echo "Path $input_database_path doesn't exist"
            fi
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}


hsa_normalize_experimental_mirecords(){
    # remove header, filter organism rows, choose non-empty refseq genes, mirnas and experiment columns, fix "has" => "hsa" typo, remove brackets from mirnas, 
    # remove version number from refseq_id, add species code to mirnas that don't have it 
    awk -F";" -v OFS=$sep 'NR > 1 && $2 ~ "Homo sapiens" && $7 ~ "Homo sapiens" && $5 != "" && $5 !~ "unknown" {sub(/has/,"hsa",$8); gsub(/\[|\]/,"", $8); print $8, $5, $12}' |
    awk -F$sep -v OFS=$sep '{ split($2, refseq,"."); $2=refseq[1]; if($1 !~ /hsa/) $1="hsa-"$1; print }'  
}
mmu_normalize_experimental_mirecords(){
    # remove header, filter organism rows, choose non-empty refseq genes, mirnas and experiment columns, remove brackets from mirnas, 
    # remove version number from refseq_id, add species code to mirnas that don't have it
    awk -F";" -v OFS=$sep 'NR > 1 && $2 ~ "Mus musculus" && $7 ~ "Mus musculus" && $5 != "" && $5 !~ "unknown" {gsub(/\[|\]/,"", $8); print $8, $5, $12}' |
    awk -F$sep -v OFS=$sep '{ split($2, refseq,"."); $2=refseq[1]; if($1 !~ /mmu/) $1="mmu-"$1; print }'
}
rno_normalize_experimental_mirecords(){
    # remove header, filter organism rows, choose non-empty refseq genes, mirnas and experiment columns, remove brackets from mirnas, 
    # remove version number from refseq_id, add species code to mirnas that don't have it
    awk -F";" -v OFS=$sep 'NR > 1 && $2 ~ "Rattus norvegicus" && $7 ~ "Rattus norvegicus" && $5 != "" && $5 !~ "unknown" {gsub(/\[|\]/,"", $8); print $8, $5, $12}' |
    awk -F$sep -v OFS=$sep '{ split($2, refseq,"."); $2=refseq[1]; if($1 !~ /rno/) $1="rno-"$1; print }'
}
dre_normalize_experimental_mirecords(){
    # remove header, filter organism rows, choose non-empty refseq genes, mirnas and experiment columns, remove brackets from mirnas, 
    # remove version number from refseq_id, add species code to mirnas that don't have it
    awk -F";" -v OFS=$sep 'NR > 1 && $2 ~ "Danio rerio" && $7 ~ "Danio rerio" && $5 != "" && $5 !~ "unknown" {gsub(/\[|\]/,"", $8); print $8, $5, $12}' |
    awk -F$sep -v OFS=$sep '{ split($2, refseq,"."); $2=refseq[1]; if($1 !~ /dre/) $1="dre-"$1; print }'
}
dme_normalize_experimental_mirecords(){
    # remove header, filter organism rows, choose non-empty refseq genes, mirnas and experiment columns, remove brackets from mirnas, 
    # remove version number from refseq_id, add species code to mirnas that don't have it
    awk -F";" -v OFS=$sep 'NR > 1 && $2 ~ "Drosophila melanogaster" && $7 ~ "Drosophila melanogaster" && $5 != "" && $5 !~ "unknown" {gsub(/\[|\]/,"", $8); print $8, $5, $12}' |
    awk -F$sep -v OFS=$sep '{ split($2, refseq,"."); $2=refseq[1]; if($1 !~ /dme/) $1="dme-"$1; print }'
}
cel_normalize_experimental_mirecords(){
    # remove header, filter organism rows, choose non-empty refseq genes, mirnas and experiment columns, remove brackets from mirnas, 
    # remove version number from refseq_id, add species code to mirnas that don't have it
    awk -F";" -v OFS=$sep 'NR > 1 && $2 ~ "Caenorhabditis elegans" && $7 ~ "Caenorhabditis elegans" && $5 != "" && $5 !~ "unknown" {gsub(/\[|\]/,"", $8); print $8, $5, $12}' |
    awk -F$sep -v OFS=$sep '{ split($2, refseq,"."); $2=refseq[1]; if($1 !~ /cel/) $1="cel-"$1; print }'
}

current_database_name="experimental_mirecords"; experimental_paths["$current_database_name"]="${normalized_organism_primary_path}${current_database_name}${database_file_extension}"
normalize_experimental_mirecords(){
    local current_database_name=$(echo ${FUNCNAME#*normalize_})
    case $organism in
        hsa|mmu|rno|dre|dme|cel)
            
            raw_paths["${current_database_name}"]="${raw_primary_path}miRecords_version3.csv" # downloaded from http://mirecords.biolead.org/download_data.php?v=2


            local input_database_path=${raw_paths["${current_database_name}"]}
            local output_database_path=${experimental_paths["$current_database_name"]}
                        
            if [ -f $input_database_path ] ; then    
                echo "Normalizing $organism $current_database_name database"
                cat $input_database_path | 
                ${organism}_normalize_${current_database_name} | 
                refseq_id_to_genbank_id |  
                filter_unique_interactions | 
                add_gene_names | 
                order_columns_add_database_name $current_database_name | 
                add_experimental_score |
                sort > $output_database_path
                
                echo_summary $output_database_path
            else
                echo "Path $input_database_path doesn't exist"
            fi
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}


choose_experimental_mirwalk_columns(){
    # choose mirnas, genbank_ids and experiment_id columns,
    awk -F"\t" -v OFS=$sep '{print $1, $5, "_"}' 
}

 
current_database_name="experimental_mirwalk"; experimental_paths["$current_database_name"]="${normalized_organism_primary_path}${current_database_name}${database_file_extension}"
normalize_experimental_mirwalk(){
    local current_database_name=$(echo ${FUNCNAME#*normalize_})
    case $organism in
        hsa|mmu|rno)
            
            raw_paths["hsa:${current_database_name}"]="${raw_primary_path}mirwalk_hsa.csv" # downloaded from http://www.ma.uni-heidelberg.de/apps/zmf/mirwalk/mirnatargetpub.html using ${normalized_organism_path}mirnas.txt as input
            raw_paths["mmu:${current_database_name}"]="${raw_primary_path}mirwalk_mmu.csv" # downloaded from http://www.ma.uni-heidelberg.de/apps/zmf/mirwalk/mirnatargetpub.html using ${normalized_organism_path}mirnas.txt as input
            raw_paths["rno:${current_database_name}"]="${raw_primary_path}mirwalk_rno.csv" # downloaded from http://www.ma.uni-heidelberg.de/apps/zmf/mirwalk/mirnatargetpub.html using ${normalized_organism_path}mirnas.txt as input

            local input_database_path=${raw_paths["${organism}:${current_database_name}"]}
            local output_database_path=${experimental_paths["$current_database_name"]}
            
            if [ -f $input_database_path ] ; then    
                echo "Normalizing $organism $current_database_name database"
                cat $input_database_path | 
                choose_${current_database_name}_columns | 
                filter_unique_interactions | 
                add_gene_names | 
                order_columns_add_database_name $current_database_name | 
                add_experimental_score |
                sort > $output_database_path
            
                echo_summary $output_database_path
            else
                echo "Path $input_database_path doesn't exist"
            fi
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}

################################################################################


########### PRIMARY DATABASES (PREDICTIVE) #####################################
# FORMAT: mirna_name;genbank_id;gene_name;score;database_name

choose_predictive_microrna_org_conserved_high_score_columns(){
    # remove header, choose mirna, genbank_id and score columns,
    awk -F"\t" -v OFS=$sep 'NR > 1 {print $2, $3, $19}'
}

current_database_name="predictive_microrna_org_conserved_high_score"; predictive_paths["$current_database_name"]="${normalized_organism_primary_path}${current_database_name}${database_file_extension}"
normalize_predictive_microrna_org_conserved_high_score(){
    local current_database_name=$(echo ${FUNCNAME#*normalize_})
    case $organism in
        hsa|mmu|rno|dme|cel)
            
            raw_paths["hsa:${current_database_name}"]="${raw_primary_path}org.microrna.human_predictions_S_C_aug2010.txt" # downloaded from http://cbio.mskcc.org/microrna_data/human_predictions_S_C_aug2010.txt.gz
            raw_paths["mmu:${current_database_name}"]="${raw_primary_path}org.microrna.mouse_predictions_S_C_aug2010.txt" # downloaded from http://cbio.mskcc.org/microrna_data/mouse_predictions_S_C_aug2010.txt.gz
            raw_paths["rno:${current_database_name}"]="${raw_primary_path}org.microrna.rat_predictions_S_C_aug2010.txt" # downloaded from http://cbio.mskcc.org/microrna_data/rat_predictions_S_C_aug2010.txt.gz
            raw_paths["dme:${current_database_name}"]="${raw_primary_path}org.microrna.fruitfly_predictions_S_C_aug2010.txt" # downloaded from http://cbio.mskcc.org/microrna_data/fruitfly_predictions_S_C_aug2010.txt.gz
            raw_paths["cel:${current_database_name}"]="${raw_primary_path}org.microrna.nematode_predictions_S_C_aug2010.txt" # downloaded from http://cbio.mskcc.org/microrna_data/nematode_predictions_S_C_aug2010.txt.gz

            local input_database_path=${raw_paths["${organism}:${current_database_name}"]}
            local output_database_path=${predictive_paths["$current_database_name"]}

            local raw_scores_path="raw_scores.tmp"
            local raw_scores_cutoff_path="raw_scores_cutoff.tmp"
            
            if [ -f $input_database_path ] ; then    
                echo "Normalizing $organism $current_database_name database"
                cat $input_database_path | 
                choose_${current_database_name}_columns | 
                filter_unique_interactions | 
                add_gene_names |  
                order_columns_add_database_name $current_database_name | 
                sort -t$sep -gk4 > $raw_scores_path

                cat $raw_scores_path |
                set_size_cutoff > $raw_scores_cutoff_path
                
                normalize_scores $raw_scores_cutoff_path "invert" |
                sort -t$sep -grk4 > $output_database_path

                rm $raw_scores_path
                
                echo_summary $output_database_path
            else
                echo "Path $input_database_path doesn't exist"
            fi
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}


hsa_normalize_predictive_microcosm(){
    # remove header, filter organism rows, choose mirna, ensembl_gene_id and score columns, 
    awk -F"\t" -v OFS=$sep 'NR>1 && $2 ~ /hsa/ {print $2, $12, $11}'
}
mmu_normalize_predictive_microcosm(){
    # remove header, filter organism rows, choose mirna, ensembl_gene_id and score columns, 
    awk -F"\t" -v OFS=$sep 'NR>1 && $2 ~ /mmu/ {print $2, $12, $11}'
}
rno_normalize_predictive_microcosm(){
    # remove header, filter organism rows, choose mirna, ensembl_gene_id and score columns, 
    awk -F"\t" -v OFS=$sep 'NR>1 && $2 ~ /rno/ {print $2, $12, $11}'
}
gga_normalize_predictive_microcosm(){
    # remove header, filter organism rows, choose mirna, ensembl_gene_id and score columns, 
    awk -F"\t" -v OFS=$sep 'NR>1 && $2 ~ /gga/ {print $2, $12, $11}'
}
dre_normalize_predictive_microcosm(){
    # remove header, filter organism rows, choose mirna, ensembl_gene_id and score columns, 
    awk -F"\t" -v OFS=$sep 'NR>1 && $2 ~ /dre/ {print $2, $12, $11}'
}

current_database_name="predictive_microcosm"; predictive_paths["$current_database_name"]="${normalized_organism_primary_path}${current_database_name}${database_file_extension}"
normalize_predictive_microcosm(){
    local current_database_name=$(echo ${FUNCNAME#*normalize_})
    case $organism in
        hsa|mmu|rno|gga|dre)

            raw_paths["hsa:${current_database_name}"]="${raw_primary_path}microcosm.v5.txt.homo_sapiens" # downloaded from ftp://ftp.ebi.ac.uk/pub/databases/microcosm/v5/arch.v5.txt.homo_sapiens.zip
            raw_paths["mmu:${current_database_name}"]="${raw_primary_path}microcosm.v5.txt.mus_musculus" # downloaded from ftp://ftp.ebi.ac.uk/pub/databases/microcosm/v5/arch.v5.txt.mus_musculus.zip
            raw_paths["rno:${current_database_name}"]="${raw_primary_path}microcosm.v5.txt.rattus_norvegicus" # downloaded from ftp://ftp.ebi.ac.uk/pub/databases/microcosm/v5/arch.v5.txt..zip
            raw_paths["gga:${current_database_name}"]="${raw_primary_path}microcosm.v5.txt.gallus_gallus" # downloaded from ftp://ftp.ebi.ac.uk/pub/databases/microcosm/v5/arch.v5.txt..zip
            raw_paths["dre:${current_database_name}"]="${raw_primary_path}microcosm.v5.txt.danio_rerio" # downloaded from ftp://ftp.ebi.ac.uk/pub/databases/microcosm/v5/arch.v5.txt..zip
            raw_paths["dme:${current_database_name}"]="${raw_primary_path}microcosm.v5.txt.drosophila_melanogaster" # downloaded from ftp://ftp.ebi.ac.uk/pub/databases/microcosm/v5/arch.v5.txt..zip
            raw_paths["cel:${current_database_name}"]="${raw_primary_path}microcosm.v5.txt.caenorhabditis_elegans" # downloaded from ftp://ftp.ebi.ac.uk/pub/databases/microcosm/v5/arch.v5.txt..zip

            local input_database_path=${raw_paths["${organism}:${current_database_name}"]}
            local output_database_path=${predictive_paths["$current_database_name"]}

            local raw_scores_path="raw_scores.tmp"
            local raw_scores_cutoff_path="raw_scores_cutoff.tmp"
            
            if [ -f $input_database_path ] ; then    
                echo "Normalizing $organism $current_database_name database"
                cat $input_database_path | 
                ${organism}_normalize_${current_database_name} | 
                ensembl_transcript_id_to_genbank_id  | 
                filter_unique_interactions | 
                add_gene_names | 
                order_columns_add_database_name $current_database_name | 
                sort -t$sep -gk4 > $raw_scores_path

                cat $raw_scores_path | 
                set_size_cutoff > $raw_scores_cutoff_path
                
                normalize_scores $raw_scores_cutoff_path "invert" |
                sort -t$sep -grk4 > $output_database_path

                rm $raw_scores_path
                
                echo_summary $output_database_path
            else
                echo "Path $input_database_path doesn't exist"
            fi
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}

hsa_normalize_predictive_targetscan_conserved(){
    # remove header, filter organism rows, choose mirna, genbank and score columns, 
    awk -F"\t" -v OFS=$sep 'NR > 1 && $3=="9606" && $8!="NULL" && $9!="NULL" && $10!="NULL" && $11!="NULL" {print $4, $1, $8 + $9 + $10 + $11}'
}

# TODO: sum all scores for mmu?
mmu_normalize_predictive_targetscan_conserved(){
    # remove header, filter organism rows, choose mirna, genbank and score columns, 
    awk -F"\t" -v OFS=$sep 'NR > 1 && $3=="10090" && $11!="NULL" {print $4, $1, $11}' 
}

current_database_name="predictive_targetscan_conserved"; predictive_paths["$current_database_name"]="${normalized_organism_primary_path}${current_database_name}${database_file_extension}"
normalize_predictive_targetscan_conserved(){
    local current_database_name=$(echo ${FUNCNAME#*normalize_})
    case $organism in
        hsa|mmu)
            
            raw_paths["hsa:${current_database_name}"]="${raw_primary_path}targetscan_Conserved_Site_Context_Scores_hsa.txt" # downloaded from http://www.targetscan.org//vert_50//vert_50_data_download/Conserved_Site_Context_scores.txt.zip
            raw_paths["mmu:${current_database_name}"]="${raw_primary_path}targetscan_Conserved_Site_Context_Scores_mmu.txt" # downloaded from http://www.targetscan.org//mmu_50//mmu_50_data_download/Conserved_Site_Context_scores.txt.zip
            raw_paths["rno:${current_database_name}"]="${raw_primary_path}targetscan_Conserved_Site_Details_mammalian.txt" # downloaded from http://www.targetscan.org//mmu_50//mmu_50_data_download/Conserved_Site_Context_scores.txt.zip
            raw_paths["gga:${current_database_name}"]="${raw_primary_path}targetscan_Conserved_Site_Details_mammalian.txt" # downloaded from http://www.targetscan.org//mmu_50//mmu_50_data_download/Conserved_Site_Context_scores.txt.zip

            local input_database_path=${raw_paths["${organism}:${current_database_name}"]}
            local output_database_path=${predictive_paths["$current_database_name"]}

            local raw_scores_path="raw_scores.tmp"
            local raw_scores_cutoff_path="raw_scores_cutoff.tmp"
            
            if [ -f $input_database_path ] ; then    
                echo "Normalizing $organism $current_database_name database"
                cat $input_database_path | 
                ${organism}_normalize_${current_database_name} | 
                filter_unique_interactions | 
                add_gene_names | 
                order_columns_add_database_name $current_database_name |
                sort -t$sep -gk4 > $raw_scores_path

                cat $raw_scores_path |
                set_size_cutoff > $raw_scores_cutoff_path
                
                normalize_scores $raw_scores_cutoff_path "invert" |
                sort -t$sep -grk4 > $output_database_path

                rm $raw_scores_path
                
                echo_summary $output_database_path
            else
                echo "Path $input_database_path doesn't exist"
            fi
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}


choose_predictive_pita_columns(){
    # remove header, choose first refseq_id, choose mirna, refseq_id and score columns
    awk -F"\t" -v OFS=$sep 'NR > 1 {split($1,a,";");$1=a[1];print $3, $1, $5}'
}

current_database_name="predictive_pita"; predictive_paths["$current_database_name"]="${normalized_organism_primary_path}${current_database_name}${database_file_extension}"
normalize_predictive_pita(){
    local current_database_name=$(echo ${FUNCNAME#*normalize_})
    case $organism in
        hsa|mmu|dme|cel)

            raw_paths["hsa:${current_database_name}"]="${raw_primary_path}PITA_targets_hg18_0_0_TOP.tab" # downloaded from http://genie.weizmann.ac.il/pubs/mir07/catalogs/PITA_sites_hg18_0_0_TOP.tab.gz
            raw_paths["mmu:${current_database_name}"]="${raw_primary_path}PITA_targets_mm9_0_0_TOP.tab" # downloaded from http://genie.weizmann.ac.il/pubs/mir07/catalogs/PITA_sites_mm9_0_0_TOP.tab.gz
            raw_paths["dme:${current_database_name}"]="${raw_primary_path}PITA_targets_dm3_0_0_TOP.tab" # downloaded from http://genie.weizmann.ac.il/pubs/mir07/catalogs/PITA_sites_mm9_0_0_TOP.tab.gz
            raw_paths["cel:${current_database_name}"]="${raw_primary_path}PITA_targets_ce6_0_0_TOP.tab" # downloaded from http://genie.weizmann.ac.il/pubs/mir07/catalogs/PITA_sites_mm9_0_0_TOP.tab.gz

            local input_database_path=${raw_paths["${organism}:${current_database_name}"]}
            local output_database_path=${predictive_paths["$current_database_name"]}

            local raw_scores_path="raw_scores.tmp"
            local raw_scores_cutoff_path="raw_scores_cutoff.tmp"
            
            if [ -f $input_database_path ] ; then    
                echo "Normalizing $organism $current_database_name database"
                cat $input_database_path | 
                choose_${current_database_name}_columns | 
                refseq_id_to_genbank_id |
                filter_unique_interactions | 
                add_gene_names | 
                order_columns_add_database_name $current_database_name | 
                sort -t$sep -gk4 > $raw_scores_path

                cat $raw_scores_path |
                set_size_cutoff > $raw_scores_cutoff_path
                
                normalize_scores $raw_scores_cutoff_path "invert" | 
                sort -t$sep -grk4 > $output_database_path

                rm $raw_scores_path
                
                echo_summary $output_database_path
            else
                echo "Path $input_database_path doesn't exist"
            fi
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}


hsa_normalize_predictive_microt(){
    # remove header, filter organism rows, choose mirna, genbank and score columns, 
    awk -F"|" -v OFS=$sep '$2 ~ /hsa/ {print $2, $3, $4}'
}
mmu_normalize_predictive_microt(){
    # remove header, filter organism rows, choose mirna, genbank and score columns, 
    awk -F"|" -v OFS=$sep '$2 ~ /mmu/ {print $2, $3, $4}'
}
dme_normalize_predictive_microt(){
    # remove header, filter organism rows, choose mirna, genbank and score columns, 
    awk -F"|" -v OFS=$sep '$2 ~ /dme/ {print $2, $3, $4}'
}

current_database_name="predictive_microt"; predictive_paths["$current_database_name"]="${normalized_organism_primary_path}${current_database_name}${database_file_extension}"
normalize_predictive_microt(){
    local current_database_name=$(echo ${FUNCNAME#*normalize_})
    case $organism in
        hsa|mmu|dme)

            raw_paths["${current_database_name}"]="${raw_primary_path}microT_v4.0.txt" # downloaded from http://diana.cslab.ece.ntua.gr/data/public/microT_v3.0.txt.gz

            local input_database_path=${raw_paths["${current_database_name}"]}
            local output_database_path=${predictive_paths["$current_database_name"]}

            local raw_scores_path="raw_scores.tmp"
            local raw_scores_cutoff_path="raw_scores_cutoff.tmp"
            
            if [ -f $input_database_path ] ; then    
                echo "Normalizing $organism $current_database_name database"
                cat $input_database_path | 
                ${organism}_normalize_${current_database_name} | 
                ensembl_gene_id_to_genbank_id | 
                filter_unique_interactions | 
                add_gene_names | 
                order_columns_add_database_name $current_database_name | 
                sort -t$sep -grk4 > $raw_scores_path

                cat $raw_scores_path |
                set_size_cutoff > $raw_scores_cutoff_path
                
                normalize_scores $raw_scores_cutoff_path |
                sort -t$sep -grk4 > $output_database_path

                rm $raw_scores_path
                
                echo_summary $output_database_path
            else
                echo "Path $input_database_path doesn't exist"
            fi
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}


choose_predictive_eimmo_columns(){
    # remove header, choose first refseq_id, choose mirna, refseq_id and score columns
    awk -F"\t" -v OFS=$sep 'NR > 1 {split($7,a,":");print a[2], a[1], $9}'
}

current_database_name="predictive_eimmo"; predictive_paths["$current_database_name"]="${normalized_organism_primary_path}${current_database_name}${database_file_extension}"
normalize_predictive_eimmo(){
    local current_database_name=$(echo ${FUNCNAME#*normalize_})
    case $organism in
        hsa|mmu|rno|dre|dme|cel)

            raw_paths["hsa:${current_database_name}"]="${raw_primary_path}eimmo_hg_targets_FullList_flat.tab" # downloaded from http://www.mirz.unibas.ch/ElMMo3/BulkDownloads/v4/hg_targets_FullList_flat.tab.gz
            raw_paths["mmu:${current_database_name}"]="${raw_primary_path}eimmo_mm_targets_FullList_flat.tab" # downloaded from http://www.mirz.unibas.ch/ElMMo3/BulkDownloads/v4/mm_targets_FullList_flat.tab.gz
            raw_paths["rno:${current_database_name}"]="${raw_primary_path}eimmo_rn_targets_FullList_flat.tab" # downloaded from http://www.mirz.unibas.ch/ElMMo3/BulkDownloads/v4/mm_targets_FullList_flat.tab.gz
            raw_paths["dre:${current_database_name}"]="${raw_primary_path}eimmo_danRer_targets_FullList_flat.tab" # downloaded from http://www.mirz.unibas.ch/ElMMo3/BulkDownloads/v4/mm_targets_FullList_flat.tab.gz
            raw_paths["dme:${current_database_name}"]="${raw_primary_path}eimmo_dm_targets_FullList_flat.tab" # downloaded from http://www.mirz.unibas.ch/ElMMo3/BulkDownloads/v4/mm_targets_FullList_flat.tab.gz
            raw_paths["cel:${current_database_name}"]="${raw_primary_path}eimmo_ce_targets_FullList_flat.tab" # downloaded from http://www.mirz.unibas.ch/ElMMo3/BulkDownloads/v4/mm_targets_FullList_flat.tab.gz


            local input_database_path=${raw_paths["${organism}:${current_database_name}"]}
            local output_database_path=${predictive_paths["$current_database_name"]}

            local raw_scores_path="raw_scores.tmp"
            local raw_scores_cutoff_path="raw_scores_cutoff.tmp"
            
            if [ -f $input_database_path ] ; then    
                echo "Normalizing $organism $current_database_name database"
                cat $input_database_path | 
                choose_${current_database_name}_columns | 
                refseq_id_to_genbank_id |
                filter_unique_interactions | 
                add_gene_names | 
                order_columns_add_database_name $current_database_name | 
                sort -t$sep -grk4 > $raw_scores_path

                cat $raw_scores_path |
                set_size_cutoff > $raw_scores_cutoff_path
                
                normalize_scores $raw_scores_cutoff_path |
                sort -t$sep -grk4 > $output_database_path

                rm $raw_scores_path
                
                echo_summary $output_database_path
            else
                echo "Path $input_database_path doesn't exist"
            fi
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}


hsa_normalize_predictive_mirtarget(){
    # remove header, filter organism, choose mirna, refseq_id and score columns
    awk -F"\t" -v OFS=$sep 'NR > 1 && $1 ~ /hsa/ {print $1, $2, $3}'
}
mmu_normalize_predictive_mirtarget(){
    # remove header, filter organism, choose mirna, refseq_id and score columns
    awk -F"\t" -v OFS=$sep 'NR > 1 && $1 ~ /mmu/ {print $1, $2, $3}'
}
rno_normalize_predictive_mirtarget(){
    # remove header, filter organism, choose mirna, refseq_id and score columns
    awk -F"\t" -v OFS=$sep 'NR > 1 && $1 ~ /rno/ {print $1, $2, $3}'
}
gga_normalize_predictive_mirtarget(){
    # remove header, filter organism, choose mirna, refseq_id and score columns
    awk -F"\t" -v OFS=$sep 'NR > 1 && $1 ~ /gga/ {print $1, $2, $3}'
}

current_database_name="predictive_mirtarget"; predictive_paths["$current_database_name"]="${normalized_organism_primary_path}${current_database_name}${database_file_extension}"
normalize_predictive_mirtarget(){
    local current_database_name=$(echo ${FUNCNAME#*normalize_})
    case $organism in
        hsa|mmu|rno|gga)
            
            raw_paths["${current_database_name}"]="${raw_primary_path}MirTarget2_v3.0_prediction_result.txt" # downloaded from http://mirdb.org/miRDB/download/MirTarget2_v3.0_prediction_result.txt.gz

            local input_database_path=${raw_paths["${current_database_name}"]}
            local output_database_path=${predictive_paths["$current_database_name"]}

            local raw_scores_path="raw_scores.tmp"
            local raw_scores_cutoff_path="raw_scores_cutoff.tmp"
            
            if [ -f $input_database_path ] ; then    
                echo "Normalizing $organism $current_database_name database"
                cat $input_database_path | 
                ${organism}_normalize_${current_database_name} | 
                refseq_id_to_genbank_id |
                filter_unique_interactions | 
                add_gene_names | 
                order_columns_add_database_name $current_database_name | 
                sort -t$sep -grk4 > $raw_scores_path

                cat $raw_scores_path |
                set_size_cutoff > $raw_scores_cutoff_path
                
                normalize_scores $raw_scores_cutoff_path |
                sort -t$sep -grk4 > $output_database_path
                
                rm $raw_scores_path
                
                echo_summary $output_database_path
            else
                echo "Path $input_database_path doesn't exist"
            fi
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}


hsa_normalize_predictive_targetspy(){
    # normalize mirna, choose refseq_id and score columns
    awk -F"\t" -v OFS=$sep '{print "hsa-"$1, $2, $10}'
}
mmu_normalize_predictive_targetspy(){
    # normalize mirna, choose refseq_id and score columns
    awk -F"\t" -v OFS=$sep '{print "mmu-"$1, $2, $10}'
}
rno_normalize_predictive_targetspy(){
    # normalize mirna, choose refseq_id and score columns
    awk -F"\t" -v OFS=$sep '{print "rno-"$1, $2, $10}'
}
gga_normalize_predictive_targetspy(){
    # normalize mirna, choose refseq_id and score columns
    awk -F"\t" -v OFS=$sep '{print "gga-"$1, $2, $10}'
}
dme_normalize_predictive_targetspy(){
    # normalize mirna, choose refseq_id and score columns
    awk -F"\t" -v OFS=$sep '{print "dme-"$1, $2, $10}'
}

current_database_name="predictive_targetspy"; predictive_paths["$current_database_name"]="${normalized_organism_primary_path}${current_database_name}${database_file_extension}"
normalize_predictive_targetspy(){
    local current_database_name=$(echo ${FUNCNAME#*normalize_})
    case $organism in
        hsa|mmu|rno|gga)

            raw_paths["hsa:${current_database_name}"]="${raw_primary_path}targetspy_hsa_refseq_seed_sens" # downloaded from http://www.targetspy.org/data/hsa_refseq_seed_sens.gz
            raw_paths["mmu:${current_database_name}"]="${raw_primary_path}targetspy_mmu_refseq_seed_sens" # downloaded from http://www.targetspy.org/data/mmu_refseq_seed_sens.gz
            raw_paths["rno:${current_database_name}"]="${raw_primary_path}targetspy_rno_refseq_seed_sens" # downloaded from http://www.targetspy.org/data/mmu_refseq_seed_sens.gz
            raw_paths["gga:${current_database_name}"]="${raw_primary_path}targetspy_gga_refseq_seed_sens" # downloaded from http://www.targetspy.org/data/mmu_refseq_seed_sens.gz
            raw_paths["dme:${current_database_name}"]="${raw_primary_path}targetspy_dme_refseq_seed_sens" # downloaded from http://www.targetspy.org/data/mmu_refseq_seed_sens.gz

            local input_database_path=${raw_paths["${organism}:${current_database_name}"]}
            local output_database_path=${predictive_paths["$current_database_name"]}

            local raw_scores_path="raw_scores.tmp"
            local raw_scores_cutoff_path="raw_scores_cutoff.tmp"

            if [ -f $input_database_path ] ; then    
                echo "Normalizing $organism $current_database_name database"
                cat $input_database_path | 
                ${organism}_normalize_${current_database_name} | 
                refseq_id_to_genbank_id |
                filter_unique_interactions | 
                add_gene_names | 
                order_columns_add_database_name $current_database_name | 
                sort -t$sep -gk4 > $raw_scores_path

                cat $raw_scores_path | 
                set_size_cutoff > $raw_scores_cutoff_path
                
                normalize_scores $raw_scores_cutoff_path "invert" |
                sort -t$sep -grk4 > $output_database_path

                rm $raw_scores_path
                
                echo_summary $output_database_path
            else
                echo "Path $input_database_path doesn't exist"
            fi
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}


# hsa_normalize_predictive_reptar(){
#     # filter human, normalize mirna, normalize refseq, normalize score, 
#     awk -F"\t" -v OFS=$sep '$2 ~ /hsa/ {split($1,a,":::"); $1=a[2]; sub(/star/,"*",$2); split($6,a,":"); $6=a[2];print $2, $1, $6}'
# }
# mmu_normalize_predictive_reptar(){
#     # filter mouse, normalize mirna, normalize refseq, normalize score, 
#     awk -F"\t" -v OFS=$sep '$2 ~ /mmu/ {split($1,a,":::"); $1=a[2]; sub(/star/,"*",$2); split($6,a,":"); $6=a[2];print $2, $1, $6}' 
# }

# current_database_name="predictive_reptar"; predictive_paths["$current_database_name"]="${normalized_organism_primary_path}${current_database_name}${database_file_extension}"
# normalize_predictive_reptar(){
#             local current_database_name=$(echo ${FUNCNAME#*normalize_})
#     case $organism in
#         hsa|mmu)

#             raw_paths["hsa:${current_database_name}"]="${raw_primary_path}reptar_human-predictions.txt" # downloaded from http://bioinformatics.ekmd.huji.ac.il/reptar/files/human-predictions.txt.gz
#             raw_paths["mmu:${current_database_name}"]="${raw_primary_path}reptar_mouse-predictions.txt" # downloaded from http://bioinformatics.ekmd.huji.ac.il/reptar/files/mouse-predictions.txt.gz            

#             local input_database_path=${raw_paths["${organism}:${current_database_name}"]}
#             local output_database_path=${predictive_paths["$current_database_name"]}

#             local raw_scores_path="raw_scores.tmp"
#             local raw_scores_cutoff_path="raw_scores_cutoff.tmp"
            
#             if [ -f $input_database_path ] ; then    
#                 echo "Normalizing $organism $current_database_name database"
#                 cat $input_database_path | 
#                 ${organism}_normalize_${current_database_name} | 
#                 refseq_id_to_genbank_id |
#                 filter_unique_interactions | 
#                 add_gene_names | 
#                 order_columns_add_database_name $current_database_name | 
#                 sort -t$sep -grk4 > $raw_scores_path

#                 cat $raw_scores_path |
#                 set_size_cutoff > $raw_scores_cutoff_path
                
#                 normalize_scores $raw_scores_cutoff_path |
#                 sort -t$sep -grk4 > $output_database_path
                
#                 rm $raw_scores_path
                
#                 echo_summary $output_database_path
#             else
#                 echo "Path $input_database_path doesn't exist"
#             fi
#             ;;
#         *)
#             error_wrong_organism
#             ;;
#     esac
# }


choose_predictive_mirwalk_columns(){
    # remove header, choose mirna, refseq and score columns, 
    awk -F"\t" -v OFS=$sep 'NR > 1 {print $1, $3, $4}'
}

current_database_name="predictive_mirwalk"; predictive_paths["$current_database_name"]="${normalized_organism_primary_path}${current_database_name}${database_file_extension}"
normalize_predictive_mirwalk(){
    local current_database_name=$(echo ${FUNCNAME#*normalize_})
    case $organism in
        hsa)
            
            raw_paths["hsa:${current_database_name}"]="${raw_primary_path}predictive_mirwalk_hsa.txt" # downloaded from http://www.ma.uni-heidelberg.de/apps/zmf/mirwalk/micrornapredictedtarget.html

            local input_database_path=${raw_paths["${organism}:${current_database_name}"]}
            local output_database_path=${predictive_paths["$current_database_name"]}

            local raw_scores_path="raw_scores.tmp"
            local raw_scores_cutoff_path="raw_scores_cutoff.tmp"
            local raw_scores_cutoff_path="raw_scores_cutoff.tmp"
            
            if [ -f $input_database_path ] ; then    
                echo "Normalizing $organism $current_database_name database"
                cat $input_database_path | 
                choose_${current_database_name}_columns | 
                refseq_id_to_genbank_id |  
                filter_unique_interactions | 
                add_gene_names | 
                order_columns_add_database_name $current_database_name | 
                sort -t$sep -gk4 > $raw_scores_path
                
                cat $raw_scores_path | 
                set_size_cutoff > $raw_scores_cutoff_path
                
                normalize_scores $raw_scores_cutoff_path "invert" |
                sort -t$sep -grk4 > $output_database_path

                rm $raw_scores_path
                rm $raw_scores_cutoff_path
                
                echo_summary $output_database_path
            else
                echo "Path $input_database_path doesn't exist"
            fi
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}

################################################################################


########### CALCULATIONS #######################################################

# Concatenates all available experimental databases for a given organism, removing duplicates 
current_database_name="combined_experimental";combined_experimental_path="${normalized_organism_calculations_path}${current_database_name}${database_file_extension}"
calculate_combined_experimental(){
    local current_database_name=$(echo ${FUNCNAME#*calculate_})
    case $organism in
        hsa|mmu|rno|gga|dre|dme|cel)
            local input_database_paths=${experimental_paths[*]}
            local output_database_path="$combined_experimental_path"
            
            echo "Calculating $organism $current_database_name database"
            cut -d$sep -f1,2 $input_database_paths | 
            sort -u > $output_database_path 
            
            echo_summary $output_database_path
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}

# Calculates the number of validated interactions that appear in each predictive database for a given organism. It doesn't start writing precisions until at least minimum_TP validated interactions have been found (this solves the initial fluctuation)
# For convenience, it allows a third argument: the predictive database name (otherwise it uses all the predictive databases)
# If a third argument is used, any fourth argument will also generate the corresponding plot (requires weighed_randomized_predictive)
minimum_TP=20
precision_prefix="precision_"
for predictive_database_name in ${!predictive_paths[@]} ; do
    precision_paths["$predictive_database_name"]="${normalized_organism_calculations_path}${precision_prefix}${predictive_database_name}${database_file_extension}"
done
calculate_precision(){
    local current_database_name=$(echo ${FUNCNAME#*calculate_})
    case $organism in
        hsa|mmu|rno|gga|dre|dme|cel)
            local input_database_paths
            if [[ -n "$specific_database_name" ]] ; then
                echo "For database $specific_database_name"
                input_database_paths=${predictive_paths["$specific_database_name"]}
            else
                echo "For all primary predictive databases"
                input_database_paths=${predictive_paths[*]}
            fi
        
            for input_database_path in $input_database_paths ; do
                if [[ -f $input_database_path ]] ; then
                    local input_database_name="$(get_name $input_database_path)"
                    local output_database_path=${precision_paths["$input_database_name"]}

                    echo "Calculating $organism $current_database_name $input_database_name database"
                    cat $input_database_path |
                    add_precision $combined_experimental_path > $output_database_path

                    echo_summary $output_database_path
                else
                    echo "Path $input_database_path doesn't exist"    
                fi
            done
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}


# It combines all the predictive databases using their individual scores and precisions 
# It also calculates the precision of the combined database to evaluate its performance
current_database_name="combined_precision";combined_precision_path="${normalized_organism_calculations_path}${current_database_name}${database_file_extension}"
calculate_combined_precision(){
    local current_database_name=$(echo ${FUNCNAME#*calculate_})
    case $organism in
        hsa|mmu|rno|gga|dre|dme|cel)
            local input_database_paths=${precision_paths[*]}
            local raw_scores_path="${normalized_organism_calculations_path}raw_combined_precision_scores.tmp"
            local output_database_path=$combined_precision_path

            echo "Calculating $organism $current_database_name database"
            cat $input_database_paths | 
            awk -F$sep -v OFS=$sep -f $combine_predictive_databases_script |
            sort -t$sep -grk4 > $raw_scores_path
            
            normalize_scores $raw_scores_path |
            sort -t$sep -grk4 |
            add_precision $combined_experimental_path |
            sort -t$sep -grk4 > $output_database_path

            rm $raw_scores_path
            
            echo_summary $output_database_path
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}


roc_curve_prefix="roc_curve_"
for predictive_database_name in ${!predictive_paths[@]} ; do
    roc_curve_paths["$predictive_database_name"]="${normalized_organism_calculations_path}${roc_curve_prefix}${predictive_database_name}${database_file_extension}"
done
calculate_roc_curve(){
    local current_database_name=$(echo ${FUNCNAME#*calculate_})
    case $organism in
        hsa|mmu|rno|gga|dre|dme|cel)
            local score_increment=0.01
            local input_database_paths
            if [[ -n "$specific_database_name" ]] ; then
                echo "For database $specific_database_name"
                input_database_paths=${predictive_paths["$specific_database_name"]}
            else
                echo "For all primary predictive databases"
                input_database_paths=${predictive_paths[*]}
            fi

            for input_database_path in $input_database_paths ; do
                if [[ -f $input_database_path ]] ; then
                    local input_database_name="$(get_name $input_database_path)"
                    local output_database_path=${roc_curve_paths["$input_database_name"]}
                    
                    echo "Calculating $organism $current_database_name $input_database_name database"
                    awk -F$sep -v increment=$score_increment -v reference=$combined_experimental_path -v input=$input_database_path -f $roc_curve_script > $output_database_path
   
                    echo_summary $output_database_path 
                else
                    echo "Path $input_database_path doesn't exist"    
                fi
            done
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}

# Since there is no combined_predictive database, we use combined_precision as input
current_database_name="combined_roc_curve";combined_roc_curve_path="${normalized_organism_calculations_path}${current_database_name}${database_file_extension}"
calculate_combined_roc_curve(){
    local current_database_name=$(echo ${FUNCNAME#*calculate_})
    case $organism in
        hsa|mmu|rno|gga|dre|dme|cel)
            local input_database_path=$combined_precision_path
            local output_database_path="$combined_roc_curve_path"

            local score_increment=0.01
            if [[ -f $input_database_path ]] ; then
                echo "Calculating $organism $current_database_name database"
                awk -F$sep -v increment=$score_increment -v reference=$combined_experimental_path -v input=$input_database_path -f $roc_curve_script > $output_database_path

                echo_summary $output_database_path 
            else
                echo "Path $input_database_path doesn't exist"
            fi
            ;;
        *)
            error_wrong_organism
            ;;
    esac
}

################################################################################


########### PRIMARY FUNCTIONS ##################################################

remove_blank_genbank_ids(){
    awk -F$sep -v OFS=$sep '$4 != "" {print $1, $4, $3}'
}

refseq_id_to_genbank_id(){
    awk -F$sep -v OFS=$sep -v translator=${auxiliary_paths["refseq_id_to_genbank_id"]} -f $translator_script |
    remove_blank_genbank_ids
}

ensembl_gene_id_to_genbank_id(){
    awk -F$sep -v OFS=$sep -v translator=${auxiliary_paths["ensembl_gene_id_to_genbank_id"]} -f $translator_script |
    remove_blank_genbank_ids
}

ensembl_transcript_id_to_genbank_id(){
   awk -F$sep -v OFS=$sep -v translator=${auxiliary_paths["ensembl_transcript_id_to_genbank_id"]} -f $translator_script |
   remove_blank_genbank_ids
}

filter_unique_interactions(){
    awk -F$sep -v OFS=$sep -f $unique_interactions_script
}

add_gene_names(){
    awk -F$sep -v OFS=$sep -v translator=${auxiliary_paths["genbank_id_to_gene_name"]} -f $translator_script
}

order_columns_add_database_name(){
    awk -F$sep -v OFS=$sep -v db_name=$1 '{print $1, $2, $4, $3, db_name}'
}

# We're adding this score to be able to sort numerically in excel, and chose 10 as an arbitrary value that would make experimental interactions score much higher than predicted interactions
experimental_score=10
add_experimental_score(){
    awk -F$sep -v OFS=$sep -v experimental_score=$experimental_score '{print $1, $2, $3, experimental_score$4, $5}'
}

set_size_cutoff(){
    awk -F$sep -v OFS=$sep -v cutoffs_file=$cutoffs_file -v current_database_name=$current_database_name -f $cutoff_script
}    

normalize_scores(){
    local invert=0
    if [[ $# > 1 ]]; then
        invert=1
    fi
    awk -F$sep -v OFS=$sep -v invert=$invert -f $score_normalization_script $1 $1
}

# Requires the reference file as first parameter
add_precision(){
    awk -F$sep -v OFS=$sep -v minimum_TP=$minimum_TP -v reference=$1 -f $precision_script
}
################################################################################


########### PLOTS ##############################################################

plot_score_distribution() {
    if [[ -z $specific_database_name ]] ; then
        for predictive_path in ${predictive_paths[@]} ; do
            Rscript "${lib_path}plot_score_distribution.R" "$plots_organism_path" "$predictive_path"
        done
    else
        Rscript "${lib_path}plot_score_distribution.R" "$plots_organism_path" "${predictive_paths[${specific_database_name}]}"
    fi
}


plot_precision() {
    if [[ -z $specific_database_name ]] ; then
        for precision_path in ${precision_paths[@]} ; do
            Rscript "${lib_path}plot_precision.R" "$plots_organism_path" "$precision_path"
        done
    else
        Rscript "${lib_path}plot_precision.R" "$plots_organism_path" "${precision_paths[${specific_database_name}]}"
    fi
}

plot_combined_precision() {
    Rscript "${lib_path}plot_precision.R" "$plots_organism_path" "$combined_precision_path"
}

plot_precision_stacked() {
Rscript "${lib_path}plot_precision_stacked.R" "$plots_organism_path" "normalized/hsa/calculations/precision_predictive_eimmo.csv" "normalized/hsa/calculations/precision_predictive_mirtarget.csv" "normalized/hsa/calculations/precision_predictive_targetspy.csv" "normalized/hsa/calculations/precision_predictive_microrna_org_conserved_high_score.csv" "normalized/hsa/calculations/precision_predictive_mirwalk.csv" "normalized/hsa/calculations/precision_predictive_pita.csv" "normalized/hsa/calculations/precision_predictive_microcosm.csv" "normalized/hsa/calculations/precision_predictive_targetscan_conserved.csv" "normalized/hsa/calculations/precision_predictive_microt.csv" "$combined_precision_path"

}

plot_roc_curve() {
    if [[ -z $specific_database_name ]] ; then
        for roc_curve_path in ${roc_curve_paths[@]} ; do
            Rscript "${lib_path}plot_roc_curve.R" "$plots_organism_path" "$roc_curve_path"
        done
    else
        Rscript "${lib_path}plot_roc_curve.R" "$plots_organism_path" "${roc_curve_paths[${specific_database_name}]}"
    fi
}

plot_combined_roc_curve() {
    Rscript "${lib_path}plot_roc_curve.R" "$plots_organism_path" "$combined_roc_curve_path"
}
################################################################################


########### AUXILIARY FUNCTIONS ################################################

get_name(){
    local name=$(basename $1)
    echo ${name%.*}
}

echo_summary(){
    echo "$(wc -l $1)"
    echo "$(head -n 5 $1)"
    echo "..."
    echo "$(tail -n 5 $1)"
    echo
}

error_wrong_organism(){
    echo
    echo "Organism '$organism' isn't valid for '$current_database_name' database"
    echo
}

do_all_auxiliary(){
    echo "########################## Generating all_auxiliary ##########################"
    normalize_genbank_id_to_gene_name
    normalize_refseq_id_to_genbank_id
    normalize_ensembl_transcript_id_to_genbank_id
    normalize_ensembl_gene_id_to_genbank_id
    normalize_mirnas
    echo "########################## Finished all_auxiliary ##########################"
    echo
}

do_all_experimental(){
    echo "########################## Generating all_experimental ##########################"
    normalize_experimental_mirtarbase
    normalize_experimental_tarbase
    normalize_experimental_mirecords
    normalize_experimental_mirwalk
    echo "########################## Finished all_experimental ##########################"
    echo
}

do_all_predictive(){
    echo "########################## Generating all_predictive ##########################"
    normalize_predictive_mirtarget
    normalize_predictive_eimmo
    normalize_predictive_microrna_org_conserved_high_score
    normalize_predictive_microcosm
    normalize_predictive_targetscan_conserved
    normalize_predictive_pita
    normalize_predictive_microt
    normalize_predictive_targetspy
    normalize_predictive_mirwalk
    # normalize_predictive_reptar
    echo "########################## Finished all_predictive ##########################"
    echo
}

do_all_calculations(){
    echo "########################## Generating all_calculations ##########################"
    calculate_combined_experimental
    calculate_precision
    calculate_combined_precision
    calculate_roc_curve
    calculate_combined_roc_curve
    echo "########################## Finished all_calculations ##########################"
    echo
}

do_all_plots(){
    echo "########################## Generating all_plots ##########################"
    plot_score_distribution
    plot_precision
    plot_combined_precision
    plot_roc_curve
    plot_combined_roc_curve
    echo "########################## Finished all_plots ##########################"
    echo
}

do_all(){
    do_all_auxiliary
    do_all_experimental
    do_all_predictive
    do_all_calculations
    do_all_plots
}

case $query_database in
    # AUXILIARY DATABASES
    genbank_id_to_gene_name)
        normalize_genbank_id_to_gene_name
        ;;
    refseq_id_to_genbank_id)
        normalize_refseq_id_to_genbank_id
        ;;
    ensembl_transcript_id_to_genbank_id)
        normalize_ensembl_transcript_id_to_genbank_id
        ;;
    ensembl_gene_id_to_genbank_id)
        normalize_ensembl_gene_id_to_genbank_id
        ;;
    mirnas)
        normalize_mirnas
        ;;
    all_auxiliary)
        calculate_all_auxiliary
        ;;
    # PRIMARY DATABASES (EXPERIMENTAL)
    experimental_mirtarbase)
        normalize_experimental_mirtarbase
        ;;
    experimental_mirwalk)
        normalize_experimental_mirwalk
        ;;
    experimental_mirecords)
        normalize_experimental_mirecords
        ;;
    experimental_tarbase)
        normalize_experimental_tarbase
        ;;    
    all_experimental)
        do_all_experimental
        ;;        
    # PRIMARY DATABASES (PREDICTIVE)
    predictive_targetspy)
        normalize_predictive_targetspy
        ;;
    predictive_mirtarget)
        normalize_predictive_mirtarget
        ;;
    predictive_eimmo)
        normalize_predictive_eimmo
        ;;
    predictive_microt)
        normalize_predictive_microt
        ;;
    predictive_pita)
        normalize_predictive_pita
        ;;
    predictive_mirwalk)
        normalize_predictive_mirwalk
        ;;
    predictive_microrna_org_conserved_high_score)
        normalize_predictive_microrna_org_conserved_high_score
        ;;
    predictive_microcosm)
        normalize_predictive_microcosm
        ;;
    predictive_targetscan_conserved)
        normalize_predictive_targetscan_conserved
        ;;
    # predictive_reptar)
    #     normalize_predictive_reptar
    #     ;;
    all_predictive)
        do_all_predictive
        ;;
    # CALCULATIONS    
    combined_experimental)       
        calculate_combined_experimental
        ;;
    precision)
        calculate_precision
        ;;
    combined_precision)       
        calculate_combined_precision
        ;;
    roc_curve)
        calculate_roc_curve
        ;;
    combined_roc_curve)
        calculate_combined_roc_curve
        ;;
    all_calculations)
        do_all_calculations
        ;;
    plot_score_distribution)
        plot_score_distribution
        ;;
    plot_precision)
        plot_precision
        ;;
    plot_combined_precision)
        plot_combined_precision
        ;;
    plot_roc_curve)
        plot_roc_curve
        ;;
    plot_combined_roc_curve)
        plot_combined_roc_curve
        ;;
    plot_precision_stacked)
        plot_precision_stacked
        ;;
    all_plots)
        do_all_plots
        ;;
    all)
        do_all
        ;;        
    *)    
        echo "Database $2 doesn't exist"
        echo "Available databases:"
        echo $all_databases
        exit
        ;;
esac

################################################################################