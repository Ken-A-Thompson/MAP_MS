#!/bin/bash
clear
# Description: This script will take a FASTA file of sequences and, if user selects "blast", match them against a BOLDistilled reference library and sort the results into BIN matches and no-BIN-matches.
#              If user selects "sintax", the script assigns probable taxonomy based on comparison to a BOLDistilled reference library.

# Author: Sean Prosser (May 2025)

# Instructions for use: 'cd' into the directory containing your FASTA file. Call the script and provide parameters like this:
#                       bash /path/to/script/BIN_Match.sh -i your_fasta.fa -o path/to/output/directory -r path/to/BOLDistilled_COI_XXXXX -m blast

# initialize variables and parameters
cores=$( [ "$(uname)" = "Darwin" ] && sysctl -n hw.ncpu || nproc )
amp_size=658
min_overlap=$(echo "$amp_size * 0.76" | bc)  #76% of amp.size
sintax_conf=0.6
input_fasta=""
reflib=""
output_directory=""
method=""

# Parse flags
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -r|--reflib)
            reflib="$2"
            shift 2
            ;;
        -i|--input)
            input_fasta="$2"
            shift 2
            ;;
        -o|--output)
            output_directory="$2"
            shift 2
            ;;
        -m|--method)
        method="$2"
        shift 2
        ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

################# VSEARCH BLAST ##############################
if [ "$method" = "blast" ]; then
    db_file=$(echo $reflib/PUBLIC/VSEARCH/*SEQUENCES_vsearch)
    tax_file=$(echo $reflib/PUBLIC/*TAXONOMY.tsv)

    # convert FASTA to single-line
    awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $input_fasta > single_line.fasta
    rm $input_fasta
    mv single_line.fasta $input_fasta

    # identify sequences using VSEARCH
    vsearch --usearch_global $input_fasta \
        --db  $db_file \
        --blast6out temp_vsearch_output.txt \
        --id 0.75 \
        --maxhits 3 \
        --maxaccepts 3 \
        --threads $((cores-1))

    # process the VSEARCH output file
    awk -v min_overlap="$min_overlap" '
    BEGIN { OFS="\t" }
    $4 >= min_overlap { print }
    ' temp_vsearch_output.txt | sort -k1,1 -k3nr -k4nr | awk '!seen[$1]++' > filtered_vsearch_output.txt
    awk -F'\t' 'BEGIN { OFS="\t" } $2 ~ /\|/ { split($2, hit_parts, "|"); print $1, hit_parts[1], hit_parts[2], $3, $4 }' filtered_vsearch_output.txt > parsed_vsearch_output.txt

    #import taxonomy of BIN hits
    awk -F'\t' 'BEGIN { OFS="\t" }
    NR==FNR { tax[$1] = $0; next } 
    { print $0, (tax[$3] ? tax[$3] : "NA") }' $tax_file parsed_vsearch_output.txt > final_output.txt

    awk -F'\t' 'BEGIN { OFS = FS } { $6 = ""; sub(/\t\t/, "\t"); print }' final_output.txt > final_output_no_col6.txt


    # add "bin_match" column
    awk -F'\t' 'BEGIN { OFS="\t" }
    {
        # Check if %id (column 4) is greater than or equal to 97.7%
        if ($4 >= 97.7) {
            bin_match = "BIN_MATCH"
        } else {
            bin_match = "NO_MATCH"
        }
        
        # Ensure 13 columns exist, add "bin_match" as column 14
        for (i = NF + 1; i <= 13; i++) {
            $i = ""  # Add empty columns if necessary
        }
        $14 = bin_match  # Set the 14th column to "bin_match"
        print $0
    }' final_output_no_col6.txt > final_output_with_bin_match.txt

    # add header and replace NA with blanks
    (echo -e "Query\tHit (PID)\tHit (BIN)\t%ID\tOverlap (bp)\tKingdom\tPhylum\tClass\tOrder\tFamily\tSubfamily\tGenus\tSpecies\tBIN Match?"; cat final_output_with_bin_match.txt) | sed 's/\tNA\t/\t\t/g' > final_output_with_bin_match_header.txt

    # reorder columns
    awk -F'\t' 'BEGIN { OFS="\t" }
    {
        print $1, $14, $3, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13
    }' final_output_with_bin_match_header.txt > "${input_fasta%.fasta}"_BIN_MATCH_RESULTS.tsv

    mkdir -m 777 "$output_directory"
    mv "${input_fasta%.fasta}"_BIN_MATCH_RESULTS.tsv "$output_directory"

    # tidy up directory
    rm temp_vsearch_output.txt filtered_vsearch_output.txt parsed_vsearch_output.txt final_output.txt final_output_no_col6.txt final_output_with_bin_match.txt final_output_with_bin_match_header.txt

elif [ "$method" = "sintax" ]; then
    db_file=$(echo $reflib/PUBLIC/SINTAX/*SEQUENCES_sintax.fasta)

    # convert FASTA to single-line
    awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $input_fasta > single_line.fasta
    rm $input_fasta
    mv single_line.fasta $input_fasta

    # identify sequences using SINTAX
    vsearch --sintax $input_fasta \
        -db $db_file \
        -tabbedout temp_sintax_output.txt \
        -strand plus \
        -sintax_cutoff 0.6 \
        -threads $((cores-1))

    # process the SINTAX output file
    awk -F'\t' '
    BEGIN {
        # Define the order of ranks
        rank_code[1] = "k"; rank_name["k"] = "Kingdom"
        rank_code[2] = "p"; rank_name["p"] = "Phylum"
        rank_code[3] = "c"; rank_name["c"] = "Class"
        rank_code[4] = "o"; rank_name["o"] = "Order"
        rank_code[5] = "f"; rank_name["f"] = "Family"
        rank_code[6] = "g"; rank_name["g"] = "Genus"
        rank_code[7] = "s"; rank_name["s"] = "Species"

        # Print header line
        printf "Sample"
        for (i = 1; i <= 7; i++) {
            printf "\t%s", rank_name[rank_code[i]]
        }
        printf "\n"
    }

    {
        sample = $1
        split("", tax)  # clear array

        # Initialize all to empty string
        for (i = 1; i <= 7; i++) {
            tax[rank_code[i]] = ""
        }

        # Only process the 4th column
        n = split($4, fields, ",")
        for (i = 1; i <= n; i++) {
            split(fields[i], parts, ":")
            code = parts[1]
            value = parts[2]
            tax[code] = value
        }

        # Print sample + taxonomic ranks
        printf "%s", sample
        for (i = 1; i <= 7; i++) {
            printf "\t%s", tax[rank_code[i]]
        }
        printf "\n"
    }
    ' temp_sintax_output.txt > "${input_fasta%.fasta}"_TAX_ASSINGMENT_RESULTS.tsv


    mkdir -m 777 "$output_directory"
    mv "${input_fasta%.fasta}"_TAX_ASSINGMENT_RESULTS.tsv "$output_directory"

    # tidy up directory
    rm temp_sintax_output.txt
else
    echo "Error: -m can only be 'blast' or 'sintax'" >&2
    exit 1
fi