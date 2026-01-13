#!/usr/bin/env bash

while getopts i: opt
do
    case $opt in
        i)  ID="$OPTARG" ;;
        ?)  printf "Usage: %s: meth_map.sh [-i] args\n" $0 exit 2;;
    esac
done

TSV_FILE="data/methylations/"$ID"/map_info.tsv"
FASTQ_FILES="data/methylations/"$ID"/samples/*fastq.gz"
REF="dbs/genomes"


for f in ${FASTQ_FILES[@]}; do
    P="${f%%_*}"
    BASE=$(basename "$P")
done

OUT_PATH=$"data/methylations/"$ID"/mapped"
mkdir -p $OUT_PATH

while IFS=$'\t' read -r -a MAP_ARRAY; do
    SAMPLE=${MAP_ARRAY[0]}
    echo -e "\nProcessing $SAMPLE"

    if [[ -n ${BASE[${MAP_ARRAY[0]}]} ]]; then
        FW="data/methylations/$ID/samples/${SAMPLE}_1.fastq.gz"
        RV="data/methylations/$ID/samples/${SAMPLE}_2.fastq.gz"
        SAMPLE_NAME=${MAP_ARRAY[1]}

        if [[ ! -f "$FW" || ! -f "$RV" ]]; then
            echo "Error: FASTQ pair missing for $SAMPLE"
            continue
        fi

        echo -e "\n\nMapping $SAMPLE_NAME"
        bismark -p 4 --genome $REF -1 $FW -2 $RV -o $OUT_PATH -B $SAMPLE_NAME
        else
        echo "Error: $SAMPLE not found for mapping"
    fi
done < "$TSV_FILE"