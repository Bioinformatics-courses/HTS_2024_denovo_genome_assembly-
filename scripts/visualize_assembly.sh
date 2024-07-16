#!/bin/bash

#Visualize assembled genome graph with Bandage

#define directories and file paths
WORKDIR="$HOME/denovo"
ASSEMBLY_OUT="$WORKDIR/assembly_output"
GRAPH_FILE="$ASSEMBLY_OUT/assembly.gfa"
BANDAGE_OUT="$WORKDIR/bandage_output"
IMAGE_OUT="$BANDAGE_OUT/assembly_graph.png"

#creating output directory
mkdir -p $BANDAGE_OUT


# Check if the assembly graph file exists
if [ ! -f $GRAPH_FILE ]; then
    echo "Assembly graph file not found: $GRAPH_FILE"
    exit 1
fi

echo -e "\n=============================="
echo "Visualizing assembly graph with Bandage"
echo -e "==============================\n"
Bandage image $GRAPH_FILE $IMAGE_OUT --lengths


echo "Visualization completed. Image stored in $IMAGE_OUT"

