#!/bin/bash

# script to create the relevant data for all experiments

GRAPH_PATH=$1 # first argument defines the path to the graph folder (e.g. /Users/user/files/europe)
OUTPUT_PATH=$2 # second argument defines the path to the output folder
NUM_ITERATIONS=$3 # defines the nunber of files per configuration

# create output directories
for folder in "different_balls" "same_balls" "sources_equals_targets"
do 
    mkdir -p "${OUTPUT_PATH}/${folder}"
done

for folder in "one_to_many_2^10" "one_to_many_2^12" "one_to_many_2^14" "many_to_one_2^10" "many_to_one_2^12" "many_to_one_2^14" 
do
    mkdir -p "${OUTPUT_PATH}/different_balls/${folder}/sources" "${OUTPUT_PATH}/different_balls/${folder}/targets" "${OUTPUT_PATH}/same_balls/${folder}/sources" "${OUTPUT_PATH}/same_balls/${folder}/targets" "${OUTPUT_PATH}/sources_equals_targets/${folder}/sources" "${OUTPUT_PATH}/sources_equals_targets/${folder}/targets"
done

for folder in "different_balls" "sources_equals_targets" "same_balls"
do
    for subfolder in "many_to_many_2^10" "many_to_many_2^12" "many_to_many_2^14" "many_to_many_2^10_to_2^14" "many_to_many_2^14_to_2^10"
    do
        mkdir -p "${OUTPUT_PATH}/${folder}/${subfolder}/sources" "${OUTPUT_PATH}/${folder}/${subfolder}/targets" "${OUTPUT_PATH}/${folder}/${subfolder}/sources" "${OUTPUT_PATH}/${folder}/${subfolder}/targets"
    done
done

# create sources data
cargo run --release --bin experiment_command create-data --graph-path=$GRAPH_PATH --output-path="${OUTPUT_PATH}" --num-iterations=$NUM_ITERATIONS --configuration="many_to_many_sources_equals_targets"
cargo run --release --bin experiment_command create-data --graph-path=$GRAPH_PATH --output-path="${OUTPUT_PATH}" --num-iterations=$NUM_ITERATIONS --configuration="many_to_many_same_ball"
cargo run --release --bin experiment_command create-data --graph-path=$GRAPH_PATH --output-path="${OUTPUT_PATH}" --num-iterations=$NUM_ITERATIONS --configuration="many_to_many_different_balls"

# create target data
cargo run --release --bin experiment_command create-data --graph-path=$GRAPH_PATH --output-path="${OUTPUT_PATH}" --num-iterations=$NUM_ITERATIONS --configuration="one_to_many_same_balls"
cargo run --release --bin experiment_command create-data --graph-path=$GRAPH_PATH --output-path="${OUTPUT_PATH}" --num-iterations=$NUM_ITERATIONS --configuration="one_to_many_different_balls"

read -rsp $'Press enter to exit...\n'