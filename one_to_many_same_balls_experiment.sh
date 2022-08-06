set -e # abort complete script on a single error

DATA_PATH=$1
GRAPH_PATH=$2
OUTPUT_PATH=$3

NUM_ITERATIONS=$4

SOURCES_TARGETS_FOLDER="same_balls"

# run all one-to-many experiments
for source_target_size in "one_to_many_2^10" "one_to_many_2^12" "one_to_many_2^14"
do
    # create output directory
    mkdir -p "${OUTPUT_PATH}/${source_target_size}/${SOURCES_TARGETS_FOLDER}"

    for algorithm in "dijkstra" "rphast" "buckets" "simultaneous_buckets" "simd_rphast" "hub_labels_partially_pruned" "hub_labels_top_down" "forward_backward_sbi" "lazy_one_to_many"
    do
        cargo run --release --bin experiment_command experiment --algorithm=$algorithm --data-path="${DATA_PATH}/${SOURCES_TARGETS_FOLDER}" --graph-path=$GRAPH_PATH --source-target-size=$source_target_size --num-iterations=$NUM_ITERATIONS > "${OUTPUT_PATH}/${source_target_size}/${SOURCES_TARGETS_FOLDER}/${algorithm}.csv"
    done
done

# run all many-to-one experiments
for source_target_size in "many_to_one_2^10" "many_to_one_2^12" "many_to_one_2^14"
do
    # create output directory
    mkdir -p "${OUTPUT_PATH}/${source_target_size}/${SOURCES_TARGETS_FOLDER}"

    for algorithm in "lazy" 
    do
        cargo run --release --bin experiment_command experiment --algorithm=$algorithm --data-path="${DATA_PATH}/${SOURCES_TARGETS_FOLDER}" --graph-path=$GRAPH_PATH --source-target-size=$source_target_size --num-iterations=$NUM_ITERATIONS > "${OUTPUT_PATH}/${source_target_size}/${SOURCES_TARGETS_FOLDER}/${algorithm}.csv"
    done
done

read -rsp $'Press enter to exit...\n'