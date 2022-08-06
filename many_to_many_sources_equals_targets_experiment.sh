set -e # abort complete script on a single error

DATA_PATH=$1
GRAPH_PATH=$2
OUTPUT_PATH=$3

NUM_ITERATIONS=$4

SOURCES_TARGETS_FOLDER="sources_equals_targets"

# run all one-to-many experiments
for source_target_size in "many_to_many_2^10" "many_to_many_2^12" "many_to_many_2^14"
do
    # create output directories
    mkdir -p "${OUTPUT_PATH}/${source_target_size}/${SOURCES_TARGETS_FOLDER}"

    for algorithm in "rphast" "buckets" "simultaneous_buckets" "rphast_multiple_trees_fwd_rank" "rphast_multiple_trees_fwd_dijkstra" "rphast_multiple_trees_fwd_buckets" "simd_rphast" "lazy_batched_advanced_buckets" "lazy_simd" "forward_backward_sbi" "hub_labels_partially_pruned" "hub_labels_top_down"
    do
        cargo run --release --bin experiment_command experiment --algorithm=$algorithm --data-path="${DATA_PATH}/${SOURCES_TARGETS_FOLDER}" --graph-path=$GRAPH_PATH --source-target-size=$source_target_size --num-iterations=$NUM_ITERATIONS > "${OUTPUT_PATH}/${source_target_size}/${SOURCES_TARGETS_FOLDER}/${algorithm}.csv"
    done
done

read -rsp $'Press enter to exit...\n'