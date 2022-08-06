set -e # abort complete script on a single error

DATA_PAH=$1
GRAPH_PATH=$2
OUTPUT_PATH=$3
NUM_ITERATIONS=$4

SOURCES_TARGETS_FOLDER="different_balls"

for source_target_size in "many_to_many_2^10" "many_to_many_2^14"
do
    mkdir -p "${OUTPUT_PATH}/${source_target_size}/custom/${SOURCES_TARGETS_FOLDER}"

    cargo run --release --bin experiment_command custom-experiment --name=rphast_bwd_dfs --data-path="${DATA_PAH}/${SOURCES_TARGETS_FOLDER}" --graph-path=$GRAPH_PATH --source-target-size=$source_target_size --num-iterations=$NUM_ITERATIONS > "${OUTPUT_PATH}/${source_target_size}/custom/${SOURCES_TARGETS_FOLDER}/rphast_bwd_dfs.csv"
    cargo run --release --bin experiment_command custom-experiment --name=rphast_bwd_ranks --data-path="${DATA_PAH}/${SOURCES_TARGETS_FOLDER}" --graph-path=$GRAPH_PATH --source-target-size=$source_target_size --num-iterations=$NUM_ITERATIONS > "${OUTPUT_PATH}/${source_target_size}/custom/${SOURCES_TARGETS_FOLDER}/rphast_bwd_ranks.csv"

    cargo run --release --bin experiment_command custom-experiment --name=simultaneous_buckets_bucket_size --data-path="${DATA_PAH}/${SOURCES_TARGETS_FOLDER}" --graph-path=$GRAPH_PATH --source-target-size=$source_target_size --num-iterations=$NUM_ITERATIONS > "${OUTPUT_PATH}/${source_target_size}/custom/${SOURCES_TARGETS_FOLDER}/simultaneous_buckets_bucket_size.csv"
    cargo run --release --bin experiment_command custom-experiment --name=buckets_bucket_size --data-path="${DATA_PAH}/${SOURCES_TARGETS_FOLDER}" --graph-path=$GRAPH_PATH --source-target-size=$source_target_size --num-iterations=$NUM_ITERATIONS > "${OUTPUT_PATH}/${source_target_size}/custom/${SOURCES_TARGETS_FOLDER}/buckets_bucket_size.csv"
    cargo run --release --bin experiment_command custom-experiment --name=buckets_no_pruning_bucket_size --data-path="${DATA_PAH}/${SOURCES_TARGETS_FOLDER}" --graph-path=$GRAPH_PATH --source-target-size=$source_target_size --num-iterations=$NUM_ITERATIONS > "${OUTPUT_PATH}/${source_target_size}/custom/${SOURCES_TARGETS_FOLDER}/buckets_no_pruning_bucket_size.csv"

    cargo run --release --bin experiment_command custom-experiment --name=rphast_restricted_graph_size --data-path="${DATA_PAH}/${SOURCES_TARGETS_FOLDER}" --graph-path=$GRAPH_PATH --source-target-size=$source_target_size --num-iterations=$NUM_ITERATIONS > "${OUTPUT_PATH}/${source_target_size}/custom/${SOURCES_TARGETS_FOLDER}/rphast_restricted_graph_size.csv"
done