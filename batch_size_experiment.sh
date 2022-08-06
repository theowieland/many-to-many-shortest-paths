set -e # abort complete script on a single error

DATA_PAH=$1
GRAPH_PATH=$2
OUTPUT_PATH=$3
NUM_ITERATIONS=$4

SOURCES_TARGETS_FOLDER="different_balls"

# run batch_size experiments 
for algorithm in "rphast_multiple_trees_fwd_buckets" "lazy_batched_advanced_buckets" "rphast_multiple_trees_fwd_rank"
do
    for source_target_size in "many_to_many_2^10"
    do
        # create output directory
        mkdir -p "${OUTPUT_PATH}/${source_target_size}/batch/${SOURCES_TARGETS_FOLDER}/${algorithm}"

        for experiment_batch_size in 2 4 8 16 32 64 128 256 512 1024
        do
            # zero pad batch size
            printf -v zero_padded_batch_size "%05d" $experiment_batch_size

            cargo run --release --bin experiment_command experiment --algorithm=$algorithm --data-path="${DATA_PATH}/${SOURCES_TARGETS_FOLDER}" --graph-path=$GRAPH_PATH --source-target-size=$source_target_size --num-iterations=$NUM_ITERATIONS --batch-size=$experiment_batch_size > "${OUTPUT_PATH}/${source_target_size}/batch/${SOURCES_TARGETS_FOLDER}/${algorithm}/${zero_padded_batch_size}.csv"
        done
    done
done

# run batch_size experiments 
for algorithm in "rphast_multiple_trees_fwd_buckets" "lazy_batched_advanced_buckets" "rphast_multiple_trees_fwd_rank"
do
    for source_target_size in "many_to_many_2^14"
    do
        # create output directory
        mkdir -p "${OUTPUT_PATH}/${source_target_size}/batch/${SOURCES_TARGETS_FOLDER}/${algorithm}"

        for experiment_batch_size in 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384
        do
            # zero pad batch size
            printf -v zero_padded_batch_size "%05d" $experiment_batch_size

            cargo run --release --bin experiment_command experiment --algorithm=$algorithm --data-path="${DATA_PATH}/${SOURCES_TARGETS_FOLDER}" --graph-path=$GRAPH_PATH --source-target-size=$source_target_size --num-iterations=$NUM_ITERATIONS --batch-size=$experiment_batch_size > "${OUTPUT_PATH}/${source_target_size}/batch/${SOURCES_TARGETS_FOLDER}/${algorithm}/${zero_padded_batch_size}.csv"
        done
    done
done


read -rsp $'Press enter to exit...\n'