set -e # abort complete script on a single error

DATA_PAH=$1
GRAPH_PATH=$2
OUTPUT_PATH=$3

NUM_ITERATIONS=$4
DEFAULT_BATCH_SIZE=16

# create output directories
mkdir -p "${OUTPUT_PATH}/rphast_variants" 

# run all many_to_many experiments
for experiment_algorithm in "rphast" "rphast_multiple_trees_fwd_rank" "rphast_multiple_trees_fwd_dijkstra" "rphast_multiple_trees_fwd_buckets" "simd_rphast"
do 
    cargo run --release --bin experiments many-to-many-experiment --algorithm=$experiment_algorithm --data-path=$DATA_PAH --graph-path=$GRAPH_PATH --num-iterations=$NUM_ITERATIONS --batch-size=$DEFAULT_BATCH_SIZE > "${OUTPUT_PATH}/rphast_variants/${experiment_algorithm}.csv"
done

read -rsp $'Press enter to exit...\n'