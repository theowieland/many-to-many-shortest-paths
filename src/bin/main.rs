use std::path::Path;
use rand::Rng;
use shortest_path_algorithms::{utils::{io, measure_time, remove_unecessary_arcs, create_restricted_graph}, graph_representation::GraphArray, contraction_hierarchies::contract_graph, graph_permutations::{create_rank_permutation, apply_permutation}, many_to_many::{rphast_simd::SIMDRPHAST, many_to_many_algorithm::ManyToManyAlgorithm}, types::NodeIds};

fn main() {
    let input_path = Path::new("/users/theo/Downloads/USA-road-t.CAL.gr");
    let output_path = Path::new("users/theo/Downloads/contracted_graph_output.CAL.gr");

    contraction_example(&input_path, &output_path);
    shortest_path_example(&input_path);
}

fn contraction_example(input_path: &dyn AsRef<Path>, output_path: &dyn AsRef<Path>) {
    // read graph data
    match io::read_graph_data(&input_path) {
        Some((first_edge, target_node, weights, _ranks)) => {
            println!("graph sucessfully loaded. num_nodes: {}, num_edges: {}", first_edge.len() - 1, target_node.len());

            // remove uneccessary arcs
            let (first_edge, target_node, weights) = remove_unecessary_arcs(&first_edge, &target_node, &weights);

            let graph = GraphArray::new(first_edge, target_node, weights);

            let (contraction_time, contracted_graph) = measure_time(|| {
                contract_graph(&graph)
            });
            
            println!("contraction done. time required: {} seconds", contraction_time.as_secs());

            let (contracted_first_edge, contracted_target_node, contracted_weight, contracted_rank) = contracted_graph;

            // rename node ids so that each node id is equal to its rank
            let num_nodes = contracted_first_edge.len() - 1;
            let rank_permutation = create_rank_permutation(num_nodes, &contracted_rank);

            let (contracted_first_edge, contracted_target_node, contracted_weight, contracted_rank) = apply_permutation(&rank_permutation, &contracted_first_edge, &contracted_target_node, &contracted_weight, &contracted_rank);

            // store result graph
            io::export_graph_data(&output_path, contracted_first_edge, contracted_target_node, contracted_weight, contracted_rank);
        },
        None => println!("unable to read the given file!"),
    }   
}

fn shortest_path_example(input_path: &dyn AsRef<Path>) {
    match io::read_graph_data(&input_path) {
        Some((first_edge, target_node, weight, rank)) => {
            println!("graph sucessfully loaded. num_nodes: {}, num_edges: {}", first_edge.len() - 1, target_node.len());

            let (fwd_first_edge, fwd_target_node, fwd_weight, bwd_first_edge, bwd_source_node, bwd_weight) = create_restricted_graph(&first_edge, &target_node, &weight, &rank);

            let mut algorithm = SIMDRPHAST::new(&fwd_first_edge, &fwd_target_node, &fwd_weight, &bwd_first_edge, &bwd_source_node, &bwd_weight, &rank);
            
            let num_nodes = rank.len() as u32;
            let mut rng = rand::thread_rng();
            let sources: NodeIds = (0..100).map(|_| rng.gen_range(0..num_nodes)).collect();
            let targets: NodeIds = (0..100).map(|_| rng.gen_range(0..num_nodes)).collect();

            algorithm.calculate(&sources, &targets);
            println!("shortest path distances: {:?}", algorithm.get_distance_table().data);
        },
        None => println!("unable to read the given file!"),
    }   
}