use std::path::Path;
use shortest_path_algorithms::{utils::{io, measure_time, remove_unecessary_arcs}, graph_representation::GraphArray, contraction_hierarchies::contract_graph, graph_permutations::{create_rank_permutation, apply_permutation}};

fn main() {
    contraction_example();
}

fn contraction_example() {
    // read graph data
    let input_path = Path::new("/users/theo/Downloads/USA-road-t.CAL.gr");
    let output_path = Path::new("/Users/theo/Downloads/contracted_graph_output.CAL.gr");

    println!("reading the graph from {:?}", input_path);

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