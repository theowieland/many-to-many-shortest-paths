use std::path::Path;

use rand::prelude::SliceRandom;
use shortest_path_algorithms::{utils::io, graph_representation::{GraphArray, Graph}, contraction_hierarchies::contract_graph};
use rand::thread_rng;

fn main() {
    // example code 

    // read graph data
    let input_path = Path::new("/users/theo/Downloads/USA-road-t.CTR.gr");
    let output_path = Path::new("/Users/theo/Downloads/graph_output.CTR.gr");

    println!("reading the graph from {:?}", input_path);

    match io::read_graph_data(&input_path) {
        Some((first_edge, target_node, weights, ranks)) => {
            let graph = GraphArray::new(first_edge, target_node, weights);
            let mut contraction_order: Vec<u32> = (0..(graph.node_ids().len() as u32)).collect();
            contraction_order.shuffle(&mut thread_rng());

            contract_graph(&graph, &contraction_order);

            //io::export_graph_data(&output_path, first_edge, target_node, weights);

        },
        None => println!("unable to read the given file!"),
    }
}