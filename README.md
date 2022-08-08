# Many to Many Shortest Paths
This repo contains the source code for my bachelor thesis "Extending, Recombining and Evaluating Contraction Hierarchy based Many-to-Many Shortest Path Algorithms".

The bachelor thesis can be found [here](https://i11www.iti.kit.edu/_media/teaching/theses/ba_wieland22.pdf).

## Note
The code is not exactly the same as the one I used for my expermients during my bachelor thesis, as some code segments were provided to me.
I have replaced these parts in order to avoid any copyright infringement.

## Usage
Graphs can be preprocessed/contracted in the following way:
The input file must be given in the same format as used by the DIMACS challenge (http://users.diag.uniroma1.it/challenge9/download.shtml)
Addtionally, contracted/exported graphs also contain an additional line that states the rank for each node id (line starts with r, an example of this can be found in the main.rs file) 

### Preprocessing/Contraction
```rust
let input_path = Path::new("path_to_input_file");
let output_path = Path::new("path_to_store_contracted_graph");

// load the graph data
match io::read_graph_data(&input_path) {
    Some((first_edge, target_node, weights, ranks)) => {
        println!("graph sucessfully loaded. num_nodes: {}, num_edges: {}", ranks.len(), target_node.len());

        // remove uneccessary arcs, i.e. loops and duplicate edges or edges between same nodes with different weights
        let (first_edge, target_node, weights) = remove_unecessary_arcs(&first_edge, &target_node, &weights);

        let graph = GraphArray::new(first_edge, target_node, weights);

        let (contraction_time, contracted_graph) = measure_time(|| {
            contract_graph(&graph)
        });
            
        println!("contraction done. time required: {} seconds", contraction_time.as_secs());

        // unwrap the contracted graph into its components
        let (contracted_first_edge, contracted_target_node, contracted_weight, contracted_rank) = contracted_graph;

        // rename node ids so that each node id is equal to its rank
        let num_nodes = contracted_first_edge.len() - 1;
        let rank_permutation = create_rank_permutation(num_nodes, &contracted_rank);

        let (contracted_first_edge, contracted_target_node, contracted_weight, contracted_rank) = apply_permutation(&rank_permutation, &contracted_first_edge, &contracted_target_node, &contracted_weight, &contracted_rank);

        // store contracted/preprocessed graph, this graph contains all the required shortcut arcs to efficently calculate all shortest paths
        io::export_graph_data(&output_path, contracted_first_edge, contracted_target_node, contracted_weight, contracted_rank);
    },
    None => println!("unable to read the given file!"),
}   
```

### Shortest Path Calculation
In order to caluclate shortest path distances a preprocessed/contracted input graph is required
```rust
match io::read_graph_data(&input_path) {
    let input_path = Path::new("path_to_preprocessed_graph");
    Some((first_edge, target_node, weight, rank)) => {
        println!("graph sucessfully loaded. num_nodes: {}, num_edges: {}", first_edge.len() - 1, target_node.len());

        let (fwd_first_edge, fwd_target_node, fwd_weight, bwd_first_edge, bwd_source_node, bwd_weight) = create_restricted_graph(&first_edge, target_node, &weight, &rank);

        let mut algorithm = SIMDRPHAST::new(&fwd_first_edge, &fwd_target_node, &fwd_weight, &bwd_first_edge, &bwd_source_node, &bwd_weight, &rank);
        
        // generate random set of source and target nodes
        let num_nodes = rank.len() as u32;
        let mut rng = rand::thread_rng();
        let sources: NodeIds = (0..100).map(|_| rng.gen_range(0..num_nodes)).collect();
        let targets: NodeIds = (0..100).map(|_| rng.gen_range(0..num_nodes)).collect();

        algorithm.calculate(&sources, &targets);
        

        println!("shortest path distances: {:?}", algorithm.get_distance_table().data);
    },
    None => println!("unable to read the given file!"),
}   
```
