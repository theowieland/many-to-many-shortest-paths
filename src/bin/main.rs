use std::path::Path;

use shortest_path_algorithms::utils::io;

fn main() {
    // example code 

    // read graph data
    let path = Path::new("/Users/theo/Downloads/USA-road-t.CTR.grs");
    println!("reading the graph from {:?}", path);

    match io::read_graph_data(&path) {
        Some((first_edge, target_nodes, weights)) => {

        },
        None => println!("unable to read the given file!"),
    }
}