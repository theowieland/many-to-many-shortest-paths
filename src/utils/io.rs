use std::{path::Path, fs::File};
use std::io::{BufRead, BufReader, Write};
use crate::types::*;

pub fn read_graph_data(path: &dyn AsRef<Path>) -> Option<(EdgeIds, NodeIds, Weights, Ranks)>{
    let input_file = File::open(path);
    if let Ok(file) = input_file {
        let reader = BufReader::new(file);

        let mut arcs: Vec<Vec<(NodeId, Weight)>> = Vec::new();
        let mut num_nodes = 0;
        let mut rank: Ranks = Vec::new();

        for line in reader.lines() {
            if let Ok(string_line) = line {
                let split = string_line.split(" ").collect::<Vec<&str>>();

                if is_graph_size_line(&string_line) {
                    num_nodes = split[2].parse().unwrap();

                    arcs.resize(num_nodes, Vec::new());
                    rank.resize(num_nodes, 0);

                }
                else if is_arc_line(&string_line) {
                    let source_node: NodeId = split[1].parse().unwrap();
                    let target_node: NodeId = split[2].parse().unwrap();
                    let arc_weight: Weight = split[3].parse().unwrap();

                    // reduce source and target node ids by one so that they are zero based
                    arcs[(source_node - 1) as usize].push((target_node - 1, arc_weight));
                }
                else if is_rank_line(&string_line) {
                    let node_id: NodeId = split[1].parse().unwrap();
                    let node_rank: usize = split[2].parse().unwrap();

                    rank[node_id as usize] = node_rank;
                }
            }
        }

        // convert arcs to adjacency array representation
        let (first_out, target_node, weight) = convert_vec_to_adjacency_array(num_nodes, &arcs);
        return Some((first_out, target_node, weight, rank));
    }

    None
}

pub fn export_graph_data(path: &dyn AsRef<Path>, first_edge: EdgeIds, target_node: NodeIds, weight: Weights, rank: Ranks) {
    let output_file = File::create(path);

    if let Ok(mut file) = output_file {
        writeln!(&mut file, "p sp {} {}", first_edge.len(), target_node.len()).unwrap();

        for node_id in 0..first_edge.len() {
            for edge_id in first_edge[node_id]..first_edge[node_id + 1] {
                writeln!(&mut file, "a {} {} {}", node_id, target_node[edge_id as usize], weight[edge_id as usize]).unwrap();
            }
        }

        for node_id in 0..first_edge.len() {
            writeln!(&mut file, "r {} {}", node_id, rank[node_id]).unwrap();
        }
    }
}

fn is_graph_size_line(line: &String) -> bool {
    line.starts_with("p sp")
}

fn is_arc_line(line: &String) -> bool {
    line.starts_with("a")
}

fn is_rank_line(line: &String) -> bool {
    line.starts_with("r")
}

fn convert_vec_to_adjacency_array(num_vertices: usize, arcs: &Vec<Vec<(NodeId, Weight)>>) -> (EdgeIds, NodeIds, Weights) {
    let mut first_edge: EdgeIds = vec![0; num_vertices + 1];
    let mut arc_target: NodeIds = Vec::new();
    let mut arc_weight: Weights = Vec::new();

    let mut first_edge_index = 0;

    for node_index in 0..num_vertices {
        first_edge[node_index] = first_edge_index;

        for (adj_node, weight) in &arcs[node_index] {
            arc_target.push(*adj_node);
            arc_weight.push(*weight);
            first_edge_index += 1;
        }
    }

    first_edge[num_vertices] = first_edge_index;

    (first_edge, arc_target, arc_weight)
}