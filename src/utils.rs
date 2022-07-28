use std::time::{Duration, Instant};

use crate::types::*;
use crate::data_structures::{ValidFlags, ArrayStructure};

pub mod binary_heap;
pub mod io;
pub mod node_picker;

/// creates a restricted graph by only adding edges whoose head node have a higher rank than the source node's rank
pub fn create_restricted_graph(first_out: &EdgeIds, head: &NodeIds, weight: &Weights, ranks: &Ranks) -> (EdgeIds, NodeIds, Weights, EdgeIds, NodeIds, Weights) {
    let num_vertices = first_out.len() - 1;
    
    let mut forward_edges: Vec<Vec<(NodeId, Weight)>> = vec![Vec::new(); num_vertices];
    let mut backward_edges: Vec<Vec<(NodeId, Weight)>> = vec![Vec::new(); num_vertices];

    for node_index in 0..num_vertices {
        for edge_index in first_out[node_index]..first_out[node_index + 1] {
            forward_edges[node_index].push((head[edge_index as usize], weight[edge_index as usize]));
            backward_edges[head[edge_index as usize] as usize].push((node_index as NodeId, weight[edge_index as usize]));
        }
    }

    for node_index in 0..num_vertices {
        forward_edges[node_index].sort_by_key(|(node_id, _)| *node_id); // sort by ascending node id
        backward_edges[node_index].sort_by_key(|(node_id, _)| *node_id); // sort by ascending node id
    }

    let (restricted_forward_first_out, restricted_forward_head, restricted_forward_weight) = convert_edges_to_arrays(&forward_edges, ranks);
    let (restricted_backward_first_out, restricted_backward_head, restricted_backward_weight) = convert_edges_to_arrays(&backward_edges, ranks);

    (restricted_forward_first_out, restricted_forward_head, restricted_forward_weight, restricted_backward_first_out, restricted_backward_head, restricted_backward_weight)
}

/// creates a full graph by combining the forward and backward search space into one search space
pub fn combine_restricted_graph(num_vertices: usize, fwd_first_out: &EdgeIds, fwd_head: &NodeIds, fwd_weight: &Weights, bwd_first_out: &EdgeIds, bwd_source: &NodeIds, bwd_weight: &Weights) -> (EdgeIds, NodeIds, Weights) {
    let mut edges: Vec<Vec<(NodeId, Weight)>> = vec![Vec::new(); num_vertices];

    // add all fwd edges
    for node_id in 0..num_vertices {
        let start = fwd_first_out[node_id] as usize;
        let end = fwd_first_out[node_id + 1] as usize;

        for (target_node, weight) in fwd_head[start..end].iter().zip(fwd_weight[start..end].iter()) {
            edges[node_id].push((*target_node, *weight));
        }
    }

    // add all bwd edges
    for node_id in 0..num_vertices {
        let start = bwd_first_out[node_id] as usize;
        let end = bwd_first_out[node_id + 1] as usize;

        for (source_node, weight) in bwd_source[start..end].iter().zip(bwd_weight[start..end].iter()) {
            edges[*source_node as usize].push((node_id as NodeId, *weight));
        }
    }

    let mut combined_first_out: Vec<EdgeId> = vec![0; num_vertices + 1];
    let mut combined_target: Vec<NodeId> = Vec::new();
    let mut combined_weight: Vec<Weight> = Vec::new();

    let mut first_edge_index = 0;

    for node_index in 0..num_vertices {
        combined_first_out[node_index] = first_edge_index;

        for (target_node, weight) in &edges[node_index] {
            combined_target.push(*target_node);
            combined_weight.push(*weight);
            first_edge_index += 1;
        }
    }

    combined_first_out[num_vertices] = first_edge_index;
 
    (combined_first_out, combined_target, combined_weight)
}

/// converts the given edges vec into another graph representation that used 3 different vecs.
/// ONLY INCLUDES EDGES TO ADJ NODES WITH HIGHER RANK
/// returns the first_out, head, weights array
fn convert_edges_to_arrays(edges: &Vec<Vec<(NodeId, Weight)>>, ranks: &Vec<usize>) -> (Vec<EdgeId>, Vec<NodeId>, Vec<Weight>) {
    let num_vertices = edges.len();

    let mut result_first_edge = vec![0; num_vertices + 1];
    let mut result_head = Vec::new();
    let mut result_weight = Vec::new();

    let mut first_edge_index = 0;

    for node_index in 0..num_vertices {
        result_first_edge[node_index] = first_edge_index;

        for (adj_node, weight) in &edges[node_index] {
            if ranks[*adj_node as usize] > ranks[node_index] {
                result_head.push(*adj_node);
                result_weight.push(*weight);

                first_edge_index += 1;
            }
        }
    }

    result_first_edge[num_vertices] = first_edge_index;

    (result_first_edge, result_head, result_weight)
}


/// performs a depth first search from the given vec of start nodes and returns all discovered nodes in order of discovery
pub fn depth_first_search(
    first_out: &Vec<EdgeId>, 
    arclist: &Vec<(NodeId, Weight)>, 
    translation_table: &mut ValidFlags<NodeId>, 
    start_nodes: &[NodeId],
    order_of_discovery: &mut Vec<NodeId>
) {
    let mut stack: Vec<(NodeId, usize, usize)> = Vec::new(); // node index and the index to continue to scan the edges from and the end index
    order_of_discovery.clear();

    translation_table.reset();

    for start_node in start_nodes {
        if translation_table.is_valid(*start_node as usize) {
            continue;
        }

        let start = first_out[*start_node as usize] as usize;
        let end = first_out[*start_node as usize + 1] as usize;
        stack.push((*start_node, start, end));

        translation_table.set(*start_node as usize, 0); // initially set to local id to 0
    
        while let Some((_last_element, start, end)) = stack.last() { //stack contains at least one element 
            let mut found_unvisited_node = false;
            
            for (target, _) in &arclist[*start..*end] {
                let len = stack.len();
                stack[len -1].1 += 1; // increase the start index for the next iterations

                if !translation_table.is_valid(*target as usize) { // node has not been previously visited by another branch
                    translation_table.set(*target as usize, 0); // initially assign id 0 to newly discovered node
                    let target_start = first_out[*target as usize] as usize;
                    let target_end = first_out[*target as usize + 1] as usize;

                    stack.push((*target, target_start, target_end));
                    found_unvisited_node = true;
                    break;
                }
            }
            
            if !found_unvisited_node {
                if let Some((last_node, _, _)) = stack.pop() {
                    translation_table.set(last_node as usize, order_of_discovery.len() as NodeId); //assign local translated id to the node
                    order_of_discovery.push(last_node);
                }
            }
        }
    }
}

pub fn measure_time<F: FnOnce()>(function: F) -> Duration {
    let start = Instant::now();
    let result = function();
    
    start.elapsed()
}