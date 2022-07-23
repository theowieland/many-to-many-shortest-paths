use crate::types::*;
use std::cmp::max;
use std::cmp::Ordering;
use rand::thread_rng;
use rand::seq::SliceRandom;
use crate::data_structures::EmptyValidFlags;

/// converts an array of unique node ids to their respective ranks
///
/// example: ordering = [5, 2, 4, 0, 1, 3]. (NodeId 5 is first, NodeId 2 is second, ...)
/// nodeid 5 is the first and therefore receives the lowest rank (0) while node id 3 is the last and receives the highest rank
/// => ranks = [3, 4, 1, 5, 2, 0]
pub fn convert_ordering_to_ranks(ordering: &Vec<NodeId>) -> Ranks {
    let mut rank: Ranks = vec![0; ordering.len()];

    for (index, node) in ordering.iter().enumerate() {
        rank[*node as usize] = index;
    }

    rank
}

pub fn create_ordering_permutation(ordering: &Vec<NodeId>) -> Vec<NodeId> {
    let mut permutation: Vec<NodeId> = vec![0; ordering.len()];

    for (index, node) in ordering.iter().enumerate() {
        permutation[*node as usize] = index as NodeId;
    }

    permutation
}

/// creates a permutation for the given node ids so that permutated NodeId = rank
/// ranks = [3, 4, 1, 5, 2, 0]
/// => permutation = [3, 4, 1, 5, 2, 0]
pub fn create_rank_permutation(num_vertices: usize, rank: &Ranks) -> Vec<NodeId> {
    let mut permutation = vec![0; num_vertices];

    for node_index in 0..num_vertices {
        permutation[node_index] = rank[node_index] as NodeId;
    }

    permutation
}

/// creates a identity permutation that does not permutate the given NodeIds
pub fn create_identity_permutation(num_vertices: usize) -> Vec<NodeId> {
    let mut permutation = vec![0; num_vertices];

    for node_index in 0..num_vertices {
        permutation[node_index] = node_index as NodeId;
    }

    permutation
}

pub fn calculate_node_levels(first_edge: &EdgeIds, target_node: &NodeIds, rank: &Ranks) -> Vec<usize> {
    let num_vertices = first_edge.len() - 1;

    let mut levels: Vec<usize> = vec![0; num_vertices];
    let mut node_ids_ordered_by_ascending_rank: Vec<NodeId> = vec![0; num_vertices];

    for node_index in 0..num_vertices {
        node_ids_ordered_by_ascending_rank[rank[node_index]] = node_index as NodeId;
    }

    for index in 0..num_vertices {
        let current_node_id = node_ids_ordered_by_ascending_rank[index];

        let start = first_edge[current_node_id as usize] as usize;
        let end = first_edge[current_node_id as usize + 1] as usize;

        for upwards_neighbor in target_node[start..end].iter() {
            if rank[*upwards_neighbor as usize] > rank[current_node_id as usize] {
                levels[*upwards_neighbor as usize] = max(levels[*upwards_neighbor as usize], levels[current_node_id as usize] + 1);
            } 
        }
    }
    
    levels
}

pub fn create_level_permutation(first_edge: &EdgeIds, target_node: &NodeIds, weight: &Weights, rank: &Ranks) -> Vec<NodeId> {
    let num_vertices = first_edge.len() - 1;

    let levels = calculate_node_levels(first_edge, target_node, rank);

    let mut order = vec![0; num_vertices]; // stores node ids sorted by ascending level

    for node_index in 0..num_vertices {
        order[node_index] = node_index;
    }

    let mut node_ids_ordered_by_ascending_rank: Vec<NodeId> = vec![0; num_vertices];

    for node_index in 0..num_vertices {
        node_ids_ordered_by_ascending_rank[rank[node_index]] = node_index as NodeId;
    }

    let mut fwd_arclist = Vec::new();
    for (t, w) in target_node.iter().zip(weight.iter()) {
        fwd_arclist.push((*t, *w));
    }

    let dfs_order = depth_first_search(rank.len(), first_edge, &fwd_arclist, &mut EmptyValidFlags::new(rank.len()), &node_ids_ordered_by_ascending_rank);

    order.sort_by(|a, b| {
        let level_ordering = levels[*a as usize].partial_cmp(&levels[*b as usize]).unwrap();

        if level_ordering == Ordering::Equal {
            return dfs_order[*b as usize].partial_cmp(&dfs_order[*a as usize]).unwrap();
        }

        level_ordering
    });

    let mut permutation: Vec<NodeId> = vec![0; num_vertices];

    for node_index in 0..num_vertices {
        permutation[order[node_index]] = node_index as NodeId;
    }

    permutation
}

fn depth_first_search(
    num_vertices: usize,
    first_out: &Vec<EdgeId>, 
    arclist: &Vec<(NodeId, Weight)>, 
    discovered: &mut EmptyValidFlags, 
    start_nodes: &[NodeId]
) -> Vec<usize> {
    let mut stack: Vec<(NodeId, usize, usize)> = Vec::new(); // node index and the index to continue to scan the edges from and the end index
    let mut order_of_discovery = vec![0; num_vertices];

    discovered.reset();
    let mut discovered_counter = 0;

    for start_node in start_nodes {
        if discovered.is_valid(*start_node as usize) {
            continue;
        }

        let start = first_out[*start_node as usize] as usize;
        let end = first_out[*start_node as usize + 1] as usize;
        stack.push((*start_node, start, end));

        discovered.set_valid(*start_node as usize);
    
        while let Some((_last_element, start, end)) = stack.last() { //stack contains at least one element 
            let mut found_unvisited_node = false;
            
            for (target, _) in &arclist[*start..*end] {
                let len = stack.len();
                stack[len -1].1 += 1; // increase the start index for the next iterations

                if !discovered.is_valid(*target as usize) {
                    discovered.set_valid(*target as usize);
                    let target_start = first_out[*target as usize] as usize;
                    let target_end = first_out[*target as usize + 1] as usize;

                    stack.push((*target, target_start, target_end));
                    found_unvisited_node = true;
                    break;
                }
            }
            
            if !found_unvisited_node {
                if let Some((last_node, _, _)) = stack.pop() {
                    order_of_discovery[last_node as usize] = discovered_counter;
                    discovered_counter += 1;
                }
            }
        }
    }

    order_of_discovery
}

pub fn create_random_permutation(num_vertices: usize) -> Vec<NodeId> {
    let mut permutation = create_identity_permutation(num_vertices);
    permutation.shuffle(&mut thread_rng());

    permutation
}

pub fn apply_permutation(permutation: &Vec<NodeId>, first_out: &EdgeIds, head: &NodeIds, weight: &Weights, rank: &Ranks) -> (EdgeIds, NodeIds, Weights, Ranks) {
    let mut permutated_edges: Vec<Vec<(NodeId, Weight)>> = Vec::new();
    let mut permutated_ranks: Vec<usize> = vec![0; rank.len()];

    let num_vertices = first_out.len() - 1;

    for _node_index in 0..num_vertices {
        permutated_edges.push(Vec::new());
    }

    for node_index in 0..num_vertices {
        let start = first_out[node_index] as usize;
        let end = first_out[node_index + 1] as usize;

        let permutated_source = permutation[node_index];

        permutated_ranks[permutated_source as usize] = rank[node_index];

        for (target_node, edge_distance) in head[start..end].iter().zip(weight[start..end].iter()) {
            let permutated_target = permutation[*target_node as usize];

            permutated_edges[permutated_source as usize].push((permutated_target, *edge_distance));
        }
    }

    //reconstruct the permutated arrays
    let mut first_edge_index = 0;
    let mut permutated_first_out = vec![0; num_vertices + 1];
    let mut permutated_head = Vec::new();
    let mut permutated_weight = Vec::new();

    for node_index in 0..num_vertices {
        permutated_first_out[node_index] = first_edge_index;

        for (target_node, edge_distance) in &permutated_edges[node_index] {
            first_edge_index += 1;
            
            permutated_head.push(*target_node);
            permutated_weight.push(*edge_distance);
        }
    }

    permutated_first_out[num_vertices] = first_edge_index;

    (permutated_first_out, permutated_head, permutated_weight, permutated_ranks)
}
