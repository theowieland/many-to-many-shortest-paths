use std::cmp::min;
use crate::types::*;
use crate::utils::binary_heap::{HeapElement, MinBinaryHeap};
use crate::utils::data_structures::{ValidFlags, ArrayStructure};
use crate::utils::depth_first_search;

#[derive(Copy, Clone, Eq, PartialEq, Debug, Ord, PartialOrd)]
pub struct RankDijkstraState {
    pub rank: usize,
    pub node_id: NodeId
}

impl HeapElement for RankDijkstraState {

    fn unique_index(&self) -> usize {
        self.node_id as usize
    }
}

pub fn fill_restricted_search_space_decreasing_rank(
    initial_nodes: &[NodeId], // the nodes to start the search from
    first_edge: &Vec<NodeId>, // the full (unrestricted) first edge graph
    arclist: &Vec<(NodeId, Weight)>, // the full (unrestricted) arclist graph
    translation_table: &mut ValidFlags<NodeId>, // stores translations from unrestricted to restricted ids
    restricted_nodes: &mut Vec<NodeId>, // stores node ids in order of discovery (unrestricted ids)
    restricted_first_edge: &mut Vec<NodeId>,
    restricted_arclist: &mut Vec<(NodeId, Weight)>
) {
    // fill the discovered nodes array with all discovered nodes
    depth_first_search(first_edge, arclist, translation_table, initial_nodes, restricted_nodes);

    // fill the provided arrays with the restricted search space data
    fill_restricted_arrays(restricted_nodes, translation_table, first_edge, arclist, restricted_first_edge, restricted_arclist);
}

pub fn fill_restricted_search_space_increasing_rank(
    initial_nodes: &[NodeId], //the nodes to start the search from
    first_edge: &Vec<NodeId>, //the full (unrestricted) first edge graph
    arclist: &Vec<(NodeId, Weight)>, //the full (unrestricted) arclist graph
    ranks: &Vec<usize>, //rank for each node
    queue: &mut MinBinaryHeap<RankDijkstraState>, //queue uesed to calculate the restricted search space
    translation_table: &mut ValidFlags<NodeId>, //stores translations from unrestricted to restricted ids
    restricted_nodes: &mut Vec<NodeId>, //stores node ids in order of discovery (unrestricted ids)
    restricted_first_edge: &mut Vec<NodeId>,
    restricted_arclist: &mut Vec<(NodeId, Weight)>
) {
    translation_table.reset();

    restricted_nodes.clear();
    for initial_node in initial_nodes {
        queue.insert(RankDijkstraState {rank: ranks[*initial_node as usize], node_id: *initial_node});
    }

    //discover all required nodes
    while let Some(RankDijkstraState { rank: _, node_id: current_node }) = queue.pop() {
        translation_table.set(current_node as usize, restricted_nodes.len() as NodeId);
        restricted_nodes.push(current_node);

        let start = first_edge[current_node as usize] as usize;
        let end = first_edge[current_node as usize + 1] as usize;

        for (target_node, _) in &arclist[start..end] {
            queue.insert_or_decrease(RankDijkstraState {rank: ranks[*target_node as usize], node_id: *target_node});
        }
    }

    //fill the provided arrays with the restricted search space data
    fill_restricted_arrays(restricted_nodes, translation_table, first_edge, arclist, restricted_first_edge, restricted_arclist);
}

fn fill_restricted_arrays(discovered_nodes: &Vec<NodeId>, translation_table: &ValidFlags<NodeId>, first_edge: &Vec<EdgeId>, arclist: &Vec<(NodeId, Weight)>, restricted_first_edge: &mut Vec<NodeId>, restricted_arclist: &mut Vec<(NodeId, Weight)>) {
    restricted_first_edge.resize(discovered_nodes.len() + 1, 0);
    restricted_arclist.clear();

    let mut first_edge_index = 0;

    //fill restricted_first_edge and restricted_arclist with the restricted search space data
    for (node_index, node) in discovered_nodes.iter().enumerate() {
        let start = first_edge[*node as usize] as usize;
        let end = first_edge[*node as usize + 1] as usize;

        restricted_first_edge[node_index] = first_edge_index;

        for (source_node, weight) in &arclist[start..end] {
            restricted_arclist.push((translation_table[*source_node as usize], *weight));
            first_edge_index += 1;
        }
    }

    restricted_first_edge[discovered_nodes.len()] = first_edge_index;
}

pub fn convert_to_arclist(head: &NodeIds, weight: &Weights) -> Vec<(NodeId, Weight)> {
    let mut arclist = Vec::new();
    for (node_id, weight) in head.iter().zip(weight.iter()) {
        arclist.push((*node_id, *weight));
    }

    arclist
}


/// split the given nodes vec into equal batch_size slices and call the exectute method for each batch
pub fn batch_nodes<F>(batch_size: usize, nodes: &[NodeId], mut execute: F) where F: FnMut(usize, &[NodeId]) {
    let num_batches = (nodes.len() / batch_size) + if nodes.len() % batch_size != 0 {1} else {0};

    for batch in 0..num_batches {
        let batch_first_source_index = batch * batch_size; // index of first source node in the current batch
        let batch_last_source_index = min((batch + 1) * batch_size, nodes.len()); // index of the last source node in the current batch

        execute(batch, &nodes[batch_first_source_index..batch_last_source_index])
    }
}

pub fn iterate_nodes_by_rank<F>(
    initial_nodes: &[NodeId], 
    queue: &mut MinBinaryHeap<RankDijkstraState>, 
    first_out: &Vec<EdgeId>,
    arclist: &Vec<(NodeId, Weight)>,
    ranks: &Vec<usize>,
    mut process_arcs: F
) where F: FnMut(NodeId, &[(NodeId, Weight)]) -> () {
    for node_id in initial_nodes {
        queue.insert(RankDijkstraState {rank: ranks[*node_id as usize], node_id: *node_id});
    }

    while let Some(RankDijkstraState {rank: _, node_id: current_node}) = queue.pop() {
        unsafe {
            let start = *first_out.get_unchecked(current_node as usize) as usize;
            let end = *first_out.get_unchecked(current_node as usize + 1) as usize;
            let current_arclist = arclist.get_unchecked(start..end);
    
            process_arcs(current_node, current_arclist);
    
            for (target_node, _edge_weight) in current_arclist {
                if !queue.contains_unique_index(*target_node as usize) {
                    queue.insert(RankDijkstraState {rank: *ranks.get_unchecked(*target_node as usize), node_id: *target_node});
                }
            }
        }
    }
}