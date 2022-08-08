use crate::{types::{NodeId, EdgeId, Weight, INFINITY}, graph_algorithms::DijkstraState, utils::{binary_heap::MinBinaryHeap, data_structures::{ResettableArray, ValidFlags, ArrayStructure}}};

use super::many_to_many_utils::RankDijkstraState;

// perform rphast forward search for multiple sources
pub fn rphast_forward_rank(
    batched_sources: &[NodeId],
    queue: &mut MinBinaryHeap<RankDijkstraState>,
    fwd_first_edge: &Vec<EdgeId>,
    fwd_arclist: &Vec<(NodeId, Weight)>,
    ranks: &Vec<usize>,
    fwd_distances: &mut Vec<Option<Box<[Weight]>>>,
    batch_size: usize,
    bwd_restricted_nodes_translation: &ValidFlags<NodeId>,
    bwd_restricted_distances: &mut Vec<Weight>
) {
    // start the fwd search by adding the batched sources with their initial 0 distance
    for (batch_source_index, source) in batched_sources.iter().enumerate() {
        queue.insert(RankDijkstraState {rank: ranks[*source as usize], node_id: *source});

        // initialize fwd distance
        if fwd_distances[*source as usize] == None {
            fwd_distances[*source as usize] = Some((0..batch_size).map(|_| INFINITY as Weight).collect());
        }
    
        fwd_distances[*source as usize].as_mut().unwrap()[batch_source_index] = 0;
    } 

    while let Some(RankDijkstraState {rank: _, node_id: current_node}) = queue.pop() { // pop current node with the lowest rank
        unsafe {
            let start = *fwd_first_edge.get_unchecked(current_node as usize) as usize;
            let end = *fwd_first_edge.get_unchecked(current_node as usize + 1) as usize;
        
            let (before_current_distances, current_and_after_current) = fwd_distances.split_at_mut(current_node as usize);
            let (current_distances, after_current_distances) = current_and_after_current.split_at_mut(1);

            let current_distances = &current_distances[0];

            for (target_node, edge_distance) in fwd_arclist.get_unchecked(start..end) {
                let target_distances = if *target_node > current_node {
                    after_current_distances.get_unchecked_mut((*target_node - current_node - 1) as usize)
                }
                else {
                    before_current_distances.get_unchecked_mut(*target_node as usize)
                };
            
                if *target_distances == None { // distances have not been initialized
                    *target_distances = Some((0..batch_size).map(|_| INFINITY as Weight).collect());
                }

                // reduce all distances by possibly using the currently discovered edges
                for batch_source_index in 0..batched_sources.len() {
                    let current_distance = current_distances.as_ref().unwrap().get_unchecked(batch_source_index);
                    let target_distance = target_distances.as_mut().unwrap().get_unchecked_mut(batch_source_index);

                    if edge_distance + current_distance < *target_distance {
                        *target_distance = edge_distance + current_distance;
                    }
                }
            
                if !queue.contains_unique_index(*target_node as usize) {
                    queue.insert(RankDijkstraState {rank: ranks[*target_node as usize], node_id: *target_node});
                }   
            }

            // copy distances to the bwd_restricted_distance array
            if bwd_restricted_nodes_translation.is_valid(current_node as usize) { // the currently visisted node is also in the target set
                let bwd_restricted_node_index = bwd_restricted_nodes_translation.get_unchecked(current_node as usize) as usize;

                bwd_restricted_distances[(bwd_restricted_node_index * batch_size)..((bwd_restricted_node_index + 1) * batch_size)].copy_from_slice(&current_distances.as_ref().unwrap());
            }

            // clear old distances that are never accessed again (rank order)
            fwd_distances[current_node as usize] = None;
        }
    }
}

// perform rphast forward search for multiple sources
pub fn rphast_forward_dijkstra(
    batched_sources: &[NodeId],
    queue: &mut MinBinaryHeap<DijkstraState>,
    fwd_first_edge: &Vec<EdgeId>,
    fwd_arclist: &Vec<(NodeId, Weight)>,
    fwd_distances: &mut ResettableArray<Weight>,
    batch_size: usize,
    bwd_restricted_nodes_translation: &ValidFlags<NodeId>,
    bwd_restricted_distances: &mut Vec<Weight>
) {
    for (source_index, source_node) in batched_sources.iter().enumerate() {
        fwd_distances.reset();
        queue.insert(DijkstraState {distance: 0, node_id: *source_node});
        fwd_distances.set(*source_node as usize, 0);

        while let Some(DijkstraState {distance: current_distance, node_id: current_node}) = queue.pop() {
            unsafe {
                let start = *fwd_first_edge.get_unchecked(current_node as usize) as usize;
                let end = *fwd_first_edge.get_unchecked(current_node as usize + 1) as usize;
    
                for (target_node, edge_distance) in fwd_arclist.get_unchecked(start..end) {
                    let target_distance = fwd_distances[*target_node as usize];

                    if current_distance + *edge_distance < target_distance {
                        fwd_distances.set(*target_node as usize, current_distance + *edge_distance);

                        queue.insert_or_decrease(DijkstraState {distance: current_distance + *edge_distance, node_id: *target_node});
                    }
                }

                if bwd_restricted_nodes_translation.is_valid(current_node as usize) { // the currently visisted node is also in the target set
                    let bwd_restricted_node_index = bwd_restricted_nodes_translation.get_unchecked(current_node as usize) as usize;
    
                    bwd_restricted_distances[bwd_restricted_node_index * batch_size + source_index] = current_distance;
                }
            }
        }
    }
}

// perform rphast downward search for multiple sources
pub fn rphast_downward(
    bwd_restricted_nodes: &Vec<NodeId>, // all nodes in the restricted search space, ordered by descending rank
    bwd_restricted_first_edge: &Vec<EdgeId>,
    bwd_restricted_arclist: &Vec<(NodeId, EdgeId)>, 
    bwd_restricted_distances: &mut Vec<Weight>, // flat 2d distance vec to store all distances
    batch_size: usize,
    batched_sources: &[NodeId]
) {
    unsafe {
        for node_index in 0..bwd_restricted_nodes.len() {
            let start = *bwd_restricted_first_edge.get_unchecked(node_index) as usize;
            let end = *bwd_restricted_first_edge.get_unchecked(node_index + 1) as usize;
    
            let (before_current_slice, current_and_after_slice) = bwd_restricted_distances.split_at_mut_unchecked(node_index as usize * batch_size);
            let (current_distances, _after_current_slice) = current_and_after_slice.split_at_mut_unchecked(batch_size);
    
            // iterate over all incoming arcs to the current_node (source_node is a local id given in the range of the bwd_restricted search space)
            for (arc_source_node, edge_distance) in bwd_restricted_arclist.get_unchecked(start..end) {
                let old_start = (*arc_source_node as usize) * batch_size;
                let old_end = old_start + batch_size;
                
                // since local ids are created from top to low nodes with low to high ids, we can assume that arc_source_node < node_index -> only access before_current_slice
                let old_distances = before_current_slice.get_unchecked_mut(old_start..old_end);
    
                // iterate over all source nodes in the current batch
                for batch_source_index in 0..batched_sources.len() {
                    let arc_source_distance = old_distances.get_unchecked(batch_source_index);
                    let old_distance = current_distances.get_unchecked_mut(batch_source_index);
    
                    if *old_distance > *edge_distance + arc_source_distance {
                        *old_distance = *edge_distance + arc_source_distance;    
                    }
                }
            }
        }
    }
}