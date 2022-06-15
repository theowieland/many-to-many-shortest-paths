use crate::bucket_containers::{BucketContainer, PreparationAccessBucketContainer};
use crate::types::*;
use crate::data_structures::Matrix;
use crate::utils::binary_heap::MinBinaryHeap;
use std::cmp::{self, min};
use crate::many_to_many::many_to_many_algorithm::ManyToManyAlgorithm;
use crate::many_to_many::many_to_many_utils::*;
use super::bucket::query_no_stopping;

pub struct AdvancedBucketManyToMany {
    ranks: Vec<usize>,
    fwd_first_edge: Vec<EdgeId>,
    fwd_arclist: Vec<(NodeId, Weight)>,
    bwd_first_edge: Vec<EdgeId>,
    bwd_arclist: Vec<(NodeId, Weight)>,
    fwd_distances: Vec<Weight>,
    fwd_queue: MinBinaryHeap<RankDijkstraState>,
    bwd_distances: Vec<Weight>,
    bwd_terminals: Vec<usize>,
    bwd_queue: MinBinaryHeap<RankDijkstraState>,
    bwd_buckets: BucketContainer<(usize, Weight)>,
    down_nodes: Vec<Vec<(NodeId, Weight)>>,
    up_nodes: Vec<Vec<(NodeId, Weight)>>,
    modified_up_nodes: Vec<NodeId>,
    distance_table: Matrix<Weight>,
}

impl AdvancedBucketManyToMany {

    pub fn new(
        fwd_first_edge: &Vec<EdgeId>, 
        fwd_head: &Vec<NodeId>, 
        fwd_weight: &Vec<Weight>, 
        bwd_first_edge: &Vec<EdgeId>, 
        bwd_head: &Vec<NodeId>, 
        bwd_weight: &Vec<Weight>, 
        ranks: &Vec<usize>
    ) -> Self {
        let num_vertices = fwd_first_edge.len() - 1;

        AdvancedBucketManyToMany {
            ranks: ranks.to_vec(),
            fwd_first_edge: fwd_first_edge.to_owned(),
            fwd_arclist: convert_to_arclist(fwd_head, fwd_weight),
            bwd_first_edge: bwd_first_edge.to_owned(),
            bwd_arclist: convert_to_arclist(bwd_head, bwd_weight),
            fwd_distances: vec![INFINITY; num_vertices],
            fwd_queue: MinBinaryHeap::new(num_vertices),
            bwd_distances: Vec::new(),
            bwd_terminals: Vec::new(),
            bwd_queue: MinBinaryHeap::new(num_vertices),
            bwd_buckets: BucketContainer::new(num_vertices),
            down_nodes: vec![Vec::new(); num_vertices],
            up_nodes: vec![Vec::new(); num_vertices],
            modified_up_nodes: Vec::new(),
            distance_table: Matrix::empty()
        }
    }

    /// initializes the buckets of all nodes discovered by a backward search starting from the given target nodes
    /// performs a full backward search for all targets in one backward search
    pub fn select_targets(&mut self, targets: &[NodeId]) {      
        populate_buckets_pruing_fwd(
            &mut self.bwd_buckets, 
            targets, 
            &mut self.bwd_terminals,
            &mut self.bwd_distances, 
            &mut self.down_nodes,
            &mut self.up_nodes,
            &mut self.modified_up_nodes,
            &mut self.bwd_queue, 
            &self.bwd_first_edge, 
            &self.bwd_arclist, 
            &self.fwd_first_edge,
            &self.fwd_arclist,
            &self.ranks,
            |_settled_node, _terminals, _distance| {}
        );
    }

    pub fn query(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        query_no_stopping(
            sources, 
            targets, 
            &self.fwd_first_edge, 
            &self.fwd_arclist, 
            &self.ranks, 
            &mut self.fwd_distances, 
            &mut self.fwd_queue, 
            &mut self.bwd_buckets, 
            &mut self.distance_table
        );
    }

    pub fn get_average_number_of_bucket_entries(&mut self) -> (usize, usize) {
        let mut total_bucket_entry_count = 0;

        for modified_bucket in self.bwd_buckets.get_modified_buckets() {
            total_bucket_entry_count += self.bwd_buckets.get_bucket(*modified_bucket as usize).len();
        }
    
        (total_bucket_entry_count, self.bwd_buckets.get_modified_buckets().len())
    }

    pub fn print_bucket(&mut self) {
        for modified_bucket in self.bwd_buckets.get_modified_buckets() {
            println!("{} - {:?}", modified_bucket, self.bwd_buckets.get_bucket(*modified_bucket));
        }
    }
}

impl ManyToManyAlgorithm for AdvancedBucketManyToMany {
    
    fn initialize(&mut self, _sources: &[NodeId], _targets: &[NodeId]) {
        //no initialization required
    }

    fn select(&mut self, _sources: &[NodeId], targets: &[NodeId]) {
        self.select_targets(targets);
    }

    fn query(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        self.query(sources, targets);
    }

    fn get_distance_table(&mut self) -> &Matrix<Weight> {
        &self.distance_table
    }
}

/// populates all the buckets of the associated nodes for a given backward search space
pub fn initialize_bwd_buckets(
    bucket_container: &mut BucketContainer<(usize, Weight)>, 
    bucket_entry_nodes: &[NodeId], 
    terminals: &mut Vec<usize>,
    distances: &mut Vec<Weight>,
    down_nodes: &mut Vec<Vec<(NodeId, Weight)>>, 
    queue: &mut MinBinaryHeap<RankDijkstraState>, 
    bwd_first_edge: &Vec<NodeId>, 
    bwd_arclist: &Vec<(NodeId, Weight)>, 
    ranks: &Vec<usize>
) {
    bucket_container.reset();

    for (target_index, target) in bucket_entry_nodes.iter().enumerate() {
        if !queue.contains_unique_index(*target as usize) {
            queue.insert(RankDijkstraState {rank: ranks[*target as usize], node_id: *target});   
        }

        bucket_container.add_bucket_entry(*target as usize, (target_index, 0));
    }
    
    distances.resize(bucket_entry_nodes.len(), INFINITY);
    distances.iter_mut().for_each(|entry| *entry = INFINITY);

    while let Some(RankDijkstraState {rank: _, node_id: current_node}) = queue.pop() {
        let start_index = bwd_first_edge[current_node as usize] as usize;
        let end_index = bwd_first_edge[current_node as usize + 1] as usize;

        for (source_node, source_distance) in &bwd_arclist[start_index..end_index] {
            down_nodes[*source_node as usize].push((current_node, *source_distance));

            if !queue.contains_unique_index(*source_node as usize) {
                queue.insert(RankDijkstraState {rank: ranks[*source_node as usize], node_id: *source_node});
            }
        }
        
        //iterate over all nodes that the current_node is connected to (edge: current_node -> down_node)
        for (down_node, down_distance) in &down_nodes[current_node as usize] {
            for (target_index, target_distance) in bucket_container.get_bucket(*down_node as usize) {
                let current_distance = &mut distances[*target_index];
                
                if *current_distance == INFINITY {
                    terminals.push(*target_index);
                }

                *current_distance = cmp::min(*current_distance, target_distance + down_distance);
            }
        }

        terminals.iter().for_each(|index| {
            bucket_container.add_bucket_entry(current_node as usize, (*index, distances[*index]));

            distances[*index] = INFINITY; //reset distance array for future iterations
        });

        terminals.clear();
        down_nodes[current_node as usize].clear();
    }
}

pub fn initialize_bwd_buckets_consecutive(
    bucket_container: &mut PreparationAccessBucketContainer<(usize, Weight)>, 
    bucket_entry_nodes: &[NodeId], 
    terminals: &mut Vec<usize>,
    distances: &mut Vec<Weight>,
    down_nodes: &mut Vec<Vec<(NodeId, Weight)>>, 
    queue: &mut MinBinaryHeap<RankDijkstraState>, 
    bwd_first_edge: &Vec<NodeId>, 
    bwd_arclist: &Vec<(NodeId, Weight)>, 
    ranks: &Vec<usize>
) {
    bucket_container.reset();

    for (target_index, target) in bucket_entry_nodes.iter().enumerate() {
        if !queue.contains_unique_index(*target as usize) {
            queue.insert(RankDijkstraState {rank: ranks[*target as usize], node_id: *target});
            bucket_container.add_preparation_entry(*target as usize, (target_index, 0));
        }
    }
    
    distances.resize(bucket_entry_nodes.len(), INFINITY);
    for index in 0..distances.len() {
        distances[index] = INFINITY;
    }

    while let Some(RankDijkstraState {rank: _, node_id: current_node}) = queue.pop() {
        let start_index = bwd_first_edge[current_node as usize] as usize;
        let end_index = bwd_first_edge[current_node as usize + 1] as usize;

        for (source_node, source_distance) in &bwd_arclist[start_index..end_index] {
            down_nodes[*source_node as usize].push((current_node, *source_distance));

            if !queue.contains_unique_index(*source_node as usize) {
                queue.insert(RankDijkstraState {rank: ranks[*source_node as usize], node_id: *source_node});
            }
        }

        for index in 0..distances.len() {
            distances[index] = INFINITY;
        }

        // iterate over all nodes that the current_node is connected to (edge: current_node -> down_node)
        for (down_node, down_distance) in &down_nodes[current_node as usize] {
            for (target_index, target_distance) in bucket_container.get_preparation_bucket(*down_node as usize) {
                let current_distance = &mut distances[*target_index];
                
                if *current_distance == INFINITY {
                    terminals.push(*target_index);
                }

                *current_distance = min(*current_distance, *target_distance + down_distance);
            }
        }

        for terminal in terminals.iter() {
            bucket_container.add_preparation_entry(current_node as usize, (*terminal, distances[*terminal]));

            distances[*terminal] = INFINITY; // reset distance array for future iterations
        }        

        terminals.clear();
        down_nodes[current_node as usize].clear();
    }

    bucket_container.prepare();
}

/// populates all the buckets of the associated nodes for a given backward search space
/// also applies pruning to reduce the number of bucket entries
pub fn populate_buckets_pruing_fwd<F>(
    bucket_container: &mut BucketContainer<(usize, Weight)>, 
    bucket_entry_nodes: &[NodeId], 
    terminals: &mut Vec<usize>,
    distances: &mut Vec<Weight>,
    down_nodes: &mut Vec<Vec<(NodeId, Weight)>>, 
    up_nodes: &mut Vec<Vec<(NodeId, Weight)>>,
    modified_up_nodes: &mut Vec<NodeId>, // keeps track of modified up nodes to clear them
    queue: &mut MinBinaryHeap<RankDijkstraState>, 
    bwd_first_edge: &Vec<NodeId>, 
    bwd_arclist: &Vec<(NodeId, Weight)>, 
    fwd_first_edge: &Vec<NodeId>,
    fwd_arclist: &Vec<(NodeId, Weight)>,
    ranks: &Vec<usize>,
    mut settled_node: F // currenlty settled node, terminals, distance to terminal
) where F: FnMut(NodeId, &Vec<usize>, &Vec<Weight>) { 
    bucket_container.reset();

    for (target_index, target) in bucket_entry_nodes.iter().enumerate() {
        if !queue.contains_unique_index(*target as usize) {
            queue.insert(RankDijkstraState {rank: ranks[*target as usize], node_id: *target});   
        }

        bucket_container.add_bucket_entry(*target as usize, (target_index, 0));
    }
    
    distances.resize(bucket_entry_nodes.len(), INFINITY);
    distances.iter_mut().for_each(|entry| *entry = INFINITY);
    
    unsafe {
        while let Some(RankDijkstraState {rank: _, node_id: current_node}) = queue.pop() {
            let start_index = *bwd_first_edge.get_unchecked(current_node as usize) as usize;
            let end_index = *bwd_first_edge.get_unchecked(current_node as usize + 1) as usize;

            for (source_node, source_distance) in bwd_arclist.get_unchecked(start_index..end_index) {
                down_nodes.get_unchecked_mut(*source_node as usize).push((current_node, *source_distance));

                if !queue.contains_unique_index(*source_node as usize) {
                    queue.insert(RankDijkstraState {rank: ranks[*source_node as usize], node_id: *source_node});
                }
            }

            let fwd_start_index = *fwd_first_edge.get_unchecked(current_node as usize) as usize;
            let fwd_end_index = *fwd_first_edge.get_unchecked(current_node as usize + 1) as usize;
            
            for (target_node, target_distance) in fwd_arclist.get_unchecked(fwd_start_index..fwd_end_index) {
                if up_nodes.get_unchecked(*target_node as usize).is_empty() {
                    modified_up_nodes.push(*target_node);
                }
                
                up_nodes.get_unchecked_mut(*target_node as usize).push((current_node, *target_distance));
            }
        
            // iterate over all nodes that the current_node is connected to via a downwards edge (edge: current_node -> down_node)
            for (down_node, down_distance) in down_nodes.get_unchecked(current_node as usize) {
                for (target_index, target_distance) in bucket_container.get_bucket(*down_node as usize) {
                    let current_distance = distances.get_unchecked_mut(*target_index);
                
                    if *current_distance == INFINITY {
                        terminals.push(*target_index);
                    }

                    *current_distance = cmp::min(*current_distance, target_distance + down_distance);
                }
            }

            // prune bucket entries of lower nodes by removing entries that can be reached faster via the current node
            for (up_node, up_distance) in up_nodes.get_unchecked(current_node as usize) {
                bucket_container.bucket_retain(*up_node as usize, |(target_index, target_distance)| *target_distance < *up_distance + distances.get_unchecked(*target_index));
            }

            settled_node(current_node, terminals, distances);

            terminals.iter().for_each(|index| {
                bucket_container.add_bucket_entry(current_node as usize, (*index, distances[*index]));

                *distances.get_unchecked_mut(*index) = INFINITY; // reset distance array for future iterations
            });

            terminals.clear();
            down_nodes.get_unchecked_mut(current_node as usize).clear();
        }

        // reset up nodes for future iterations, we cannot clear them during the iteration as we have initialized it for more vertices than we settle
        for up_node_index in modified_up_nodes.iter() {
            up_nodes.get_unchecked_mut(*up_node_index as usize).clear();
        }
        modified_up_nodes.clear();
    } 
}