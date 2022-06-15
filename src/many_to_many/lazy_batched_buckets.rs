use crate::bucket_containers::BucketContainer;
use crate::types::*;
use crate::data_structures::{Matrix, LazyBatchedValidFlags};
use crate::utils::binary_heap::MinBinaryHeap;
use std::cmp::min;
use crate::many_to_many::many_to_many_algorithm::ManyToManyAlgorithm;
use crate::many_to_many::many_to_many_utils::*;

use super::advanced_bucket::{populate_buckets_pruing_fwd};

pub struct LazyBatchedBuckets {
    batch_size: usize,
    ranks: Vec<usize>,
    fwd_first_edge: Vec<EdgeId>,
    fwd_arclist: Vec<(NodeId, Weight)>,
    bwd_first_edge: Vec<EdgeId>,
    bwd_arclist: Vec<(NodeId, Weight)>,
    bwd_buckets: BucketContainer<(usize, Weight)>,
    bwd_distances: Vec<Weight>,
    bwd_terminals: Vec<usize>,
    bwd_down_nodes: Vec<Vec<(NodeId, Weight)>>,
    bwd_up_nodes: Vec<Vec<(NodeId, Weight)>>,
    bwd_modified_up_nodes: Vec<NodeId>,
    fwd_distances: LazyBatchedValidFlags,
    bwd_queue: MinBinaryHeap<RankDijkstraState>,
    distance_table: Matrix<Weight>
}

impl LazyBatchedBuckets { 

    pub fn new(
        batch_size: usize, 
        fwd_first_edge: &Vec<EdgeId>, 
        fwd_head: &Vec<NodeId>, 
        fwd_weight: &Vec<Weight>, 
        bwd_first_edge: &Vec<EdgeId>, 
        bwd_head: &Vec<NodeId>, 
        bwd_weight: &Vec<Weight>, 
        ranks: &Vec<usize>
    ) -> Self {
        let num_vertices = fwd_first_edge.len() - 1;

        LazyBatchedBuckets {
            batch_size,
            ranks: ranks.to_vec(),
            fwd_first_edge: fwd_first_edge.to_owned(),
            fwd_arclist: convert_to_arclist(fwd_head, fwd_weight),
            bwd_first_edge: bwd_first_edge.to_owned(),
            bwd_arclist: convert_to_arclist(bwd_head, bwd_weight),
            bwd_buckets: BucketContainer::new(num_vertices),
            bwd_distances: Vec::new(),
            bwd_terminals: Vec::new(),
            bwd_down_nodes: vec![Vec::new(); num_vertices],
            bwd_up_nodes: vec![Vec::new(); num_vertices],
            bwd_modified_up_nodes: Vec::new(),
            fwd_distances: LazyBatchedValidFlags::new(num_vertices, batch_size),
            bwd_queue: MinBinaryHeap::new(num_vertices),
            distance_table: Matrix::<Weight>::empty(),
        }
    }

    /// initializes the distances of all nodes discovered by a backward search starting from the given target nodes
    pub fn select_targets(&mut self, targets: &[NodeId]) {
        populate_buckets_pruing_fwd(
            &mut self.bwd_buckets, 
            targets, 
            &mut self.bwd_terminals, 
            &mut self.bwd_distances, 
            &mut self.bwd_down_nodes, 
            &mut self.bwd_up_nodes,
            &mut self.bwd_modified_up_nodes,
            &mut self.bwd_queue, 
            &self.bwd_first_edge, 
            &self.bwd_arclist, 
            &self.fwd_first_edge,
            &self.fwd_arclist,
            &self.ranks,
            |_settled_node, _terminals, _distances| {}
        );
    }

    pub fn compute_and_memoize_dist(&mut self, node: NodeId) {
        rust_mut_borrow_bypass(self.batch_size, node, &self.fwd_first_edge, &self.fwd_arclist, &mut self.fwd_distances, &mut self.bwd_buckets);
    }
}

fn rust_mut_borrow_bypass(
    batch_size: usize, 
    node: NodeId, 
    fwd_first_edge: &Vec<EdgeId>, 
    fwd_arclist: &Vec<(NodeId, Weight)>, 
    fwd_distances: &mut LazyBatchedValidFlags,
    bwd_buckets: &mut BucketContainer<(usize, Weight)>    
) {
    if !fwd_distances.is_valid(node as usize) {
        fwd_distances.initialize_batch(node as usize);

        for (target_index, target_distance) in bwd_buckets.get_bucket(node as usize) {
            fwd_distances.set(node as usize, *target_index, *target_distance);
        }
        
        let start = fwd_first_edge[node as usize] as usize;
        let end = fwd_first_edge[node as usize + 1] as usize;
        
        for (target_node, edge_distance) in &fwd_arclist[start..end] {
            rust_mut_borrow_bypass(batch_size, *target_node, fwd_first_edge, fwd_arclist, fwd_distances, bwd_buckets);
            fwd_distances.reduce_batch(node as usize, *edge_distance, *target_node as usize);
        }
    }
}

impl ManyToManyAlgorithm for LazyBatchedBuckets {

    fn initialize(&mut self, _sources: &[NodeId], _targets: &[NodeId]) {
        // no initialization required
    }

    fn select(&mut self, _sources: &[NodeId], _targets: &[NodeId]) {
        // no selection required
    }

    fn query(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        self.distance_table.resize(sources.len(), targets.len(), INFINITY);

        let num_batches = (targets.len() / self.batch_size) + if targets.len() % self.batch_size != 0 {1} else {0};

        for batch in 0..num_batches {
            let batch_first_target_index = batch * self.batch_size; //index of first source node in the current batch
            let batch_last_target_index = min((batch + 1) * self.batch_size, targets.len()); //index of the last source node in the current batch

            self.select_targets(&targets[batch_first_target_index..batch_last_target_index]);
        
            self.fwd_distances.reset();

            for (source_index, source) in sources.iter().enumerate() {
                self.compute_and_memoize_dist(*source);
                let distances = self.fwd_distances.get_row(*source as usize).as_ref().unwrap();

                for (target_index_in_batch, _target) in targets[batch_first_target_index..batch_last_target_index].iter().enumerate() {
                    self.distance_table.set(source_index, batch * self.batch_size + target_index_in_batch, distances[target_index_in_batch]);
                }
            }
        }
    }

    fn get_distance_table(&mut self) -> &Matrix<Weight> {
        &self.distance_table
    }
}