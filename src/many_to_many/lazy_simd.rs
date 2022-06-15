use core_simd::{u32x16, SimdPartialOrd};

use crate::bucket_containers::BucketContainer;
use crate::types::*;
use crate::data_structures::{Matrix, ValidFlags, ArrayStructure};
use crate::many_to_many::many_to_many_algorithm::ManyToManyAlgorithm;
use crate::many_to_many::many_to_many_utils::*;
use crate::utils::binary_heap::MinBinaryHeap;

use super::advanced_bucket::{populate_buckets_pruing_fwd};

type SIMDType = u32x16;
pub const BATCH_SIZE: usize = 16;

pub struct LazySIMD {
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
    fwd_distances: ValidFlags<SIMDType>,
    bwd_queue: MinBinaryHeap<RankDijkstraState>,
    distance_table: Matrix<Weight>
}

impl LazySIMD { 

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

        LazySIMD {
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
            fwd_distances: ValidFlags::new(num_vertices, SIMDType::splat(INFINITY)),
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
        rust_mut_borrow_bypass(node, &self.fwd_first_edge, &self.fwd_arclist, &mut self.fwd_distances, &mut self.bwd_buckets);
    }
}

fn rust_mut_borrow_bypass(
    node: NodeId, 
    fwd_first_edge: &Vec<EdgeId>, 
    fwd_arclist: &Vec<(NodeId, Weight)>, 
    fwd_distances: &mut ValidFlags<SIMDType>,
    bwd_buckets: &mut BucketContainer<(usize, Weight)>    
) {
    if !fwd_distances.is_valid(node as usize) {
        fwd_distances.set(node as usize, SIMDType::splat(INFINITY));

        for (target_index, target_distance) in bwd_buckets.get_bucket(node as usize) {    
            let mut current_values = fwd_distances[node as usize];
            current_values[*target_index] = *target_distance;
            
            fwd_distances.set(node as usize, current_values);
        }
        
        let start = fwd_first_edge[node as usize] as usize;
        let end = fwd_first_edge[node as usize + 1] as usize;
        
        for (target_node, edge_distance) in &fwd_arclist[start..end] {
            rust_mut_borrow_bypass(*target_node, fwd_first_edge, fwd_arclist, fwd_distances, bwd_buckets);
            
            let edge_distance_splat = SIMDType::splat(*edge_distance);
            let target_distances = fwd_distances[*target_node as usize];

            let mask = fwd_distances[node as usize].simd_gt(target_distances + edge_distance_splat);
            
            fwd_distances.set(node as usize, mask.select(target_distances + edge_distance_splat, fwd_distances[node as usize]));
        }
    }
}

impl ManyToManyAlgorithm for LazySIMD {

    fn initialize(&mut self, _sources: &[NodeId], _targets: &[NodeId]) {
        // no initialization required
    }

    fn select(&mut self, _sources: &[NodeId], _targets: &[NodeId]) {
        // no selection required
    }

    fn query(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        self.distance_table.resize(sources.len(), targets.len(), INFINITY);

        batch_nodes(BATCH_SIZE, targets, |batch_index, batched_nodes| {
            self.select_targets(batched_nodes);
            self.fwd_distances.reset();

            for (source_index, source) in sources.iter().enumerate() {
                self.compute_and_memoize_dist(*source);
                let distances = self.fwd_distances[*source as usize];

                for (target_index_in_batch, _target) in batched_nodes.iter().enumerate() {
                    self.distance_table.set(source_index, (batch_index * BATCH_SIZE) + target_index_in_batch, distances[target_index_in_batch]);
                }
            }
        });
    }

    fn get_distance_table(&mut self) -> &Matrix<Weight> {
        &self.distance_table
    }
}