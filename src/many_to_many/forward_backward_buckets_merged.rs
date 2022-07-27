use crate::bucket_containers::BucketContainer;
use crate::types::*;
use crate::data_structures::{Matrix, EmptyValidFlags};
use crate::many_to_many::many_to_many_algorithm::ManyToManyAlgorithm;
use crate::many_to_many::many_to_many_utils::*;
use crate::utils::binary_heap::MinBinaryHeap;

use super::advanced_bucket::populate_buckets_pruing_fwd;

pub struct ForwardBackwardBucketsMergedManyToMany {
    ranks: Vec<usize>,
    fwd_first_edge: Vec<EdgeId>,
    fwd_arclist: Vec<(NodeId, Weight)>,
    bwd_first_edge: Vec<EdgeId>,
    bwd_arclist: Vec<(NodeId, Weight)>,
    distances: Vec<Weight>,
    queue: MinBinaryHeap<RankDijkstraState>,
    terminals: Vec<usize>,
    fwd_buckets: BucketContainer<(usize, Weight)>,
    bwd_buckets: BucketContainer<(usize, Weight)>,
    down_nodes: Vec<Vec<(NodeId, Weight)>>,
    up_nodes: Vec<Vec<(NodeId, Weight)>>,
    modified_up_nodes: Vec<NodeId>,
    distance_table: Matrix<Weight>,

    mutual_buckets: EmptyValidFlags
}

impl ForwardBackwardBucketsMergedManyToMany {

    pub fn new(
        fwd_first_edge: &EdgeIds, 
        fwd_head: &NodeIds, 
        fwd_weight: &Weights, 
        bwd_first_edge: &EdgeIds, 
        bwd_head: &NodeIds, 
        bwd_weight: &Weights, 
        ranks: &Ranks
    ) -> Self {
        let num_vertices = fwd_first_edge.len() - 1;

        ForwardBackwardBucketsMergedManyToMany {
            ranks: ranks.to_vec(),
            fwd_first_edge: fwd_first_edge.to_owned(),
            fwd_arclist: convert_to_arclist(fwd_head, fwd_weight),
            bwd_first_edge: bwd_first_edge.to_owned(),
            bwd_arclist: convert_to_arclist(bwd_head, bwd_weight),
            distances: Vec::new(),
            queue: MinBinaryHeap::new(num_vertices),
            terminals: Vec::new(),
            fwd_buckets: BucketContainer::new(num_vertices),
            bwd_buckets: BucketContainer::new(num_vertices),
            down_nodes: vec![Vec::new(); num_vertices],
            up_nodes: vec![Vec::new(); num_vertices],
            modified_up_nodes: Vec::new(),
            distance_table: Matrix::empty(),
            mutual_buckets: EmptyValidFlags::new(num_vertices)
        }
    }

    pub fn populate_buckets(&mut self, sources: &[NodeId], targets: &[NodeId]) {

        // initialize the forward buckets
        populate_buckets_pruing_fwd(
            &mut self.fwd_buckets,
            sources,
            &mut self.terminals,
            &mut self.distances,
            &mut self.down_nodes,
            &mut self.up_nodes,
            &mut self.modified_up_nodes,
            &mut self.queue,
            &self.fwd_first_edge,
            &self.fwd_arclist,
            &self.bwd_first_edge,
            &self.bwd_arclist,
            &self.ranks,
            |_settled_node, _terminals, _distances| {}
        );

        // initialize the backward buckets
        populate_buckets_pruing_fwd(
            &mut self.bwd_buckets,
            targets,
            &mut self.terminals,
            &mut self.distances,
            &mut self.down_nodes,
            &mut self.up_nodes,
            &mut self.modified_up_nodes,
            &mut self.queue,
            &self.bwd_first_edge,
            &self.bwd_arclist,
            &self.fwd_first_edge,
            &self.fwd_arclist,
            &self.ranks,
            |_settled_node, _terminals, _distances| {}
        );
    }

    pub fn calculate_mutual_buckets(&mut self) -> Vec<NodeId> {

        // discover buckets visited by both the forward and backward search
        self.mutual_buckets.reset();
        let mut visited: Vec<NodeId> = Vec::new();
        for forward_visited in self.fwd_buckets.get_modified_buckets() {
            self.mutual_buckets.set_valid(*forward_visited as usize);
        }

        for backward_visited in self.bwd_buckets.get_modified_buckets() {
            if self.mutual_buckets.is_valid(*backward_visited as usize) { // bucket was visited by the fwd search
                visited.push(*backward_visited as NodeId);
            }
        }

        visited
    }

    pub fn merge_buckets(&mut self, mutual_buckets: &Vec<NodeId>) {
        for visited_bucket in mutual_buckets {
            for (source_index, source_distance) in self.fwd_buckets.get_bucket(*visited_bucket as usize) {
                for (target_index, target_distance) in self.bwd_buckets.get_bucket(*visited_bucket as usize) {
                    self.distance_table.reduce_value(*source_index, *target_index, *source_distance + *target_distance);
                }
            }
        }
    }

    pub fn get_total_bucket_size(&mut self, mutual_buckets: &Vec<NodeId>) -> (usize, usize) {
        let mut fwd_bucket_total_size = 0;
        let mut bwd_bucket_total_size = 0;

        for mutual_bucket in mutual_buckets {
            fwd_bucket_total_size += self.fwd_buckets.get_bucket(*mutual_bucket as usize).len();
            bwd_bucket_total_size += self.bwd_buckets.get_bucket(*mutual_bucket as usize).len();
        }

        (fwd_bucket_total_size, bwd_bucket_total_size)
    }
}

impl ManyToManyAlgorithm for ForwardBackwardBucketsMergedManyToMany {
    
    fn initialize(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        // no initialization required
        self.distance_table.resize(sources.len(), targets.len(), INFINITY);
    }

    fn select(&mut self, _sources: &[NodeId], _targets: &[NodeId]) {
        // no selection required
    }

    fn query(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        self.populate_buckets(sources, targets);
        let mutual_buckets = self.calculate_mutual_buckets();
        self.merge_buckets(&mutual_buckets);
    }

    fn get_distance_table(&mut self) -> &Matrix<Weight> {
        &self.distance_table
    }
}