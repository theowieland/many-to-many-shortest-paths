use crate::types::*;
use crate::many_to_many::many_to_many_algorithm::ManyToManyAlgorithm;
use crate::many_to_many::many_to_many_utils::*;
use crate::utils::binary_heap::MinBinaryHeap;
use crate::utils::bucket_containers::BucketContainer;
use crate::utils::data_structures::Matrix;
use super::advanced_bucket::populate_buckets_pruing_fwd;
use super::bucket::query_no_stopping;

pub struct AdvancedBucketBatchedManyToMany {
    batch_size: usize,
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
    distance_table: Matrix<Weight>
}

impl AdvancedBucketBatchedManyToMany {

    pub fn new(
        batch_size: usize, 
        fwd_first_edge: &EdgeIds, 
        fwd_head: &NodeIds, 
        fwd_weight: &Weights, 
        bwd_first_edge: &EdgeIds, 
        bwd_head: &NodeIds, 
        bwd_weight: &Weights, 
        ranks: &Ranks
    ) -> Self {
        let num_vertices = fwd_first_edge.len() - 1;

        AdvancedBucketBatchedManyToMany {
            batch_size,
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
}

impl ManyToManyAlgorithm for AdvancedBucketBatchedManyToMany {
    
    fn initialize(&mut self, _sources: &[NodeId], _targets: &[NodeId]) {
        // no initialization required
    }

    fn select(&mut self, _sources: &[NodeId], _targets: &[NodeId]) {
        // no selection required
    }

    fn query(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        self.distance_table.resize(sources.len(), targets.len(), INFINITY);
        
        batch_nodes(self.batch_size, targets, |batch_index, batched_nodes| {
            self.select_targets(batched_nodes);

            for (source_index, source) in sources.iter().enumerate() {
                self.fwd_distances[*source as usize] = 0;
                self.fwd_queue.insert(RankDijkstraState {rank: self.ranks[*source as usize], node_id: *source});
        
                while let Some(RankDijkstraState {rank: _rank, node_id: current_node}) = self.fwd_queue.pop() {
                    unsafe {
                        let current_distance = *self.fwd_distances.get_unchecked(current_node as usize);
                        let current_bucket_entries = self.bwd_buckets.get_bucket(current_node as usize);
                        self.distance_table.reduce_buckets_batch_offset(source_index, current_bucket_entries, current_distance, batch_index * self.batch_size);
                        
                        let start_index = *self.fwd_first_edge.get_unchecked(current_node as usize) as usize;
                        let end_index = *self.fwd_first_edge.get_unchecked(current_node as usize + 1) as usize;
                    
                        for (target_node, edge_distance) in self.fwd_arclist.get_unchecked(start_index..end_index) {
                            let old_distance = self.fwd_distances.get_unchecked_mut(*target_node as usize);
        
                            if edge_distance + current_distance < *old_distance {
                                *old_distance = edge_distance + current_distance;
                            } 
        
                            if !self.fwd_queue.contains_unique_index(*target_node as usize) {
                                self.fwd_queue.insert(RankDijkstraState {rank: self.ranks[*target_node as usize], node_id: *target_node});
                            } 
                        }
        
                        *self.fwd_distances.get_unchecked_mut(current_node as usize) = INFINITY;
                    }
                }
            }                
        });
    }

    fn get_distance_table(&mut self) -> &Matrix<Weight> {
        &self.distance_table
    }
}