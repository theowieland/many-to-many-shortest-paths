use core_simd::SimdPartialOrd;
use core_simd::u32x16;
use crate::types::*;
use crate::data_structures::*;
use crate::many_to_many::many_to_many_algorithm::ManyToManyAlgorithm;
use crate::many_to_many::many_to_many_utils::*;
use crate::many_to_many::rphast::target_selection;
use crate::utils::binary_heap::MinBinaryHeap;

type SIMDType = u32x16;
pub const BATCH_SIZE: usize = 16;

pub struct SIMDRPHAST {
    ranks: Vec<usize>,
    dfs_discoverd: ValidFlags<bool>,
    fwd_first_edge: Vec<EdgeId>,
    fwd_arclist: Vec<(NodeId, Weight)>,
    bwd_first_edge: Vec<EdgeId>,
    bwd_arclist: Vec<(NodeId, Weight)>,
    bwd_visited_nodes: Vec<NodeId>,
    bwd_restricted_nodes_translation: ValidFlags<NodeId>,
    bwd_restricted_first_edge: Vec<EdgeId>,
    bwd_restricted_arclist: Vec<(NodeId, Weight)>,
    bwd_restricted_distances: Vec<SIMDType>,
    fwd_distances: Vec<SIMDType>,
    fwd_queue: MinBinaryHeap<RankDijkstraState>,
    distance_table: Matrix<Weight>
}

impl SIMDRPHAST {

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

        SIMDRPHAST {
            ranks: ranks.to_vec(),
            dfs_discoverd: ValidFlags::new(num_vertices, false),
            fwd_first_edge: fwd_first_edge.to_owned(),
            fwd_arclist: convert_to_arclist(fwd_head, fwd_weight),
            bwd_first_edge: bwd_first_edge.to_owned(),
            bwd_arclist: convert_to_arclist(bwd_head, bwd_weight),
            bwd_visited_nodes: Vec::new(),
            bwd_restricted_nodes_translation: ValidFlags::new(num_vertices, INFINITY), //translate from node id to local id
            bwd_restricted_first_edge: Vec::new(),
            bwd_restricted_arclist: Vec::new(),
            bwd_restricted_distances: Vec::new(),
            fwd_distances: vec![SIMDType::splat(INFINITY); num_vertices],
            fwd_queue: MinBinaryHeap::new(num_vertices),
            distance_table: Matrix::empty()
        }
    }

    /// build the restricted target search space
    pub fn select_targets(&mut self, targets: &[NodeId]) {
        target_selection(&self.bwd_first_edge, &self.bwd_arclist, &mut self.bwd_restricted_nodes_translation, &mut self.dfs_discoverd, targets, &mut self.bwd_visited_nodes, &mut self.bwd_restricted_first_edge, &mut self.bwd_restricted_arclist);

        // resize bwd_restricted_distances to fit current targets
        self.bwd_restricted_distances.resize(self.bwd_visited_nodes.len(), SIMDType::splat(INFINITY));
    }

    pub fn query(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        self.distance_table.resize(sources.len(), targets.len(), INFINITY);

        batch_nodes(BATCH_SIZE, sources, |batch_index, batched_nodes| {

            // clear the bwd restricted distances from previous queries
            self.bwd_restricted_distances.iter_mut().for_each(|entry| *entry = SIMDType::splat(INFINITY));

            // start the fwd search by adding the batched sources with their initial 0 distance
            for (batch_source_index, source) in batched_nodes.iter().enumerate() {
                self.fwd_queue.insert(RankDijkstraState {rank: self.ranks[*source as usize], node_id: *source});
                self.fwd_distances[*source as usize][batch_source_index] = 0;
            } 

            // fwd search
            while let Some(RankDijkstraState {rank: _, node_id: current_node}) = self.fwd_queue.pop() { // pop current node with the lowest rank
                unsafe {
                    let start = *self.fwd_first_edge.get_unchecked(current_node as usize) as usize;
                    let end = *self.fwd_first_edge.get_unchecked(current_node as usize + 1) as usize;

                    let current_distances = *self.fwd_distances.get_unchecked(current_node as usize);

                    for (target_node, edge_distance) in self.fwd_arclist.get_unchecked(start..end) {
                        let old_distances = *self.fwd_distances.get_unchecked_mut(*target_node as usize);
                        
                        // vec that stores edge_offset for all lanes
                        let edge_offset = SIMDType::splat(*edge_distance);

                        // true if distance can be reduced
                        let mask = old_distances.simd_gt(current_distances + edge_offset);    
                    
                        // asign reduced distances
                        self.fwd_distances[*target_node as usize] = mask.select(current_distances + edge_offset, old_distances);

                        if !self.fwd_queue.contains_unique_index(*target_node as usize) {
                            self.fwd_queue.insert(RankDijkstraState {rank: *self.ranks.get_unchecked(*target_node as usize), node_id: *target_node});
                        }   
                    }

                    // copy distances to the bwd_restricted_distance array
                    if self.bwd_restricted_nodes_translation.is_valid(current_node as usize) { // the currently visisted node is also in the target set
                        let bwd_restricted_node_index = self.bwd_restricted_nodes_translation.get_unchecked(current_node as usize) as usize;

                        *self.bwd_restricted_distances.get_unchecked_mut(bwd_restricted_node_index) = current_distances;
                    }

                    // clear old distances that are never accessed again (rank order)
                    *self.fwd_distances.get_unchecked_mut(current_node as usize) = SIMDType::splat(INFINITY);
                }
            }

            unsafe {
                // iterate over the bwd nodes in order of discovery. top to bottom
                for node_index in 0..self.bwd_visited_nodes.len() {
                    let start = *self.bwd_restricted_first_edge.get_unchecked(node_index) as usize;
                    let end = *self.bwd_restricted_first_edge.get_unchecked(node_index + 1) as usize;

                    // iterate over all incoming arcs to the current_node (source_node is a local id given in the range of the bwd_restricted search space)
                    for (arc_source_node, edge_distance) in self.bwd_restricted_arclist.get_unchecked(start..end) {
                        let current_distances = *self.bwd_restricted_distances.get_unchecked(node_index);
                        
                        let old_distances = *self.bwd_restricted_distances.get_unchecked(*arc_source_node as usize);
                
                        // vec that stores edge_offset for all lanes
                        let edge_offset = SIMDType::splat(*edge_distance);
                    
                        // true if distance can be reduced
                        let mask = current_distances.simd_gt(old_distances + edge_offset);

                        // asign reduced distances
                        self.bwd_restricted_distances[node_index] = mask.select(old_distances + edge_offset, current_distances);
                    }
                }
            }
            
            //update the distance table
            unsafe {
                for (batch_source_index, _batch_source) in batched_nodes.iter().enumerate() {
                    for (target_index, target) in targets.iter().enumerate() {
                        self.distance_table.set_unsafe(batch_source_index + (batch_index * BATCH_SIZE), target_index, self.bwd_restricted_distances[(self.bwd_restricted_nodes_translation[*target as usize] as usize)][batch_source_index]);
                    }
                }
            }
        });
    }
}

impl ManyToManyAlgorithm for SIMDRPHAST {

    fn initialize(&mut self, _sources: &[NodeId], _targets: &[NodeId]) {
        // no initialization required
    }

    fn select(&mut self, _sources: &[NodeId], targets: &[NodeId]) {
        self.select_targets(targets);
    }

    fn query(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        self.query(sources, targets);
    }

    /// return the previously calculated distance table
    fn get_distance_table(&mut self) -> &Matrix<Weight> {
        &self.distance_table
    }
}