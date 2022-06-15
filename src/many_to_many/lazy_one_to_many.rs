use crate::types::*;
use crate::graph_algorithms::DijkstraState;
use crate::data_structures::{ResettableArray, ArrayStructure, ValidFlags, Matrix};
use crate::utils::binary_heap::MinBinaryHeap;
use std::cmp;
use crate::many_to_many::many_to_many_algorithm::ManyToManyAlgorithm;

use super::many_to_many_utils::convert_to_arclist;

/// lazy rphast for the one to many problem
/// selects a single source vertex and caluclate distances from the target vertices
/// to that selected source
pub struct LazyRPHASTOneToMany {
    fwd_first_edge: Vec<EdgeId>, 
    fwd_arclist: Vec<(NodeId, Weight)>,
    bwd_first_edge: Vec<EdgeId>,
    bwd_arclist: Vec<(NodeId, Weight)>,
    fwd_distances: ResettableArray<Weight>,
    bwd_distances: ValidFlags<Weight>,
    fwd_queue: MinBinaryHeap<DijkstraState>,
    distance_table: Matrix<Weight>
}

impl LazyRPHASTOneToMany { 

    pub fn new(
        fwd_first_edge: &Vec<EdgeId>, 
        fwd_head: &Vec<NodeId>, 
        fwd_weight: &Vec<Weight>, 
        bwd_first_edge: &Vec<EdgeId>, 
        bwd_head: &Vec<NodeId>, 
        bwd_weight: &Vec<Weight>
    ) -> Self {
        let num_vertices = fwd_first_edge.len() - 1;

        LazyRPHASTOneToMany {
            fwd_first_edge: fwd_first_edge.to_owned(),
            fwd_arclist: convert_to_arclist(fwd_head, fwd_weight),
            bwd_first_edge: bwd_first_edge.to_owned(),
            bwd_arclist: convert_to_arclist(bwd_head, bwd_weight),
            fwd_distances: ResettableArray::<Weight>::new(num_vertices, INFINITY),
            bwd_distances: ValidFlags::<Weight>::new(num_vertices, INFINITY),
            fwd_queue: MinBinaryHeap::new(num_vertices),
            distance_table: Matrix::<Weight>::empty()
        }
    }


    /// select a single source vertex for one-to-many queries
    pub fn select_source(&mut self, source: NodeId) {
        self.fwd_distances.reset();
        self.bwd_distances.reset();
        self.fwd_queue.clear();

        self.fwd_distances.set(source as usize, 0);
        self.fwd_queue.insert(DijkstraState {distance: 0, node_id: source});

        while let Some(DijkstraState { distance: current_distance, node_id: current_node }) = self.fwd_queue.pop() {
            unsafe {
                let start = *self.fwd_first_edge.get_unchecked(current_node as usize) as usize;
                let end = *self.fwd_first_edge.get_unchecked(current_node as usize + 1) as usize;

                for (target_node, edge_distance) in self.fwd_arclist.get_unchecked(start..end) {
                    let old_distance = self.fwd_distances.get_unchecked(*target_node as usize);

                    if edge_distance + current_distance < *old_distance {
                        self.fwd_queue.insert_or_decrease(DijkstraState {distance: edge_distance + current_distance, node_id: *target_node});
                        self.fwd_distances.set(*target_node as usize, edge_distance + current_distance);
                    }
                }
            }
        }
    }

    pub fn compute_and_memoize_dist(&mut self, node: NodeId) -> Weight {
        LazyRPHASTOneToMany::rust_mut_borrow_bypass(node, &self.bwd_first_edge, &self.bwd_arclist, &mut self.fwd_distances, &mut self.bwd_distances)
    }

    fn rust_mut_borrow_bypass(node: NodeId, bwd_first_edge: &Vec<EdgeId>, bwd_arclist: &Vec<(NodeId, Weight)>, fwd_distances: &mut ResettableArray<Weight>, bwd_distances: &mut ValidFlags<Weight>) -> Weight {
        if !bwd_distances.is_valid(node as usize) {
            bwd_distances.set(node as usize, fwd_distances[node as usize]);
            
            unsafe {
                let start = *bwd_first_edge.get_unchecked(node as usize) as usize;
                let end = *bwd_first_edge.get_unchecked(node as usize + 1) as usize;
                
                for (source_node, edge_distance) in bwd_arclist.get_unchecked(start..end) {
                    let recursive_distance = LazyRPHASTOneToMany::rust_mut_borrow_bypass(*source_node, bwd_first_edge, bwd_arclist, fwd_distances, bwd_distances);
    
                    bwd_distances.set(node as usize, cmp::min(
                        bwd_distances.get_unchecked(node as usize),
                        edge_distance + recursive_distance
                    ));
                }
            }
        }

        bwd_distances[node as usize]
    }
}

impl ManyToManyAlgorithm for LazyRPHASTOneToMany {

    fn initialize(&mut self, _sources: &[NodeId], _targets: &[NodeId]) {
        // no initialization required
    }

    fn select(&mut self, _sources: &[NodeId], _targets: &[NodeId]) {
        // no selection required
    }

    fn query(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        self.distance_table.resize(sources.len(), targets.len(), INFINITY);

        for (source_index, source) in sources.iter().enumerate() {
            self.select_source(*source);

            for (target_index, target) in targets.iter().enumerate() {
                let distance = self.compute_and_memoize_dist(*target);
                self.distance_table.set(source_index, target_index, distance);
            }
        }
    }

    fn get_distance_table(&mut self) -> &Matrix<Weight> {
        &self.distance_table
    }
}