use crate::graph_algorithms::DijkstraState;
use crate::types::*;
use crate::utils::binary_heap::MinBinaryHeap;
use crate::utils::data_structures::{ResettableArray, ValidFlags, Matrix, ArrayStructure};
use std::cmp;
use crate::many_to_many::many_to_many_algorithm::ManyToManyAlgorithm;

use super::many_to_many_utils::convert_to_arclist;

pub struct LazyRPHAST {
    fwd_first_edge: Vec<EdgeId>,
    fwd_arclist: Vec<(NodeId, Weight)>,
    bwd_first_edge: Vec<EdgeId>,
    bwd_arclist: Vec<(NodeId, Weight)>,
    bwd_distances: ResettableArray<Weight>,
    fwd_distances: ValidFlags<Weight>,
    bwd_queue: MinBinaryHeap<DijkstraState>,
    distance_table: Matrix<Weight>
}

impl LazyRPHAST { 

    pub fn new(
        fwd_first_edge: &EdgeIds, 
        fwd_head: &NodeIds, 
        fwd_weight: &Weights, 
        bwd_first_edge: &EdgeIds, 
        bwd_head: &NodeIds, 
        bwd_weight: &Weights
    ) -> Self {
        let num_vertices = fwd_first_edge.len() - 1;

        LazyRPHAST {
            fwd_first_edge: fwd_first_edge.to_owned(),
            fwd_arclist: convert_to_arclist(fwd_head, fwd_weight),
            bwd_first_edge: bwd_first_edge.to_owned(),
            bwd_arclist: convert_to_arclist(bwd_head, bwd_weight),
            bwd_distances: ResettableArray::<Weight>::new(num_vertices, INFINITY), // distance array used to calculate the backward tentative distances
            fwd_distances: ValidFlags::<Weight>::new(num_vertices, INFINITY), 
            bwd_queue: MinBinaryHeap::new(num_vertices),
            distance_table: Matrix::<Weight>::empty()
        }
    }

    /// selects a single target for future many to one calculations
    pub fn select_target(&mut self, target: NodeId) {
        self.bwd_distances.reset();
        self.fwd_distances.reset();
        self.bwd_queue.clear();

        self.bwd_distances.set(target as usize, 0);
        self.bwd_queue.insert(DijkstraState {distance: 0, node_id: target});

        while let Some(DijkstraState {distance: current_distance, node_id: current_node}) = self.bwd_queue.pop() {
            unsafe {
                let start = *self.bwd_first_edge.get_unchecked(current_node as usize) as usize;
                let end = *self.bwd_first_edge.get_unchecked(current_node as usize + 1) as usize;
    
                for (source_node, edge_distance) in self.bwd_arclist.get_unchecked(start..end) {
                    let old_distance = self.bwd_distances.get_unchecked(*source_node as usize);
    
                    if edge_distance + current_distance < *old_distance {
                        self.bwd_queue.insert_or_decrease(DijkstraState {distance: edge_distance + current_distance, node_id: *source_node});
                        self.bwd_distances.set(*source_node as usize, edge_distance + current_distance);
                    }
                }
            }
        }
    }

    pub fn compute_and_memoize_dist(&mut self, node: NodeId) -> Weight {
        LazyRPHAST::rust_mut_borrow_bypass(node, &self.fwd_first_edge, &self.fwd_arclist, &mut self.fwd_distances, &mut self.bwd_distances)
    }

    fn rust_mut_borrow_bypass(node: NodeId, fwd_first_edge: &Vec<EdgeId>, fwd_arclist: &Vec<(NodeId, Weight)>, fwd_distances: &mut ValidFlags<Weight>, bwd_distances: &mut ResettableArray<Weight>) -> Weight {
        if !fwd_distances.is_valid(node as usize) {
            fwd_distances.set(node as usize, bwd_distances[node as usize]);
            
            unsafe {
                let start = *fwd_first_edge.get_unchecked(node as usize) as usize;
                let end = *fwd_first_edge.get_unchecked(node as usize + 1) as usize;
                
                for (target_node, edge_distance) in fwd_arclist.get_unchecked(start..end) {
                    let recursive_distance = LazyRPHAST::rust_mut_borrow_bypass(*target_node, fwd_first_edge, fwd_arclist, fwd_distances, bwd_distances);
    
                    fwd_distances.set(node as usize, cmp::min(
                        fwd_distances.get_unchecked(node as usize),
                        edge_distance + recursive_distance
                    ));
                }
            }
        }

        fwd_distances[node as usize]
    }
}

impl ManyToManyAlgorithm for LazyRPHAST {

    fn initialize(&mut self, _sources: &[NodeId], _targets: &[NodeId]) {
        // no initialization required
    }

    fn select(&mut self, _sources: &[NodeId], _targets: &[NodeId]) {
        // no selection required
    }

    fn query(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        self.distance_table.resize(sources.len(), targets.len(), INFINITY);

        for (target_index, target) in targets.iter().enumerate() {
            self.select_target(*target);

            for (source_index, source) in sources.iter().enumerate() {
                let distance = self.compute_and_memoize_dist(*source);
                self.distance_table.set(source_index, target_index, distance);
            }
        }
    }

    fn get_distance_table(&mut self) -> &Matrix<Weight> {
        &self.distance_table
    }
}