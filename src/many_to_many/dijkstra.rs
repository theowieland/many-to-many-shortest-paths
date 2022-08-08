use crate::types::*;
use crate::graph_algorithms::DijkstraState;
use crate::many_to_many::many_to_many_algorithm::ManyToManyAlgorithm;
use crate::utils::binary_heap::MinBinaryHeap;
use crate::utils::data_structures::{ValidFlags, ResettableArray, Matrix, ArrayStructure};

use super::many_to_many_utils::convert_to_arclist;

pub struct DijkstraManyToMany {
    first_out: Vec<EdgeId>,
    arclist: Vec<(NodeId, Weight)>,
    target_indicies: ValidFlags<usize>, //stores the target index in the result matrix for all valid targets
    distances: ResettableArray<Weight>,
    queue: MinBinaryHeap<DijkstraState>,
    distance_table: Matrix<Weight>
}

impl DijkstraManyToMany {

    pub fn new(first_out: &Vec<EdgeId>, head: &Vec<NodeId>, weight: &Vec<Weight>) -> Self {
        let num_vertices = first_out.len() - 1;
        
        DijkstraManyToMany {
            first_out: first_out.clone(),
            arclist: convert_to_arclist(head, weight),
            target_indicies: ValidFlags::new(num_vertices, 0),
            distances: ResettableArray::new(num_vertices, INFINITY),
            queue: MinBinaryHeap::new(num_vertices),
            distance_table: Matrix::empty()
        }
    }

    fn target_selection(&mut self, targets: &[NodeId]) {
        self.target_indicies.reset();
        
        for (target_index, target) in targets.iter().enumerate() {
            self.target_indicies.set(*target as usize, target_index);
        }
    }

    fn calculate_distances(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        self.distance_table.resize(sources.len(), targets.len(), INFINITY);

        for (source_index, source) in sources.iter().enumerate() {
            self.distances.reset();
            self.queue.insert(DijkstraState {distance: 0, node_id: *source});
            let mut visited_targets = 0;

            while let Some(DijkstraState {distance: current_distance, node_id: current_node}) = self.queue.pop() {
                if self.target_indicies.is_valid(current_node as usize) {
                    let target_index = self.target_indicies.get_unchecked(current_node as usize);

                    self.distance_table.set(source_index, target_index, current_distance);
                    visited_targets += 1;

                    if visited_targets >= targets.len() {
                        self.queue.clear();
                        break;
                    }
                }

                unsafe {
                    let start = *self.first_out.get_unchecked(current_node as usize) as usize;
                    let end = *self.first_out.get_unchecked(current_node as usize + 1) as usize;
                    
                    for (target_node, edge_distance) in self.arclist.get_unchecked(start..end) {
                        let old_distance = self.distances[*target_node as usize];
                        let new_distance = current_distance + *edge_distance;
    
                        if old_distance > new_distance {
                            self.distances.set(*target_node as usize, new_distance);
                            self.queue.insert_or_decrease(DijkstraState {distance: new_distance, node_id: *target_node});
                        }
                    }
                }
            }
        }
    }
}

impl ManyToManyAlgorithm for DijkstraManyToMany {

    fn initialize(&mut self, _sources: &[NodeId], _targets: &[NodeId]) {
        // no initialization required
    }

    fn select(&mut self, _sources: &[NodeId], targets: &[NodeId]) {
        self.target_selection(targets);
    }

    fn query(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        self.calculate_distances(sources, targets);
    }

    fn get_distance_table(&mut self) -> &Matrix<Weight> {
        &self.distance_table
    }
}