use crate::types::*;
use crate::utils::binary_heap::HeapElement;
use crate::utils::binary_heap::MinBinaryHeap;
use crate::utils::data_structures::ArrayStructure;
use crate::utils::data_structures::ValidFlags;

use super::contraction_hierarchies::ContractionHierarchiesGraph;

const WITNESS_SEARCH_TIMEOUT: usize = 50;

#[derive(Copy, Clone, Eq, PartialEq, Debug, Ord, PartialOrd)]
pub struct WitnessSearchState {
    
    distance: Weight,
    node_id: NodeId,
}

impl HeapElement for WitnessSearchState {

    fn unique_index(&self) -> usize {
        self.node_id as usize
    }
}

pub struct WitnessSearch {
    start: NodeId,
    node_to_contract: NodeId,
    distance_table: ValidFlags<Weight>,
    to_visit: MinBinaryHeap<WitnessSearchState>,
}

impl WitnessSearch {

    pub fn new(num_nodes: usize) -> WitnessSearch {
        let to_visit: MinBinaryHeap<WitnessSearchState> = MinBinaryHeap::new(num_nodes);

        WitnessSearch {
            start: 0,
            node_to_contract: 0,
            distance_table: ValidFlags::new(num_nodes, INFINITY),
            to_visit,
        }
    }

    pub fn prepare(&mut self, start: NodeId, node_to_contract: NodeId) {
        self.distance_table.reset();
        self.to_visit.clear();
        self.start = start;
        self.node_to_contract = node_to_contract;

        self.to_visit.insert(WitnessSearchState {distance: 0, node_id: start});
        self.distance_table.set(start as usize, 0);
    }

    pub fn continue_search(&mut self, goal: NodeId, max_weight: Weight, graph: &ContractionHierarchiesGraph) -> Weight {
        if goal == self.start {
            return 0;
        }

        if self.distance_table[goal as usize] < max_weight {
            return self.distance_table[goal as usize]; 
        }

        let mut search_depth = WITNESS_SEARCH_TIMEOUT;

        while let Some(WitnessSearchState {distance, node_id}) = self.to_visit.pop() {
            search_depth -= 1;

            for (target_node, edge_weight) in &graph.remaining_out_edges[node_id as usize] {
                if *target_node == self.node_to_contract {
                    continue;
                }

                let new_distance = distance + edge_weight;
                let old_distance = self.distance_table[*target_node as usize];

                if new_distance < old_distance { // found new shortest path
                    self.to_visit.insert_or_decrease(WitnessSearchState {distance: new_distance, node_id: *target_node});
                    self.distance_table.set(*target_node as usize, new_distance);
                }    
            }

            if node_id == goal {
                return distance;
            }

            if distance > max_weight {
                return max_weight + 1;
            }

            if search_depth == 0 {
                return max_weight + 1;
            }
        }

        self.distance_table[goal as usize]
    }
}