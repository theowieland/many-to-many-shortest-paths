use crate::types::*;
use crate::graph_representation::Graph;
use crate::utils::binary_heap::{HeapElement, MinBinaryHeap};

#[derive(Copy, Clone, Eq, PartialEq, Debug, Ord, PartialOrd)]
pub struct DijkstraState {
    pub distance: Weight,
    pub node_id: NodeId,
}

impl HeapElement for DijkstraState {

    fn unique_index(&self) -> usize {
        self.node_id as usize
    }
}

pub fn dijkstra(start: NodeId, goal: NodeId, graph: &impl Graph) -> Option<Weight> {
    let mut to_visit: MinBinaryHeap<DijkstraState> = MinBinaryHeap::new(graph.node_ids().len());
    to_visit.insert(DijkstraState {
        distance: 0,
        node_id: start,
    });

    let mut distance_table: Vec<Weight> = vec![INFINITY; graph.node_ids().len()];
    distance_table[start as usize] = 0;

    while let Some(DijkstraState {distance, node_id}) = to_visit.pop() {
        if node_id == goal {
            return Some(distance_table[goal as usize]);
        }

        for edge_id in graph.edge_ids(node_id) {
            let target_node = graph.target_id(node_id, edge_id);
            let new_distance = distance + graph.weight(node_id, edge_id);

            if distance_table[target_node as usize] > new_distance {
                distance_table[target_node as usize] = new_distance;

                to_visit.insert_or_decrease(DijkstraState {distance: new_distance, node_id: target_node});
            }
        }
    }

    None
}