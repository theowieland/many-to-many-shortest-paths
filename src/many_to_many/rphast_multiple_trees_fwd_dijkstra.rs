use crate::graph_algorithms::DijkstraState;
use crate::types::*;
use crate::many_to_many::many_to_many_algorithm::ManyToManyAlgorithm;
use crate::many_to_many::many_to_many_utils::*;
use crate::many_to_many::rphast::target_selection;
use crate::utils::binary_heap::MinBinaryHeap;
use crate::utils::data_structures::Matrix;
use crate::utils::data_structures::ResettableArray;
use crate::utils::data_structures::ValidFlags;

use super::rphast_utils::rphast_downward;
use super::rphast_utils::rphast_forward_dijkstra;

pub struct RPHASTMultipleTreesFwdDijkstra {
    batch_size: usize,
    dfs_discoverd: ValidFlags<bool>,
    fwd_first_edge: Vec<EdgeId>,
    fwd_arclist: Vec<(NodeId, Weight)>,
    bwd_first_edge: Vec<EdgeId>,
    bwd_arclist: Vec<(NodeId, Weight)>,
    bwd_visited_nodes: Vec<NodeId>,
    bwd_restricted_nodes_translation: ValidFlags<NodeId>,
    bwd_restricted_first_edge: Vec<EdgeId>,
    bwd_restricted_arclist: Vec<(NodeId, Weight)>,
    bwd_restricted_distances: Vec<Weight>, //flat 2d distance array
    fwd_distances: ResettableArray<Weight>,
    fwd_queue: MinBinaryHeap<DijkstraState>,
    distance_table: Matrix<Weight>
}

impl RPHASTMultipleTreesFwdDijkstra {

    pub fn new(
        batch_size: usize, 
        fwd_first_edge: &EdgeIds, 
        fwd_head: &NodeIds, 
        fwd_weight: &Weights,
        bwd_first_edge: &EdgeIds, 
        bwd_head: &NodeIds, 
        bwd_weight: &Weights
    ) -> Self {
        let num_vertices = fwd_first_edge.len() - 1;

        RPHASTMultipleTreesFwdDijkstra {
            batch_size,
            dfs_discoverd: ValidFlags::new(num_vertices, false),
            fwd_first_edge: fwd_first_edge.to_owned(),
            fwd_arclist: convert_to_arclist(fwd_head, fwd_weight),
            bwd_first_edge: bwd_first_edge.to_owned(),
            bwd_arclist: convert_to_arclist(bwd_head, bwd_weight),
            bwd_visited_nodes: Vec::new(),
            bwd_restricted_nodes_translation: ValidFlags::new(num_vertices, INFINITY), // translate from node id to local id
            bwd_restricted_first_edge: Vec::new(),
            bwd_restricted_arclist: Vec::new(),
            bwd_restricted_distances: Vec::new(),
            fwd_distances: ResettableArray::new(num_vertices, INFINITY),
            fwd_queue: MinBinaryHeap::new(num_vertices),
            distance_table: Matrix::empty()
        }
    }

    /// build the restricted target search space
    pub fn select_targets(&mut self, targets: &[NodeId]) {
        target_selection(&self.bwd_first_edge, &self.bwd_arclist, &mut self.bwd_restricted_nodes_translation, &mut self.dfs_discoverd, targets, &mut self.bwd_visited_nodes, &mut self.bwd_restricted_first_edge, &mut self.bwd_restricted_arclist);

        // resize bwd_restricted_distances to fit current targets
        self.bwd_restricted_distances.resize(self.bwd_visited_nodes.len() * self.batch_size, INFINITY);
    }

    pub fn query(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        self.distance_table.resize(sources.len(), targets.len(), INFINITY);
        

        batch_nodes(self.batch_size, sources, |batch_index, batched_nodes| {
            
            // clear the bwd restricted distances from previous queries
            self.bwd_restricted_distances.iter_mut().for_each(|entry| *entry = INFINITY);

            rphast_forward_dijkstra(
                batched_nodes, 
                &mut self.fwd_queue, 
                &self.fwd_first_edge, 
                &self.fwd_arclist, 
                &mut self.fwd_distances, 
                self.batch_size, 
                &self.bwd_restricted_nodes_translation, 
                &mut self.bwd_restricted_distances
            );

            rphast_downward(
                &self.bwd_visited_nodes, 
                &self.bwd_restricted_first_edge, 
                &self.bwd_restricted_arclist, 
                &mut self.bwd_restricted_distances, 
                self.batch_size, 
                batched_nodes
            );
            
            //update the distance table
            unsafe {
                for (batch_source_index, _batch_source) in batched_nodes.iter().enumerate() {
                    for (target_index, target) in targets.iter().enumerate() {
                        self.distance_table.set_unsafe(batch_source_index + (batch_index * self.batch_size), target_index, self.bwd_restricted_distances[(self.bwd_restricted_nodes_translation[*target as usize] as usize) * self.batch_size + batch_source_index]);
                    }
                }
            }
        });
    }
}

impl ManyToManyAlgorithm for RPHASTMultipleTreesFwdDijkstra {

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