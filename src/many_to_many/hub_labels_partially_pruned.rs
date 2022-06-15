use crate::graph_algorithms::DijkstraState;
use crate::types::*;
use crate::data_structures::{ValidFlags, Matrix, ArrayStructure};
use crate::many_to_many::many_to_many_algorithm::ManyToManyAlgorithm;
use crate::many_to_many::many_to_many_utils::*;
use crate::utils::binary_heap::MinBinaryHeap;

pub struct HubLabelsPartiallyPruned {
    fwd_first_edge: Vec<EdgeId>,
    fwd_arclist: Vec<(NodeId, Weight)>,
    bwd_first_edge: Vec<EdgeId>,
    bwd_arclist: Vec<(NodeId, Weight)>,
    fwd_labels: Vec<Vec<(NodeId, Weight)>>,
    bwd_labels: Vec<Vec<(NodeId, Weight)>>,
    queue: MinBinaryHeap<DijkstraState>,
    distances: ValidFlags<Weight>,
    distance_table: Matrix<Weight>,
}

impl HubLabelsPartiallyPruned {

    pub fn new(
        fwd_first_edge: &Vec<EdgeId>, 
        fwd_head: &Vec<NodeId>, 
        fwd_weight: &Vec<Weight>, 
        bwd_first_edge: &Vec<EdgeId>, 
        bwd_head: &Vec<NodeId>, 
        bwd_weight: &Vec<Weight>
    ) -> Self {
        let num_vertices = fwd_first_edge.len() - 1;

        HubLabelsPartiallyPruned {
            fwd_first_edge: fwd_first_edge.to_owned(),
            fwd_arclist: convert_to_arclist(fwd_head, fwd_weight),
            bwd_first_edge: bwd_first_edge.to_owned(),
            bwd_arclist: convert_to_arclist(bwd_head, bwd_weight),
            fwd_labels: Vec::new(),
            bwd_labels: Vec::new(),
            queue: MinBinaryHeap::new(num_vertices),
            distances: ValidFlags::new(num_vertices, INFINITY),
            distance_table: Matrix::empty(),
        }
    }

    pub fn calculate_labels(&mut self, sources: &[NodeId], targets: &[NodeId]) { 
        self.fwd_labels.resize(sources.len(), Vec::new());
        self.bwd_labels.resize(targets.len(), Vec::new());

        self.fwd_labels.iter_mut().for_each(|hubs| hubs.clear());
        self.bwd_labels.iter_mut().for_each(|hubs| hubs.clear());    

        for (source_index, source) in sources.iter().enumerate() {
            HubLabelsPartiallyPruned::fill_labels(
                source_index,
                *source, 
                &mut self.fwd_labels, 
                &mut self.queue, 
                &mut self.distances, 
                &self.fwd_first_edge, 
                &self.fwd_arclist, 
                &self.bwd_first_edge, 
                &self.bwd_arclist
            );
        }

        for (target_index, target) in targets.iter().enumerate() {
            HubLabelsPartiallyPruned::fill_labels(
                target_index,
                *target, 
                &mut self.bwd_labels, 
                &mut self.queue, 
                &mut self.distances, 
                &self.bwd_first_edge, 
                &self.bwd_arclist, 
                &self.fwd_first_edge, 
                &self.fwd_arclist
            );
        }
    }

    /// calculate hub labels for the given node
    /// runs a dijkstra search upwards - with full stall on demand applied
    fn fill_labels(
        label_index: usize, // index in the labels array to store the labels
        start_node: NodeId, 
        labels: &mut Vec<Vec<(NodeId, Weight)>>, 
        queue: &mut MinBinaryHeap<DijkstraState>, 
        distances: &mut ValidFlags<Weight>, 
        first_edge: &Vec<EdgeId>, 
        arclist: &Vec<(NodeId, Weight)>, 
        prune_first_edge: &Vec<EdgeId>, 
        prune_arclist: &Vec<(NodeId, Weight)>
    ) {
        distances.reset();
        distances.set(start_node as usize, 0);

        queue.insert(DijkstraState {distance: 0, node_id: start_node});

        while let Some(DijkstraState {distance: current_distance, node_id: current_node}) = queue.pop() {
            unsafe {
                //check for prune edges
                let prune_start = *prune_first_edge.get_unchecked(current_node as usize) as usize;
                let prune_end = *prune_first_edge.get_unchecked(current_node as usize + 1) as usize;
                let mut prune = false;
                for (prune_node, prune_distance) in prune_arclist.get_unchecked(prune_start..prune_end) {
                    if distances[*prune_node as usize] + *prune_distance < current_distance { //found a shorter path to the current node
                        prune = true;
                        break;
                    }
                }

                if !prune {
                    labels.get_unchecked_mut(label_index).push((current_node, current_distance));

                    let start = *first_edge.get_unchecked(current_node as usize) as usize;
                    let end = *first_edge.get_unchecked(current_node as usize + 1) as usize;

                    for (target_node, arc_distance) in arclist.get_unchecked(start..end) {
                        let target_current_distance = distances[*target_node as usize];

                        if current_distance + *arc_distance < target_current_distance {
                            distances.set(*target_node as usize, current_distance + *arc_distance);

                            queue.insert_or_decrease(DijkstraState {distance: distances[*target_node as usize], node_id: *target_node});
                        }
                    }
                }
            }
        }

        labels[label_index].sort_by_key(|(hub, _distance)| *hub);
    }

    pub fn get_average_label_size(&self, sources: &[NodeId], targets: &[NodeId]) -> (usize, usize) {
        let mut total_fwd_labels = 0;
        let mut total_bwd_labels = 0;

        for source in sources {
            total_fwd_labels += self.fwd_labels[*source as usize].len();
        }

        for target in targets {
            total_bwd_labels += self.bwd_labels[*target as usize].len();
        }

        ((total_fwd_labels as f64 / sources.len() as f64) as usize, (total_bwd_labels as f64 / targets.len() as f64) as usize)
    }

    /// calculates the shortest distance together with the used hub node id
    /// returns (hub, distance)
    pub fn query_shortest_path(fwd_labels: &[(NodeId, Weight)], bwd_labels: &[(NodeId, Weight)]) -> (NodeId, Weight) {
        let mut used_hub = 0;
        let mut shortest_distance = INFINITY; //currently discovered shortest path between the fwd_labels and bwd_labels

        let mut fwd_index = 0;
        let mut bwd_index = 0;

        while fwd_index < fwd_labels.len() && bwd_index < bwd_labels.len() {
            let (fwd_hub, fwd_distance) = fwd_labels[fwd_index];
            let (bwd_hub, bwd_distance) = bwd_labels[bwd_index];

            if fwd_hub == bwd_hub {
                if fwd_distance + bwd_distance < shortest_distance {
                    shortest_distance = fwd_distance + bwd_distance;
                    used_hub = fwd_hub;
                }

                fwd_index += 1;
                bwd_index += 1;
            }
            else if fwd_hub < bwd_hub {
                fwd_index += 1;
            }
            else {
                bwd_index += 1;
            }
        }
        
        (used_hub, shortest_distance)
    }

    pub fn query(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        HubLabelsPartiallyPruned::query_all(sources, targets, &mut self.distance_table, &self.fwd_labels, &self.bwd_labels);
    }

    /// calculate shortest path distances for all source target pairs
    /// stores the calculated shortest path distances in the given distance_table matrix
    pub fn query_all(
        sources: &[NodeId], 
        targets: &[NodeId],
        distance_table: &mut Matrix<Weight>,
        fwd_labels: &Vec<Vec<(NodeId, Weight)>>,
        bwd_labels: &Vec<Vec<(NodeId, Weight)>>
    ) {
        distance_table.resize(sources.len(), targets.len(), INFINITY);

        unsafe {
            for (source_index, _source) in sources.iter().enumerate() {
                let fwd_labels = fwd_labels.get_unchecked(source_index);

                for (target_index, _target) in targets.iter().enumerate() {
                    let bwd_labels = bwd_labels.get_unchecked(target_index);

                    let (_hub, distance) = HubLabelsPartiallyPruned::query_shortest_path(fwd_labels, bwd_labels);
                    distance_table.set(source_index, target_index, distance);
                }
            }
        }
    }
}

impl ManyToManyAlgorithm for HubLabelsPartiallyPruned {

    fn initialize(&mut self, _sources: &[NodeId], _targets: &[NodeId]) {
        // no initialization required
    }

    fn select(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        self.calculate_labels(sources, targets);
    }

    fn query(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        self.query(sources, targets);
    }

    fn get_distance_table(&mut self) -> &Matrix<Weight> {
        &self.distance_table
    }
}