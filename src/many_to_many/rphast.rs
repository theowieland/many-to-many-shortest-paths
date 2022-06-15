use crate::types::*;
use crate::data_structures::*;
use crate::many_to_many::many_to_many_algorithm::ManyToManyAlgorithm;
use crate::many_to_many::many_to_many_utils::*;
use crate::utils::binary_heap::MinBinaryHeap;

pub struct RPHAST {
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
    bwd_restricted_distances: Vec<Weight>,
    fwd_distances: Vec<Weight>,
    fwd_queue: MinBinaryHeap<RankDijkstraState>,
    distance_table: Matrix<Weight>
}

impl RPHAST {

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

        RPHAST {
            ranks: ranks.to_vec(),
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
            fwd_distances: vec![INFINITY; num_vertices],
            fwd_queue: MinBinaryHeap::new(num_vertices),
            distance_table: Matrix::empty()
        }
    }

    pub fn select_targets(&mut self, targets: &[NodeId]) {
        target_selection(&self.bwd_first_edge, &self.bwd_arclist, &mut self.bwd_restricted_nodes_translation, &mut self.dfs_discoverd, targets, &mut self.bwd_visited_nodes, &mut self.bwd_restricted_first_edge, &mut self.bwd_restricted_arclist);

        // resize bwd distance array to fit all nodes of the selected search space (reset happens during the query phase)
        self.bwd_restricted_distances.resize(self.bwd_visited_nodes.len(), INFINITY);
    }

    pub fn select_targets_decreasing_rank(&mut self, targets: &[NodeId]) {
        target_selection_decreasing_rank(
            &self.bwd_first_edge, 
            &self.bwd_arclist, 
            &self.ranks, 
            &mut self.bwd_restricted_nodes_translation, 
            &mut MinBinaryHeap::new(self.ranks.len()), 
            targets, 
            &mut self.bwd_visited_nodes, 
            &mut self.bwd_restricted_first_edge, 
            &mut self.bwd_restricted_arclist
        );

        // resize bwd distance array to fit all nodes of the selected search space (reset happens during the query phase)
        self.bwd_restricted_distances.resize(self.bwd_visited_nodes.len(), INFINITY);
    }

    pub fn query(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        self.distance_table.resize(sources.len(), targets.len(), INFINITY);

        // perform the forward search for all given source nodes
        for (source_index, source) in sources.iter().enumerate() {

            // perform forward search
            self.bwd_restricted_distances.iter_mut().for_each(|entry| *entry = INFINITY);

            self.fwd_queue.insert(RankDijkstraState {rank: self.ranks[*source as usize], node_id: *source});
            self.fwd_distances[*source as usize] = 0;

            while let Some(RankDijkstraState {rank: _current_rank, node_id: current_node}) = self.fwd_queue.pop() {
                unsafe {
                    let start = *self.fwd_first_edge.get_unchecked(current_node as usize) as usize;
                    let end = *self.fwd_first_edge.get_unchecked(current_node as usize + 1) as usize;

                    let current_distance = self.fwd_distances[current_node as usize];

                    for (target_node, edge_distance) in self.fwd_arclist.get_unchecked(start..end) {
                        let old_distance = self.fwd_distances.get_unchecked_mut(*target_node as usize);

                        if edge_distance + current_distance < *old_distance {
                            *old_distance = *edge_distance + current_distance;
                        }

                        if !self.fwd_queue.contains_unique_index(*target_node as usize) {
                            self.fwd_queue.insert(RankDijkstraState {rank: self.ranks[*target_node as usize], node_id: *target_node});
                        }
                    }

                    // copy distance to the bwd_restricted_distance array
                    if self.bwd_restricted_nodes_translation.is_valid(current_node as usize) {
                        self.bwd_restricted_distances[self.bwd_restricted_nodes_translation.get_unchecked(current_node as usize) as usize] = current_distance;
                    }

                    // reset old distances entry for future queries
                    *self.fwd_distances.get_unchecked_mut(current_node as usize) = INFINITY;
                }
            }
        
            // iterate over the backward nodes in order of discovery (highest to lowest rank)
            for node_index in 0..self.bwd_visited_nodes.len() {
                unsafe {
                    let start = *self.bwd_restricted_first_edge.get_unchecked(node_index) as usize;
                    let end  = *self.bwd_restricted_first_edge.get_unchecked(node_index + 1) as usize;     

                    let mut current_distance = *self.bwd_restricted_distances.get_unchecked(node_index);

                    // source_node is the local bwd_id, source_node is the node of the incoming arc into current_node
                    for (source_node, edge_distance) in self.bwd_restricted_arclist.get_unchecked(start..end) {
                        let source_distance = self.bwd_restricted_distances[*source_node as usize];
                        
                        if current_distance > *edge_distance + source_distance {
                            current_distance = *edge_distance + source_distance;
                        }
                    }
    
                    self.bwd_restricted_distances[node_index] = current_distance;
                }
            }    

            // update the distance_table
            unsafe {
                for (target_index, target) in targets.iter().enumerate() {
                    self.distance_table.set_unsafe(source_index, target_index, *self.bwd_restricted_distances.get_unchecked(self.bwd_restricted_nodes_translation[*target as usize] as usize))
                }
            }
        }
    }

    pub fn get_restricted_graph_size(&self) -> (usize, usize) {
        (self.bwd_restricted_first_edge.len() - 1, self.bwd_restricted_arclist.len())
    }
}

impl ManyToManyAlgorithm for RPHAST {

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

/// build the restricted backwards search space by only including nodes required to calculate the distance to any given target node
/// stores the restricted search space in the provided arguments
/// vertices reordered in dfs order
pub fn target_selection(
    bwd_first_edge: &Vec<EdgeId>,
    bwd_arclist: &Vec<(NodeId, Weight)>,
    bwd_restricted_nodes_translation: &mut ValidFlags<NodeId>,
    dfs_discoverd: &mut ValidFlags<bool>,
    targets: &[NodeId],
    bwd_visited_nodes: &mut Vec<NodeId>, // stores the restricted nodes in order of discovery
    bwd_restricted_first_edge: &mut Vec<EdgeId>,
    bwd_restricted_arclist: &mut Vec<(NodeId, Weight)>
) {
    bwd_restricted_nodes_translation.reset();

    depth_first_search(&bwd_first_edge, &bwd_arclist, dfs_discoverd, bwd_visited_nodes, &targets);

    // resize the reusable arrays
    bwd_restricted_first_edge.resize(bwd_visited_nodes.len() + 1, 0);
    bwd_restricted_arclist.clear();

    let mut first_edge_index = 0;

    for (node_index, node) in bwd_visited_nodes.iter().enumerate() {
        bwd_restricted_nodes_translation.set(*node as usize, node_index as NodeId);
    }

    // iterate over the visited nodes from highest to lowest rank
    for (node_index, node) in bwd_visited_nodes.iter().enumerate() {
        let start = bwd_first_edge[*node as usize] as usize;
        let end = bwd_first_edge[*node as usize + 1] as usize;

        bwd_restricted_first_edge[node_index] = first_edge_index;

        for (source_node, weight) in &bwd_arclist[start..end] {
            bwd_restricted_arclist.push((bwd_restricted_nodes_translation[*source_node as usize], *weight));
            first_edge_index += 1;
        }
    }

    bwd_restricted_first_edge[bwd_visited_nodes.len()] = first_edge_index;
}

/// build the restricted backwards search space by only including vertices required to calculate the distance to any given target
/// stores the restricted search space in tge given arguments
/// vertices reorderd in decreasing rank order
pub fn target_selection_decreasing_rank(
    bwd_first_edge: &Vec<EdgeId>,
    bwd_arclist: &Vec<(NodeId, Weight)>,
    ranks: &Vec<usize>,
    bwd_restricted_nodes_translation: &mut ValidFlags<NodeId>,
    queue: &mut MinBinaryHeap<RankDijkstraState>,
    targets: &[NodeId],
    bwd_visited_nodes: &mut Vec<NodeId>,
    bwd_restricted_first_edge: &mut Vec<EdgeId>,
    bwd_restricted_arclist: &mut Vec<(NodeId, Weight)>
) {
    bwd_restricted_nodes_translation.reset();
    bwd_visited_nodes.clear();

    for target in targets {
        queue.insert(RankDijkstraState {rank: ranks[*target as usize], node_id: *target});
    }

    while let Some(RankDijkstraState { rank: _, node_id: current_node }) = queue.pop() {
        let start = bwd_first_edge[current_node as usize] as usize;
        let end = bwd_first_edge[current_node as usize + 1] as usize;

        for (target_node, _weight) in &bwd_arclist[start..end] {
            if !queue.contains_unique_index(*target_node as usize) {
                queue.insert(RankDijkstraState {rank: ranks[*target_node as usize], node_id: *target_node});
            }
        }

        bwd_visited_nodes.push(current_node);
    }

    bwd_visited_nodes.reverse();
    bwd_restricted_first_edge.resize(bwd_visited_nodes.len() + 1, 0);
    bwd_restricted_arclist.clear();

    for (index, restricted_node) in bwd_visited_nodes.iter().enumerate() {
        bwd_restricted_nodes_translation.set(*restricted_node as usize, index as NodeId);
    }

    let mut first_arc_index = 0;
    for (index, restricted_node) in bwd_visited_nodes.iter().enumerate() {
        bwd_restricted_first_edge[index] = first_arc_index;

        let start = bwd_first_edge[*restricted_node as usize] as usize;
        let end = bwd_first_edge[*restricted_node as usize + 1] as usize;

        for (target_node, weight) in &bwd_arclist[start..end] {
            bwd_restricted_arclist.push((bwd_restricted_nodes_translation[*target_node as usize] as NodeId, *weight));
            first_arc_index += 1;
        }
    }

    bwd_restricted_first_edge[bwd_visited_nodes.len()] = first_arc_index;
}

/// performs a depth first search starting from the given vec of start nodes
/// returns all discovered nodes in order of discovery
pub fn depth_first_search(
    first_out: &Vec<EdgeId>, 
    arclist: &Vec<(NodeId, Weight)>, 
    discovered: &mut ValidFlags<bool>, 
    order_of_discovery: &mut Vec<NodeId>, // stores the discovered nodes
    start_nodes: &[NodeId]
) {
    let mut stack: Vec<(NodeId, usize, usize)> = Vec::new(); // node index and the index to continue to scan the edges from and the end index
    order_of_discovery.clear();

    discovered.reset();

    for start_node in start_nodes {
        if discovered.is_valid(*start_node as usize) {
            continue;
        }

        let start = first_out[*start_node as usize] as usize;
        let end = first_out[*start_node as usize + 1] as usize;
        stack.push((*start_node, start, end));

        discovered.set(*start_node as usize, true);
    
        while let Some((_last_element, start, end)) = stack.last() {
            let mut found_unvisited_node = false;
            
            for (target, _) in &arclist[*start..*end] {
                let len = stack.len();
                stack[len -1].1 += 1; // increase the arc start index for the next iterations

                if !discovered.is_valid(*target as usize) {
                    discovered.set(*target as usize, true);
                    let target_start = first_out[*target as usize] as usize;
                    let target_end = first_out[*target as usize + 1] as usize;

                    stack.push((*target, target_start, target_end));
                    found_unvisited_node = true;
                    break;
                }
            }
            
            if !found_unvisited_node {
                if let Some((last_node, _child_start_index, _child_end_index)) = stack.pop() {
                    order_of_discovery.push(last_node);
                }
            }
        }
    }
}