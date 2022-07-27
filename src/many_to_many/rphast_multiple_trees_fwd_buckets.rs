use crate::types::*;
use crate::data_structures::*;
use crate::bucket_containers::*;
use crate::utils::binary_heap::MinBinaryHeap;
use std::cmp::min;
use crate::many_to_many::many_to_many_algorithm::ManyToManyAlgorithm;
use crate::many_to_many::many_to_many_utils::*;
use super::advanced_bucket::populate_buckets_pruing_fwd;
use super::rphast::target_selection;
use super::rphast_utils::rphast_downward;

pub struct RPHASTMultipleTreesFwdBuckets {
    batch_size: usize,
    ranks: Vec<usize>, 
    fwd_first_edge: Vec<EdgeId>,
    fwd_arclist: Vec<(NodeId, Weight)>,
    bwd_first_edge: Vec<EdgeId>,
    bwd_arclist: Vec<(NodeId, Weight)>,
    terminals: Vec<usize>,
    fwd_buckets: BucketContainer<(usize, Weight)>,
    fwd_distances: Vec<Weight>,
    down_nodes: Vec<Vec<(NodeId, Weight)>>,
    up_nodes: Vec<Vec<(NodeId, Weight)>>,
    modified_up_nodes: Vec<NodeId>,
    fwd_queue: MinBinaryHeap<RankDijkstraState>,
    dfs_discovered: ValidFlags<bool>,
    bwd_visited_nodes: Vec<NodeId>,
    bwd_restricted_nodes_translation: ValidFlags<NodeId>,
    bwd_restricted_first_edge: Vec<EdgeId>,
    bwd_restricted_arclist: Vec<(NodeId, Weight)>,
    bwd_restricted_distances: Vec<Weight>, // flat 2d array that stores compute distances of the restricted search space

    distance_table: Matrix<Weight>
}

impl RPHASTMultipleTreesFwdBuckets {

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

        RPHASTMultipleTreesFwdBuckets {
            batch_size,
            ranks: ranks.to_owned(),
            fwd_first_edge: fwd_first_edge.to_owned(),
            fwd_arclist: convert_to_arclist(fwd_head, fwd_weight),
            bwd_first_edge: bwd_first_edge.to_owned(),
            bwd_arclist: convert_to_arclist(bwd_head, bwd_weight),
            terminals: Vec::new(),
            fwd_buckets: BucketContainer::new(num_vertices),
            fwd_distances: Vec::new(),
            down_nodes: vec![Vec::new(); num_vertices],
            up_nodes: vec![Vec::new(); num_vertices],
            modified_up_nodes: Vec::new(),
            fwd_queue: MinBinaryHeap::new(num_vertices),
            dfs_discovered: ValidFlags::new(num_vertices, false),
            bwd_visited_nodes: Vec::new(),
            bwd_restricted_nodes_translation: ValidFlags::new(num_vertices, INFINITY),
            bwd_restricted_first_edge: Vec::new(),
            bwd_restricted_arclist: Vec::new(),
            bwd_restricted_distances: Vec::new(),
            distance_table: Matrix::empty()
        }
    }

    pub fn print_buckets(&mut self) {
        for node_id in 0..self.ranks.len() {
            println!("{} {:?}", node_id, self.fwd_buckets.get_bucket(node_id));
        }
    }

    //build the restricted target search space
    pub fn select_targets(&mut self, targets: &[NodeId]) {
        target_selection(
            &self.bwd_first_edge, 
            &self.bwd_arclist, 
            &mut self.bwd_restricted_nodes_translation, 
            &mut self.dfs_discovered, 
            targets, 
            &mut self.bwd_visited_nodes, 
            &mut self.bwd_restricted_first_edge, 
            &mut self.bwd_restricted_arclist
        );

        //resize bwd_restricted_distances to fit current targets
        self.bwd_restricted_distances.resize(self.bwd_visited_nodes.len() * self.batch_size, INFINITY);
    }    

    pub fn query(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        self.distance_table.resize(sources.len(), targets.len(), INFINITY);

        batch_nodes(self.batch_size, sources, |batch_index, batched_nodes| {
            let bwd_restricted_distances = &mut self.bwd_restricted_distances;
            let bwd_restricted_nodes_translation = &mut self.bwd_restricted_nodes_translation;
            let batch_size = self.batch_size;

            bwd_restricted_distances.iter_mut().for_each(|entry| *entry = INFINITY);
            
            populate_buckets_pruing_fwd(
                &mut self.fwd_buckets,
                batched_nodes, 
                &mut self.terminals,
                &mut self.fwd_distances,
                &mut self.down_nodes,
                &mut self.up_nodes,
                &mut self.modified_up_nodes,
                &mut self.fwd_queue, 
                &self.fwd_first_edge, 
                &self.fwd_arclist,
                
                &mut self.bwd_first_edge, 
                &self.bwd_arclist, 
                &self.ranks, 
                
                |settled_node, _terminals, distances| {
                    if bwd_restricted_nodes_translation.is_valid(settled_node as usize) {
                        let bwd_restricted_id = bwd_restricted_nodes_translation.get_unchecked(settled_node as usize);
                        let start = bwd_restricted_id as usize * batch_size;
                        let end = min(start + batch_size, start + sources.len()); // if sources < batch_size we only need to copy remaining sources
            
                        bwd_restricted_distances[start..end].copy_from_slice(distances);
                    }
                }
            );

            // fix distances that were not copied during settling a vertex
            for (source_index, source_node) in batched_nodes.iter().enumerate() {
                if bwd_restricted_nodes_translation.is_valid(*source_node as usize) {
                    self.bwd_restricted_distances[(bwd_restricted_nodes_translation.get_unchecked(*source_node as usize) as usize) * self.batch_size + source_index] = 0;
                }
            }

            rphast_downward(
                &self.bwd_visited_nodes, 
                &self.bwd_restricted_first_edge, 
                &self.bwd_restricted_arclist, 
                &mut self.bwd_restricted_distances, 
                self.batch_size, 
                batched_nodes
            );

            // update the distance table
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

impl ManyToManyAlgorithm for RPHASTMultipleTreesFwdBuckets {

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

/// perform foward bucket fill algorithm and copy calculated distances to the bwd restricted distance array
pub fn fill_forward_buckets(
    sources: &[NodeId],
    batch_size: usize,
    bucket_container: &mut BucketContainer<(usize, Weight)>,
    queue: &mut MinBinaryHeap<RankDijkstraState>,
    bwd_restricted_distances: &mut Vec<Weight>, // flat distance vec that stores the discovered distances
    bwd_restricted_nodes_translation: &ValidFlags<NodeId>,
    distances: &mut Vec<Weight>, // used to copy distances from lower buckets to higher buckets
    down_nodes: &mut Vec<Vec<(NodeId, Weight)>>,
    terminals: &mut Vec<usize>,
    ranks: &Vec<usize>,
    fwd_first_edge: &Vec<EdgeId>,
    fwd_arclist: &Vec<(NodeId, Weight)>,

) {
    bucket_container.reset();
    bwd_restricted_distances.iter_mut().for_each(|entry| *entry = INFINITY);

    for (target_index, target) in sources.iter().enumerate() {
        if !queue.contains_unique_index(*target as usize) {
            queue.insert(RankDijkstraState {rank: ranks[*target as usize], node_id: *target});   
        }

        bucket_container.add_bucket_entry(*target as usize, (target_index, 0));
    }

    distances.resize(sources.len(), INFINITY);
    distances.iter_mut().for_each(|entry| *entry = INFINITY);

    while let Some(RankDijkstraState {rank: _, node_id: current_node}) = queue.pop() {
        let start_index = fwd_first_edge[current_node as usize] as usize;
        let end_index = fwd_first_edge[current_node as usize + 1] as usize;

        for (source_node, source_distance) in &fwd_arclist[start_index..end_index] {
            down_nodes[*source_node as usize].push((current_node, *source_distance));

            if !queue.contains_unique_index(*source_node as usize) {
                queue.insert(RankDijkstraState {rank: ranks[*source_node as usize], node_id: *source_node})
            }
        }
    
        //iterate over all nodes that the current_node is connected to (edge: current_node -> down_node)
        for (down_node, down_distance) in &down_nodes[current_node as usize] {
            for (target_index, target_distance) in bucket_container.get_bucket(*down_node as usize) {
                let current_distance = &mut distances[*target_index];
            
                if *current_distance == INFINITY {
                    terminals.push(*target_index);
                }

                *current_distance = min(*current_distance, target_distance + down_distance);
            }
        }

        // copy distance to the bwd restricted distance array in case node is in the bwd search space
        if bwd_restricted_nodes_translation.is_valid(current_node as usize) {
            let bwd_restricted_id = bwd_restricted_nodes_translation.get_unchecked(current_node as usize);
            let start = bwd_restricted_id as usize * batch_size;
            let end = min(start + batch_size, start + sources.len()); // if sources < batch_size we only need to copy remaining sources

            bwd_restricted_distances[start..end].copy_from_slice(distances);
        }

        terminals.iter().for_each(|index| {
            bucket_container.add_bucket_entry(current_node as usize, (*index, distances[*index]));

            distances[*index] = INFINITY; //reset distance array for future iterations
        });

        terminals.clear();
        down_nodes[current_node as usize].clear();
    }
}