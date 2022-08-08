use crate::utils::bucket_containers::BucketContainer;
use crate::utils::data_structures::{Matrix, ArrayStructure, EmptyValidFlags, LazyMaxHeap};
use crate::{types::*, utils::data_structures::ValidFlags};
use crate::graph_algorithms::DijkstraState;
use crate::many_to_many::many_to_many_algorithm::ManyToManyAlgorithm;
use crate::utils::binary_heap::MinBinaryHeap;
use super::many_to_many_utils::{convert_to_arclist, RankDijkstraState};

pub struct BucketManyToMany {
    ranks: Vec<usize>,
    fwd_first_edge: Vec<EdgeId>,
    fwd_arclist: Vec<(NodeId, Weight)>,
    fwd_distances: ValidFlags<Weight>,
    fwd_queue: MinBinaryHeap<RankDijkstraState>,
    bwd_first_edge: Vec<EdgeId>,
    bwd_arclist: Vec<(NodeId, Weight)>,
    bwd_distances: ValidFlags<Weight>,
    bwd_queue: MinBinaryHeap<DijkstraState>,
    bwd_buckets: BucketContainer<(usize, Weight)>,
    distance_table: Matrix<Weight>,
    fwd_distance_vec: Vec<Weight>
}

impl BucketManyToMany {

    pub fn new(
        fwd_first_edge: &EdgeIds, 
        fwd_head: &NodeIds, 
        fwd_weight: &Weights, 
        bwd_first_edge: &NodeIds, 
        bwd_head: &NodeIds, 
        bwd_weight: &Weights, 
        ranks: &Ranks
    ) -> Self {
        let num_vertices = fwd_first_edge.len() - 1;

        BucketManyToMany {
            ranks: ranks.to_owned(),
            fwd_first_edge: fwd_first_edge.to_owned(),
            fwd_arclist: convert_to_arclist(fwd_head, fwd_weight),
            fwd_distances: ValidFlags::<Weight>::new(num_vertices, INFINITY),
            fwd_queue: MinBinaryHeap::new(num_vertices),
            bwd_first_edge: bwd_first_edge.to_owned(),
            bwd_arclist: convert_to_arclist(bwd_head, bwd_weight),
            bwd_distances: ValidFlags::<Weight>::new(num_vertices, INFINITY),
            bwd_queue: MinBinaryHeap::new(num_vertices),
            bwd_buckets: BucketContainer::new(num_vertices),
            distance_table: Matrix::empty(),
            fwd_distance_vec: vec![INFINITY; num_vertices]
        }
    }

    pub fn init(&mut self) {
        self.fwd_queue.clear();
        self.bwd_queue.clear();
        self.fwd_distances.reset();
        self.bwd_buckets.reset();
    }

    pub fn select_targets(&mut self, targets: &[NodeId]) {
        initialize_bwd_buckets(
            &mut self.bwd_buckets, 
            targets, 
            &mut self.bwd_distances, 
            &mut self.bwd_queue, 
            &mut self.bwd_first_edge, 
            &mut self.bwd_arclist, 
            &mut self.fwd_first_edge, 
            &mut self.fwd_arclist
        );
    }

    pub fn select_targets_no_stall_on_demand(&mut self, targets: &[NodeId]) {
        initialize_bwd_buckets_no_stall_on_demand(
            &mut self.bwd_buckets, 
            targets, 
            &mut self.bwd_distances, 
            &mut self.bwd_queue, 
            &mut self.bwd_first_edge, 
            &mut self.bwd_arclist, 
        );
    }

    pub fn query(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        query_no_stopping(
            sources, 
            targets, 
            &self.fwd_first_edge, 
            &self.fwd_arclist, 
            &self.ranks, 
            &mut self.fwd_distance_vec, 
            &mut self.fwd_queue, 
            &mut self.bwd_buckets, 
            &mut self.distance_table
        );
    }

    pub fn get_average_number_of_bucket_entries(&mut self) -> (usize, usize) {
        let mut total_bucket_entry_count = 0;

        for modified_bucket in self.bwd_buckets.get_modified_buckets() {
            total_bucket_entry_count += self.bwd_buckets.get_bucket(*modified_bucket as usize).len();
        }
    
        (total_bucket_entry_count, self.bwd_buckets.get_modified_buckets().len())
    }
}

impl ManyToManyAlgorithm for BucketManyToMany {

    fn initialize(&mut self, _sources: &[NodeId], _targets: &[NodeId]) {
        self.init();
    }
    
    fn select(&mut self, _sources: &[NodeId], targets: &[NodeId]) {
        self.select_targets(targets);
    }

    fn query(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        self.query(sources, targets);
    }

    fn get_distance_table(&mut self) -> &Matrix<Weight> {
        &self.distance_table
    }
}

/// initializes the buckets of all nodes discovered by a backward search starting from the given target nodes
/// performs a full backward search for all targets individually
/// uses stall-on-demand to avoid unnecessary bucket entries
pub fn initialize_bwd_buckets(
    bucket_container: &mut BucketContainer<(usize, Weight)>,
    bucket_entry_nodes: &[NodeId],
    distances: &mut ValidFlags<Weight>,
    queue: &mut MinBinaryHeap<DijkstraState>, 
    bwd_first_edge: &Vec<NodeId>, 
    bwd_arclist: &Vec<(NodeId, Weight)>, 
    fwd_first_edge: &Vec<NodeId>,
    fwd_arclist: &Vec<(NodeId, Weight)>
) {
    for (entry_index, entry_node) in bucket_entry_nodes.iter().enumerate() { // start bwd search for each given entry_node
        distances.reset();
        distances.set(*entry_node as usize, 0);
        queue.insert(DijkstraState {distance: 0, node_id: *entry_node});

        while let Some(DijkstraState {distance: current_distance, node_id: current_node}) = queue.pop() {
            unsafe {
                // stall on demand, check if current node can be reached faster by adj fwd edge
                let fwd_start = *fwd_first_edge.get_unchecked(current_node as usize) as usize;
                let fwd_end = *fwd_first_edge.get_unchecked(current_node as usize + 1) as usize;
                let mut stall = false;

                // iterate over fwd edges leading away from the current node, to another node that is of higher rank but may have a shorter distance than the current node
                for (shortcut_node, shortcut_edge) in fwd_arclist.get_unchecked(fwd_start..fwd_end) {
                    if distances[*shortcut_node as usize] + *shortcut_edge <= current_distance {
                        stall = true;
                        break;
                    }
                }

                if !stall { //if not add the bucket entry and releax bwd edges from the current node
                    bucket_container.add_bucket_entry(current_node as usize, (entry_index, current_distance));

                    let start_index = *bwd_first_edge.get_unchecked(current_node as usize) as usize;
                    let end_index = *bwd_first_edge.get_unchecked(current_node as usize + 1) as usize;
    
                    for (source_node, edge_distance) in &bwd_arclist[start_index..end_index] {
                        let old_distance = distances[*source_node as usize];
                        
                        if edge_distance + current_distance < old_distance {
                            distances.set(*source_node as usize, edge_distance + current_distance);
                            queue.insert_or_decrease(DijkstraState {distance: current_distance + edge_distance, node_id: *source_node});
                        }
                    }
                }
            }
        }
    }
}

/// initializes the buckets of all nodes discovered by a backward search starting from the given target nodes
/// performs a full backward search for all targets individually
/// uses stall-on-demand to avoid unnecessary bucket entries
pub fn initialize_bwd_buckets_no_stall_on_demand(
    bucket_container: &mut BucketContainer<(usize, Weight)>,
    bucket_entry_nodes: &[NodeId],
    distances: &mut ValidFlags<Weight>,
    queue: &mut MinBinaryHeap<DijkstraState>, 
    bwd_first_edge: &Vec<NodeId>, 
    bwd_arclist: &Vec<(NodeId, Weight)>
) {
    for (entry_index, entry_node) in bucket_entry_nodes.iter().enumerate() { // start bwd search for each given entry_node
        distances.reset();
        distances.set(*entry_node as usize, 0);
        queue.insert(DijkstraState {distance: 0, node_id: *entry_node});

        while let Some(DijkstraState {distance: current_distance, node_id: current_node}) = queue.pop() {
            unsafe {
                bucket_container.add_bucket_entry(current_node as usize, (entry_index, current_distance));

                let start_index = *bwd_first_edge.get_unchecked(current_node as usize) as usize;
                let end_index = *bwd_first_edge.get_unchecked(current_node as usize + 1) as usize;
    
                for (source_node, edge_distance) in &bwd_arclist[start_index..end_index] {
                    let old_distance = distances[*source_node as usize];
                        
                    if edge_distance + current_distance < old_distance {
                        distances.set(*source_node as usize, edge_distance + current_distance);
                        queue.insert_or_decrease(DijkstraState {distance: current_distance + edge_distance, node_id: *source_node});
                    }
                }    
            }
        }
    }
}

pub fn query_no_stopping(
    sources: &[NodeId],
    targets: &[NodeId],
    fwd_first_edge: &Vec<EdgeId>,
    fwd_arclist: &Vec<(NodeId, Weight)>,
    ranks: &Vec<usize>,
    fwd_distances: &mut Vec<Weight>,
    fwd_queue: &mut MinBinaryHeap<RankDijkstraState>,
    bwd_buckets: &mut BucketContainer<(usize, Weight)>,
    distance_table: &mut Matrix<Weight>
) {
    distance_table.resize(sources.len(), targets.len(), INFINITY);

    for (source_index, source) in sources.iter().enumerate() {
        fwd_distances[*source as usize] = 0;
        fwd_queue.insert(RankDijkstraState {rank: ranks[*source as usize], node_id: *source});

        while let Some(RankDijkstraState {rank: _rank, node_id: current_node}) = fwd_queue.pop() {
            unsafe {
                let current_distance = *fwd_distances.get_unchecked(current_node as usize);
                distance_table.reduce_buckets(source_index, bwd_buckets.get_bucket(current_node as usize), current_distance);
                
                let start_index = *fwd_first_edge.get_unchecked(current_node as usize) as usize;
                let end_index = *fwd_first_edge.get_unchecked(current_node as usize + 1) as usize;
            
                for (target_node, edge_distance) in fwd_arclist.get_unchecked(start_index..end_index) {
                    let old_distance = fwd_distances.get_unchecked_mut(*target_node as usize);

                    if edge_distance + current_distance < *old_distance {
                        *old_distance = edge_distance + current_distance;
                    } 

                    if !fwd_queue.contains_unique_index(*target_node as usize) {
                        fwd_queue.insert(RankDijkstraState {rank: ranks[*target_node as usize], node_id: *target_node});
                    } 
                }

                *fwd_distances.get_unchecked_mut(current_node as usize) = INFINITY;
            }
        }
    }    
}

pub fn query(
    sources: &[NodeId],
    targets: &[NodeId],
    fwd_first_edge: &Vec<EdgeId>,
    fwd_arclist: &Vec<(NodeId, Weight)>,
    fwd_distances: &mut ValidFlags<Weight>,
    fwd_queue: &mut MinBinaryHeap<DijkstraState>,
    bwd_buckets: &mut BucketContainer<(usize, Weight)>,
    distance_table: &mut Matrix<Weight>
) {
    distance_table.resize(sources.len(), targets.len(), INFINITY);

    for (source_index, source) in sources.iter().enumerate() {
        fwd_distances.reset();
        fwd_distances.set(*source as usize, 0);
        fwd_queue.insert(DijkstraState {distance: 0, node_id: *source});
        
        let mut current_max_target_distance = INFINITY; // used for the stopping criteria, stores the currently calculated upper bound for any source-target distance
        let mut iteration_counter = 0; // used for the stopping criteria

        while let Some(DijkstraState {distance: current_distance, node_id: current_node}) = fwd_queue.pop() {
            // stopping criteria: abort if current_distance is larger than the current max of all targets reached from the current source
            // once the current_distance is lager than the max to all targets we cannot improve any distance because current_distance+(some bucket entry) would be larger than any previously discovered distance
            if current_distance > current_max_target_distance {
                fwd_queue.clear();
                break;
            }

            // scan all bucket entries and reduce the distance table entries
            for (target_index, target_distance) in bwd_buckets.get_bucket(current_node as usize) {
                if *target_distance > current_max_target_distance {
                    break;
                }

                if distance_table.reduce_value(source_index, *target_index, current_distance + *target_distance) {
                    
                    iteration_counter += 1;

                    // check if the max among all targets has changed and update the current_max_target_distance accordingly
                    if iteration_counter > 10000 { // only compute updated max after a given amount of iterations to reduce computational load
                        current_max_target_distance = 0; // temporarily set max to 0 and increase afterwards
                        for possible_max_index in 0..targets.len() { // iterate over all targets and increase max
                            if distance_table.get(source_index, possible_max_index) > current_max_target_distance {
                                current_max_target_distance = distance_table.get(source_index, possible_max_index);
                            }
                        }

                        iteration_counter = 0;
                    }
                }
            }

            let start_index = fwd_first_edge[current_node as usize] as usize;
            let end_index = fwd_first_edge[current_node as usize + 1] as usize;
            
            for (target_node, edge_distance) in fwd_arclist[start_index..end_index].iter() {
                let old_distance = fwd_distances[*target_node as usize];

                if edge_distance + current_distance < old_distance {
                    fwd_queue.insert_or_decrease(DijkstraState {distance: current_distance + edge_distance, node_id: *target_node});
                    fwd_distances.set(*target_node as usize, edge_distance + current_distance);
                }
            }
        }
    }    
}

//query by using a stopping criteria that counts the remaining targets that are above the current min key
pub fn query_with_counting_stopping_criteria(
    sources: &[NodeId],
    targets: &[NodeId],
    fwd_first_edge: &Vec<EdgeId>,
    fwd_arclist: &Vec<(NodeId, Weight)>,
    fwd_distances: &mut ValidFlags<Weight>,
    fwd_queue: &mut MinBinaryHeap<DijkstraState>,
    bwd_buckets: &mut BucketContainer<(usize, Weight)>,
    distance_table: &mut Matrix<Weight>,
    upgradeable_targets: &mut EmptyValidFlags 
) {
    distance_table.resize(sources.len(), targets.len(), INFINITY);

    upgradeable_targets.resize(targets.len()); //reset stopping criteria targets and resize

    for (source_index, source) in sources.iter().enumerate() {
        fwd_distances.reset();
        fwd_distances.set(*source as usize, 0);
        fwd_queue.insert(DijkstraState {distance: 0, node_id: *source});
        
        let mut num_upgradeable_targets = targets.len(); //used for the stopping criteria - e.g. we can stop once this counter reaches zero

        upgradeable_targets.reset();

        while let Some(DijkstraState {distance: current_distance, node_id: current_node}) = fwd_queue.pop() {
            //stopping criteria abort if no target can be upgraded
            if num_upgradeable_targets == 0 {
                fwd_queue.clear();
                break;
            }

            //scan all bucket entries and reduce the distance table entries
            for (target_index, target_distance) in bwd_buckets.get_bucket(current_node as usize) {
                distance_table.reduce_value(source_index, *target_index, current_distance + *target_distance);
            
                if !upgradeable_targets.is_valid(*target_index as usize) {
                    if distance_table.get(source_index, *target_index) <= current_distance {
                        upgradeable_targets.set_valid(*target_index as usize);
                        num_upgradeable_targets -= 1;
                    }
                }
            }
    
            let start_index = fwd_first_edge[current_node as usize] as usize;
            let end_index = fwd_first_edge[current_node as usize + 1] as usize;
            
            for (target_node, edge_distance) in fwd_arclist[start_index..end_index].iter() {
                let old_distance = fwd_distances[*target_node as usize];

                if edge_distance + current_distance < old_distance {
                    fwd_queue.insert_or_decrease(DijkstraState {distance: current_distance + edge_distance, node_id: *target_node});
                    fwd_distances.set(*target_node as usize, edge_distance + current_distance);
                }
            }
        }
    }    
}

pub fn query_stopping_criteria_heap(
    sources: &[NodeId],
    targets: &[NodeId],
    fwd_first_edge: &Vec<EdgeId>,
    fwd_arclist: &Vec<(NodeId, Weight)>,
    fwd_distances: &mut ValidFlags<Weight>,
    fwd_queue: &mut MinBinaryHeap<DijkstraState>,
    bwd_buckets: &mut BucketContainer<(usize, Weight)>,
    distance_table: &mut Matrix<Weight>,
    bucket_max_heap: &mut LazyMaxHeap<Weight>
) {
    distance_table.resize(sources.len(), targets.len(), INFINITY);

    bwd_buckets.sort_all_buckets(|(_id, dist)| *dist as usize);

    for (source_index, source) in sources.iter().enumerate() {
        fwd_distances.reset();
        fwd_distances.set(*source as usize, 0);
        fwd_queue.insert(DijkstraState {distance: 0, node_id: *source});

        bucket_max_heap.reset(targets.len(), INFINITY);

        while let Some(DijkstraState {distance: current_distance, node_id: current_node}) = fwd_queue.pop() {
            //stopping criteria: abort if current_distance is larger than the current max of all targets reached from the current source
            //once the current_distance is lager than the max to all targets we cannot improve any distance because current_distance+(some bucket entry) would be larger than any previously discovered distance
            let (current_max_index, current_max_distance) = bucket_max_heap.get_current_max();
            
            if current_distance > current_max_distance {
                fwd_queue.clear();
                break;
            }

            //scan all bucket entries and reduce the distance table entries
            for (target_index, target_distance) in bwd_buckets.get_bucket(current_node as usize) {
                if *target_distance > current_max_distance {
                    break;
                }

                if distance_table.reduce_value(source_index, *target_index, current_distance + *target_distance) {
                    bucket_max_heap.decrease_value(*target_index, current_distance + *target_distance);

                    if *target_index == current_max_index {
                        bucket_max_heap.recompute_max_element(); //recompute max on change
                    }
                }
            }

            let start_index = fwd_first_edge[current_node as usize] as usize;
            let end_index = fwd_first_edge[current_node as usize + 1] as usize;
            
            for (target_node, edge_distance) in &fwd_arclist[start_index..end_index] {
                let old_distance = fwd_distances[*target_node as usize];

                if edge_distance + current_distance < old_distance {
                    fwd_queue.insert_or_decrease(DijkstraState {distance: current_distance + edge_distance, node_id: *target_node});
                    fwd_distances.set(*target_node as usize, edge_distance + current_distance);
                }
            }
        }
    }    
}