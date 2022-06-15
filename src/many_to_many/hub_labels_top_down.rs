use crate::types::*;
use crate::data_structures::Matrix;
use crate::many_to_many::many_to_many_algorithm::ManyToManyAlgorithm;
use crate::many_to_many::many_to_many_utils::*;
use crate::utils::binary_heap::MinBinaryHeap;

use super::hub_labels_partially_pruned::HubLabelsPartiallyPruned;

pub struct HubLabelsTopDown {
    ranks: Vec<usize>,
    fwd_first_edge: Vec<EdgeId>,
    fwd_arclist: Vec<(NodeId, Weight)>,
    bwd_first_edge: Vec<EdgeId>,
    bwd_arclist: Vec<(NodeId, Weight)>,
    fwd_labels: Vec<Vec<(NodeId, Weight)>>,
    bwd_labels: Vec<Vec<(NodeId, Weight)>>,
    modified_labels: Vec<NodeId>, // contains the indices of previously modified labels - fwd and bwd labels are modified for the same vertices
    queue: MinBinaryHeap<RankDijkstraState>,
    merge_buffer: Vec<Vec<(NodeId, Weight)>>,
    retainable_hubs: Vec<bool>,
    distance_table: Matrix<Weight>,
}

/// compute hub labels for a restricted rphast search space top down - fully pruned
impl HubLabelsTopDown {

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

        HubLabelsTopDown {
            ranks: ranks.to_vec(),
            fwd_first_edge: fwd_first_edge.to_owned(),
            fwd_arclist: convert_to_arclist(fwd_head, fwd_weight),
            bwd_first_edge: bwd_first_edge.to_owned(),
            bwd_arclist: convert_to_arclist(bwd_head, bwd_weight),
            fwd_labels: vec![Vec::new(); num_vertices],
            bwd_labels: vec![Vec::new(); num_vertices],
            modified_labels: Vec::new(),
            queue: MinBinaryHeap::new(num_vertices),
            merge_buffer: Vec::new(),
            retainable_hubs: Vec::new(),
            distance_table: Matrix::empty(),
        }
    }

    pub fn calculate_labels(&mut self, sources: &[NodeId], targets: &[NodeId]) {

        // clear previously initialized labels
        for modified_label in &self.modified_labels {
            self.fwd_labels[*modified_label as usize].clear();
            self.bwd_labels[*modified_label as usize].clear();
        }
        self.modified_labels.clear();

        let mut discovered_nodes: Vec<NodeId> = Vec::new();
        
        HubLabelsTopDown::fill_hub_labels_top_down(
            sources,
            targets,
            &mut self.fwd_labels,
            &mut self.bwd_labels,
            &mut self.modified_labels,
            &mut self.queue,
            &mut discovered_nodes,
            &mut self.merge_buffer,
            &mut self.retainable_hubs,
            &self.ranks,
            &self.fwd_first_edge,
            &self.fwd_arclist,
            &self.bwd_first_edge,
            &self.bwd_arclist
        ); 
    }

    /// initialize all labels for all nodes in the used search space
    pub fn fill_hub_labels_top_down(
        sources: &[NodeId], 
        targets: &[NodeId],
        fwd_labels: &mut Vec<Vec<(NodeId, Weight)>>, 
        bwd_labels: &mut Vec<Vec<(NodeId, Weight)>>,
        modified_labels: &mut Vec<NodeId>,
        queue: &mut MinBinaryHeap<RankDijkstraState>, 
        discovered_nodes: &mut Vec<NodeId>,
        merge_buffer: &mut Vec<Vec<(NodeId, Weight)>>,
        retainable_hubs: &mut Vec<bool>,
        ranks: &Vec<usize>,
        fwd_first_edge: &Vec<EdgeId>, 
        fwd_arclist: &Vec<(NodeId, Weight)>, 
        bwd_first_edge: &Vec<EdgeId>, 
        bwd_arclist: &Vec<(NodeId, Weight)>
    ) {
        discovered_nodes.clear();

        for source_node in sources {
            if !queue.contains_unique_index(*source_node as usize) {
                queue.insert(RankDijkstraState {rank: ranks[*source_node as usize], node_id: *source_node});
            }
        }

        while let Some(RankDijkstraState {rank: _current_rank, node_id: current_node}) = queue.pop() {
            discovered_nodes.push(current_node);

            let start = fwd_first_edge[current_node as usize] as usize;
            let end = fwd_first_edge[current_node as usize + 1] as usize;

            for (target_node, _edge_distance) in &fwd_arclist[start..end] {
                if !queue.contains_unique_index(*target_node as usize) {
                    queue.insert(RankDijkstraState {rank: ranks[*target_node as usize], node_id: *target_node});
                }
            }
        }

        for target_node in targets {
            if !queue.contains_unique_index(*target_node as usize) {
                queue.insert(RankDijkstraState {rank: ranks[*target_node as usize], node_id: *target_node});
            }
        }

        while let Some(RankDijkstraState {rank: _current_rank, node_id: current_node}) = queue.pop() {
            discovered_nodes.push(current_node);

            let start = bwd_first_edge[current_node as usize] as usize;
            let end = bwd_first_edge[current_node as usize + 1] as usize;

            for (target_node, _edge_distance) in &bwd_arclist[start..end] {
                if !queue.contains_unique_index(*target_node as usize) {
                    queue.insert(RankDijkstraState {rank: ranks[*target_node as usize], node_id: *target_node});
                }
            }
        }

        discovered_nodes.sort_by_key(|node_id| ranks[*node_id as usize]); //sort from highest to lowest rank
        discovered_nodes.dedup(); //remove duplicate node ids

        //iterate over nodes top to bottom (high rank to low rank)
        for node in discovered_nodes.iter().rev() {
            modified_labels.push(*node);

            k_way_merge(*node, fwd_labels, fwd_first_edge, fwd_arclist, merge_buffer);
            k_way_merge(*node, bwd_labels, bwd_first_edge, bwd_arclist, merge_buffer);

            retainable_hubs.resize(fwd_labels[*node as usize].len(), false);
            for hub_index in 0..fwd_labels[*node as usize].len() {
                let (hub, _hub_distance) = fwd_labels[*node as usize][hub_index];

                if hub == *node {
                    retainable_hubs[hub_index] = true;
                }
                else {
                    let (used_hub, _used_hub_distance) = HubLabelsPartiallyPruned::query_shortest_path(&fwd_labels[*node as usize], &bwd_labels[hub as usize]);
                    
                    if used_hub == hub {
                        retainable_hubs[hub_index] = true;
                    }
                    else {
                        retainable_hubs[hub_index] = false;
                    }
                }
            }

            let mut index = 0;
            fwd_labels[*node as usize].retain(|_| { index+=1; retainable_hubs[index-1] } );

            retainable_hubs.resize(bwd_labels[*node as usize].len(), false);
            for hub_index in 0..bwd_labels[*node as usize].len() {
                let (hub, _hub_distance) = bwd_labels[*node as usize][hub_index];

                if hub == *node {
                    retainable_hubs[hub_index] = true;
                }
                else {
                    let (used_hub, _used_hub_distance) = HubLabelsPartiallyPruned::query_shortest_path(&fwd_labels[hub as usize], &bwd_labels[*node as usize]);
                    
                    if used_hub == hub {
                        retainable_hubs[hub_index] = true;
                    }
                    else {
                        retainable_hubs[hub_index] = false;
                    }
                }
            }

            let mut index = 0;
            bwd_labels[*node as usize].retain(|_| { index+=1; retainable_hubs[index-1] } );
        }    
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
}

impl ManyToManyAlgorithm for HubLabelsTopDown {

    fn initialize(&mut self, _sources: &[NodeId], _targets: &[NodeId]) {
        // no initilization required
    }

    fn select(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        self.calculate_labels(sources, targets);
    }

    fn query(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        self.distance_table.resize(sources.len(), targets.len(), INFINITY);
    
        unsafe {
            for (source_index, source) in sources.iter().enumerate() {
                let fwd_labels = self.fwd_labels.get_unchecked(*source as usize);
    
                for (target_index, target) in targets.iter().enumerate() {
                    let bwd_labels = self.bwd_labels.get_unchecked(*target as usize);
    
                    let (_hub, distance) = HubLabelsPartiallyPruned::query_shortest_path(fwd_labels, bwd_labels);
                    self.distance_table.set(source_index, target_index, distance);
                }
            }
        }
    }

    fn get_distance_table(&mut self) -> &Matrix<Weight> {
        &self.distance_table
    }
}

/// given two sorted (by ascending NodeId, each NodeId has to be unique) vecs merges them so that the result vec 
/// adds a weight offset to each hub distance
/// contains every unique NodeId with shortest weight in ascending order
fn merge(first_offset: Weight, first: &Vec<(NodeId, Weight)>, second_offset: Weight, second: &Vec<(NodeId, Weight)>, result: &mut Vec<(NodeId, Weight)>) {
    let mut first_index = 0;
    let mut second_index = 0;

    while first_index < first.len() && second_index < second.len() {
        let (first_id, first_distance) = first[first_index];
        let (second_id, second_distance) = second[second_index];

        if first_id < second_id {
            result.push((first_id, first_offset + first_distance));
            first_index += 1;
        }
        else if second_id < first_id {
            result.push((second_id, second_offset + second_distance));
            second_index += 1;
        }
        else { //ids are equal -> exclusively add node id with shortest distance
            if first_offset + first_distance <= second_offset + second_distance {
                result.push((first_id, first_offset + first_distance));
            }
            else {
                result.push((second_id, second_offset + second_distance));
            }

            first_index += 1;
            second_index += 1;
        }
    }

    //add remaining values in case one of the both vecs has been scanned to the end
    for index in first_index..first.len() {
        result.push((first[index].0, first_offset + first[index].1));
    }  

    for index in second_index..second.len() {
        result.push((second[index].0, second_offset + second[index].1));
    }          
}

fn k_way_merge(
    node_id: NodeId,
    labels: &mut Vec<Vec<(NodeId, Weight)>>,
    first_out: &Vec<EdgeId>,
    arclist: &Vec<(NodeId, Weight)>,

    merge_buffer: &mut Vec<Vec<(NodeId, Weight)>>
) {
    let start = first_out[node_id as usize] as usize;
    let end = first_out[node_id as usize + 1] as usize;

    // resize merge_buffer to hold enough vecs to perform merges
    merge_buffer.resize(end - start + 1, Vec::new()); // number of merges for k arrays: (end - start - 1) +1 added for adding the own node id with distance 0
    merge_buffer.iter_mut().for_each(|entry| entry.clear());

    let mut merge_buffer_result_index = 0; // index of first unused vec in merge_buffer

    // initial merges
    for initial_merge in 0..((end - start) / 2) {
        let (first, first_distance) = arclist[start + (initial_merge * 2)];
        let (second, second_distance) = arclist[start + (initial_merge * 2) + 1];

        merge(first_distance, &labels[first as usize], second_distance, &labels[second as usize], &mut merge_buffer[merge_buffer_result_index]);
    
        merge_buffer_result_index += 1;
    }

    if (end - start) % 2 != 0 { // uneven initial length -> add missing vec to merge_buffer
        let (target, target_distance) = arclist[end - 1];
        merge_buffer[merge_buffer_result_index].append(&mut labels[target as usize].iter().map(|(hub, distance)| (*hub, *distance + target_distance)).collect());
       
        merge_buffer_result_index += 1;
    }

    // add own id as hub
    merge_buffer[merge_buffer_result_index].push((node_id, 0));
    merge_buffer_result_index += 1;

    let mut first_index_to_process = 0;

    while merge_buffer_result_index - first_index_to_process > 2 {
        let (to_process, result) = merge_buffer.split_at_mut(merge_buffer_result_index);

        let first = &to_process[first_index_to_process];
        let second = &to_process[first_index_to_process + 1];
        first_index_to_process += 2;

        merge(0, first, 0, second, &mut result[0]);
        merge_buffer_result_index += 1;
    }

    // directly merge result into labels array, last two remaining vecs to merge
    if merge_buffer_result_index - first_index_to_process == 2 {
        let first = &merge_buffer[first_index_to_process];
        let second = &merge_buffer[first_index_to_process + 1];

        merge(0, first, 0, second, &mut labels[node_id as usize]);
    }
    else if merge_buffer_result_index - first_index_to_process == 1 { // assign only remaing result
        let result = &mut merge_buffer[first_index_to_process];

        labels[node_id as usize].append(result);
    }
}