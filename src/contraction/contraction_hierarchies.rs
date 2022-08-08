use std::cmp::max;

use crate::{graph_representation::Graph, utils::binary_heap::{HeapElement, MinBinaryHeap}, types::{NodeId, Weight, Ranks, EdgeIds, NodeIds, Weights}};

use super::witness_search::WitnessSearch;

pub struct ContractionHierarchiesGraph {
    pub remaining_in_edges: Vec<Vec<(NodeId, Weight)>>,
    pub remaining_out_edges: Vec<Vec<(NodeId, Weight)>>,
    node_level: Vec<usize>,
    node_rank: Ranks,

    result_in_edges: Vec<Vec<(NodeId, Weight)>>,
    result_out_edges: Vec<Vec<(NodeId, Weight)>>
}

impl ContractionHierarchiesGraph {

    fn new(graph: &impl Graph) -> Self {
        let num_nodes = graph.node_ids().len();

        let mut remaining_in_edges: Vec<Vec<(NodeId, Weight)>> = vec![Vec::new(); num_nodes];
        let mut remaining_out_edges: Vec<Vec<(NodeId, Weight)>> = vec![Vec::new(); num_nodes];

        let node_level: Vec<usize> = vec![0; num_nodes];
        let node_rank: Ranks = vec![0; num_nodes];

        let mut result_in_edges: Vec<Vec<(NodeId, Weight)>> = vec![Vec::new(); num_nodes];
        let mut result_out_edges: Vec<Vec<(NodeId, Weight)>> = vec![Vec::new(); num_nodes];
        
        for source_node in graph.node_ids() {
            for edge_id in graph.edge_ids(source_node) {
                let target_node = graph.target_id(source_node, edge_id);
                let edge_weight = graph.weight(source_node, edge_id);

                remaining_in_edges[target_node as usize].push((source_node, edge_weight));
                remaining_out_edges[source_node as usize].push((target_node, edge_weight));

                result_in_edges[target_node as usize].push((source_node, edge_weight));
                result_out_edges[source_node as usize].push((target_node, edge_weight));
            }
        }

        Self {  
            remaining_in_edges,
            remaining_out_edges,
            node_level,
            node_rank,
            result_in_edges,
            result_out_edges
            
        }
    }

    fn export_result(&self) -> (EdgeIds, NodeIds, Weights, Ranks) {
        let num_nodes = self.result_out_edges.len();

        let mut first_edge: EdgeIds = vec![0; num_nodes + 1];
        let mut target_node: NodeIds = Vec::new();
        let mut weight: Weights = Vec::new();

        let mut first_edge_index = 0;

        for node_index in 0..num_nodes {
            first_edge[node_index] = first_edge_index;

            for (adj_node, arc_weight) in &self.result_out_edges[node_index] {
                target_node.push(*adj_node);
                weight.push(*arc_weight);

                first_edge_index += 1;
            }
        }

        first_edge[num_nodes] = first_edge_index;

        (first_edge, target_node, weight, self.node_rank.to_vec())
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, Ord, PartialOrd)]
struct PriorityQueueNode {
    score: usize,
    node_id: NodeId
}

impl HeapElement for PriorityQueueNode {
    
    fn unique_index(&self) -> usize {
        self.node_id as usize
    }
}

impl PriorityQueueNode {

    fn new(node_id: NodeId) -> Self {
        PriorityQueueNode { 
            score: 0, 
            node_id: node_id
        }
    }

    fn update_priority(&mut self, num_required_shortcuts: usize, remaining_graph: &ContractionHierarchiesGraph) {
        let in_edge_count = remaining_graph.remaining_in_edges[self.node_id as usize].len();
        let out_edge_count = remaining_graph.remaining_out_edges[self.node_id as usize].len();

        let level = remaining_graph.node_level[self.node_id as usize];
        let edge_count = in_edge_count + out_edge_count;

        self.score = (((0.05 * level as f32) + (1.0 * num_required_shortcuts as f32 + 1.0) / (edge_count as f32 + 1.0)) * 1000.0) as usize;
    }
}

pub fn contract_graph(graph: &impl Graph) -> (EdgeIds, NodeIds, Weights, Ranks) {
    let num_nodes = graph.node_ids().len();

    let mut remaining_graph = ContractionHierarchiesGraph::new(graph);
    let mut priority_queue: MinBinaryHeap<PriorityQueueNode> = MinBinaryHeap::new(graph.node_ids().len());
    let mut witness_search = WitnessSearch::new(num_nodes);

    // initially add all nodes to the priority queue
    for node_id in graph.node_ids() {
        let mut queue_node = PriorityQueueNode::new(node_id);
        let required_shortcuts = get_required_shortcuts(node_id, &remaining_graph, &mut witness_search);
        queue_node.update_priority(required_shortcuts.len(), &remaining_graph);

        priority_queue.insert(queue_node);
    }

    let mut num_contracted_nodes = 0;

    while let Some(mut entry) = priority_queue.pop() {
        let shortcuts = get_required_shortcuts(entry.node_id, &mut remaining_graph, &mut witness_search);

        entry.update_priority(shortcuts.len(), &remaining_graph);

        if !priority_queue.is_empty() && priority_queue.get_min().unwrap().score < entry.score {
            priority_queue.insert(entry);
            continue;
        }

        if num_contracted_nodes % 1000 == 0 {
            println!("finished nodes: {} remaining nodes: {}, progress: {}%", num_contracted_nodes, priority_queue.len(), (num_contracted_nodes as f64 / num_nodes as f64) * 100.0);
        }

        contract_node(entry.node_id, shortcuts, &mut remaining_graph);
        remaining_graph.node_rank[entry.node_id as usize] = num_contracted_nodes;
        num_contracted_nodes += 1;
    }

    remaining_graph.export_result()
}

fn add_shortcut_edge(graph: &mut ContractionHierarchiesGraph, source_node: NodeId, target_node: NodeId, shortcut_weight: Weight) {
    graph.remaining_out_edges[source_node as usize].push((target_node, shortcut_weight));
    graph.remaining_in_edges[target_node as usize].push((source_node, shortcut_weight));

    graph.result_out_edges[source_node as usize].push((target_node, shortcut_weight));
    graph.result_in_edges[target_node as usize].push((source_node, shortcut_weight));

    remove_duplicated_edges(&mut graph.remaining_out_edges[source_node as usize]);
    remove_duplicated_edges(&mut graph.remaining_in_edges[target_node as usize]);
    remove_duplicated_edges(&mut graph.result_out_edges[source_node as usize]);
    remove_duplicated_edges(&mut graph.result_in_edges[target_node as usize]);
}

fn remove_duplicated_edges(arcs: &mut Vec<(NodeId, Weight)>) {
    arcs.sort_by_key(|(target_node, weight)| (*target_node, *weight));
    arcs.dedup_by_key(|(target_node, _weight)| *target_node);
}

fn contract_node(node_id: NodeId, shortcuts_to_add: Vec<(NodeId, NodeId, Weight)>, remaining_graph: &mut ContractionHierarchiesGraph) {
    // add the shortcut edges
    for (source_node, target_node, weight) in shortcuts_to_add {
        add_shortcut_edge(remaining_graph, source_node, target_node, weight);
    }

    let current_node_level = remaining_graph.node_level[node_id as usize];

    // remove edges that are now longer needed
    for (source_node, _weight) in &remaining_graph.remaining_in_edges[node_id as usize] {
        remaining_graph.remaining_out_edges[*source_node as usize].retain(|(remove_target, _weight)| *remove_target != node_id);
        remaining_graph.node_level[*source_node as usize] = max(remaining_graph.node_level[*source_node as usize], current_node_level + 1);
    }

    for (target_node, _weight) in &remaining_graph.remaining_out_edges[node_id as usize] {
        remaining_graph.remaining_in_edges[*target_node as usize].retain(|(remove_source, _weight)| *remove_source != node_id);
        remaining_graph.node_level[*target_node as usize] = max(remaining_graph.node_level[*target_node as usize], current_node_level + 1);
    }

    remaining_graph.remaining_in_edges[node_id as usize].clear();
    remaining_graph.remaining_out_edges[node_id as usize].clear();
}

fn get_required_shortcuts(node_to_contract: NodeId, remaining_graph: &ContractionHierarchiesGraph, witness_search: &mut WitnessSearch) -> Vec<(NodeId, NodeId, Weight)> {
    let mut shortcuts = Vec::new();
    
    for (source_node, source_weight) in &remaining_graph.remaining_in_edges[node_to_contract as usize] {
        witness_search.prepare(*source_node, node_to_contract);
        
        for (target_node, target_weight) in &remaining_graph.remaining_out_edges[node_to_contract as usize] {
            let shortcut_weight = *source_weight + *target_weight;
            let witness_search_distance = witness_search.continue_search(*target_node, shortcut_weight, remaining_graph);

            if shortcut_weight < witness_search_distance {
                shortcuts.push((*source_node, *target_node, *source_weight + *target_weight));
            }
        }
    }

    shortcuts
}