use crate::types::*;
use std::ops::Range;

pub trait Graph {

    fn add_edge(&mut self, start: NodeId, end: NodeId, weight: Weight);
    fn node_ids(&self) -> Range<NodeId>;
    fn edge_ids(&self, node: NodeId) -> Range<EdgeId>;
    fn target_id(&self, node: NodeId, edge_id: EdgeId) -> NodeId;
    fn weight(&self, node: NodeId, edge_id: EdgeId) -> Weight;
    
    fn print(&self);
}

pub struct GraphArray {

    first_edge: Vec<EdgeId>,
    target_node: Vec<NodeId>,
    weights: Vec<Weight>,
}

impl GraphArray {

    pub fn new(first_edge: EdgeIds, target_node: NodeIds, weights: Weights) -> GraphArray {
        GraphArray {
            first_edge,
            target_node,
            weights,
        }
    }
}

pub struct GraphList {

    edges: Vec<Vec<NodeId>>,
    weights: Vec<Vec<Weight>>,
}

impl GraphList {

    pub fn new(edges: Vec<Vec<NodeId>>, weights: Vec<Vec<Weight>>) -> GraphList {
        GraphList {
            edges,
            weights,
        }
    }

    pub fn from_array(first_edge: &Vec<EdgeId>, target_node: &Vec<NodeId>, weights: &Vec<Weight>) -> GraphList {
        let mut list_edges: Vec<Vec<NodeId>> = Vec::new();
        let mut list_weights: Vec<Vec<Weight>> = Vec::new();

        for node_index in 0..(first_edge.len() - 1) {
            let mut edges: Vec<NodeId> = Vec::new();
            let mut node_weights: Vec<Weight> = Vec::new();

            for edge_id in first_edge[node_index]..first_edge[node_index + 1] {
                edges.push(target_node[edge_id as usize]);
                node_weights.push(weights[edge_id as usize]);
            }

            list_edges.push(edges);
            list_weights.push(node_weights);
        } 

        GraphList {
            edges: list_edges,
            weights: list_weights,
        }
    }
}

impl Graph for GraphArray {

    fn add_edge(&mut self, start: NodeId, end: NodeId, weight: Weight) {
        let current_first_edge: EdgeId = self.first_edge[start as usize];

        // insert new target node and weight
        self.target_node.insert(current_first_edge as usize, end);
        self.weights.insert(current_first_edge as usize, weight);

        // increase all subsequent first edge ids
        for node_id in (start as usize + 1)..self.first_edge.len() {
            self.first_edge[node_id] += 1;
        }
    }

    fn node_ids(&self) -> Range<NodeId> {
        0..((self.first_edge.len() as NodeId) - 1)
    }

    fn edge_ids(&self, node: NodeId) -> Range<EdgeId> {
        self.first_edge[node as usize]..self.first_edge[(node as usize) + 1]
    }

    fn target_id(&self, _node: NodeId, edge_id: EdgeId) -> NodeId {
        self.target_node[edge_id as usize]
    }

    fn weight(&self, _node: NodeId, edge_id: EdgeId) -> Weight {
        self.weights[edge_id as usize]
    }

    fn print(&self) {
        for node_index in 0..(self.first_edge.len() - 1) {
            for edge_id in self.first_edge[node_index]..self.first_edge[node_index + 1] {
                println!("edge from {} to {} with weight: {}", node_index, self.target_node[edge_id as usize], self.weights[edge_id as usize]);
            }
        }
    }
}

impl Graph for GraphList {

    fn add_edge(&mut self, start: NodeId, end: NodeId, weight: Weight) {
        self.edges[start as usize].push(end);
        self.weights[start as usize].push(weight);
    }

    fn node_ids(&self) -> Range<NodeId> {
        0..(self.edges.len() as NodeId)
    }

    fn edge_ids(&self, node: NodeId) -> Range<EdgeId> {
        0..(self.edges[node as usize].len() as EdgeId)
    }

    fn target_id(&self, node: NodeId, edge_id: EdgeId) -> NodeId {
        self.edges[node as usize][edge_id as usize]
    }

    fn weight(&self, node: NodeId, edge_id: EdgeId) -> Weight {
        self.weights[node as usize][edge_id as usize]
    }

    fn print(&self) {
        for start_index in 0..self.edges.len() {
            for target_index in 0..self.edges[start_index].len() {
                println!("edge from {} to {} with weight: {}", start_index, self.edges[start_index][target_index], self.weights[start_index][target_index]);
            }
        }
    }
}