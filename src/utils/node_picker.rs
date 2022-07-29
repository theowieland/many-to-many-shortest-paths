use crate::types::*;
use crate::graph_algorithms::*;
use crate::utils::binary_heap::MinBinaryHeap;
use crate::utils::io::export_node_ids;
use crate::utils::io::read_node_ids;
use rand::thread_rng;
use rand::seq::SliceRandom;
use rand::Rng;
use std::path::Path;
use crate::data_structures::{ArrayStructure, ValidFlags};

/// used to store reusable objects used by the node picker
pub struct NodePicker {
    to_visit: MinBinaryHeap<DijkstraState>,
    distances: ValidFlags<Weight>,

    first_edge: Vec<EdgeId>,
    target_node: Vec<NodeId>,
    weights: Vec<Weight>
}

/// loads node ids from disk
pub struct CachedNodePicker<'a> {
    path: &'a Path
}

impl<'a> CachedNodePicker<'a> {

    pub fn new(path: &'a Path) -> Self {
        CachedNodePicker {
            path
        }
    }

    pub fn get_sources(&self, iteration: usize, ball_size: usize, source_target_size: &str) -> Vec<NodeId> {
        let file_name = format!("{}/{}/{}/iteration_{}_ball_size_{}", self.path.display(), source_target_size, "sources", iteration, ball_size);

        let nodes: NodeIds = read_node_ids(&Path::new(&file_name)).unwrap();

        nodes
    }

    pub fn get_targets(&self, iteration: usize, ball_size: usize, source_target_size: &str) -> Vec<NodeId> {
        let file_name = format!("{}/{}/{}/iteration_{}_ball_size_{}", self.path.display(), source_target_size, "targets", iteration, ball_size);
        
        let nodes: NodeIds = read_node_ids(&Path::new(&file_name)).unwrap();

        nodes
    }
}

/// the node picker can be used to pick random nodes from a given graph
impl NodePicker {

    pub fn new(num_vertices: usize, first_edge: Vec<EdgeId>, target_node: Vec<NodeId>, weights: Vec<Weight>) -> Self {
        NodePicker {
            to_visit: MinBinaryHeap::new(num_vertices),
            distances: ValidFlags::new(num_vertices, INFINITY),

            first_edge,
            target_node,
            weights,
        }
    }

    /// picks a random start node to run dijkstra's algorithm from until a given amount of nodes has been settled
    /// returns all settled nodes
    pub fn grow_random_dijkstras_ball(&mut self, center: NodeId, ball_size: usize) -> Vec<NodeId> {
        let mut settled: Vec<NodeId> = Vec::new();

        self.to_visit.clear();
        self.to_visit.insert(DijkstraState {distance: 0, node_id: center});

        self.distances.reset();
        self.distances.set(center as usize, 0);

        while let Some(DijkstraState { distance: current_distance, node_id: current_node }) = self.to_visit.pop() {
            settled.push(current_node);

            // abort if enough nodes have been visited or only one sample is needed and we have already seen more than one node
            if settled.len() >= ball_size { 
                self.to_visit.clear();
                break;
            }

            // relax edges
            let start_index = self.first_edge[current_node as usize] as usize;
            let end_index = self.first_edge[current_node as usize + 1] as usize;
            for (target_node, weight) in self.target_node[start_index..end_index].iter().zip(self.weights[start_index..end_index].iter()) {
                let new_distance = current_distance + *weight;
    
                if self.distances[*target_node as usize] > new_distance {
                    self.distances.set(*target_node as usize, new_distance);
    
                    self.to_visit.insert_or_decrease(DijkstraState {distance: new_distance, node_id: *target_node});
                }
            }
        }

        settled
    }

    /// runs the dijkstra algorithm starting from the given center node until ball_size nodes are settled
    /// then picks sample_size, unique, random nodes from the given settled nodes
    /// the ball size has to be larger or equal to the sample size
    pub fn get_random_nodes(&mut self, center: NodeId, ball_size: usize, sample_size: usize) -> Vec<NodeId> {
        assert!(ball_size >= sample_size);

        let mut visited: Vec<NodeId> = self.grow_random_dijkstras_ball(center, ball_size);
    
        //shuffle the result and extract sample_size nodes
        visited.shuffle(&mut thread_rng());
        let result_nodes = visited[0..sample_size].to_vec();

        result_nodes
    }

    /// runs the dijkstra algorithm starting from the given center node until ball_size nodes are settled
    /// then picks sample_size, unique, random nodes from the given settled nodes
    /// the ball size has to be larger or equal to the sample size
    pub fn get_two_random_sets(&mut self, center: NodeId, ball_size: usize, sources: usize, targets: usize) -> (Vec<NodeId>, Vec<NodeId>) {
        assert!(ball_size >= sources && ball_size >= targets);

        let mut visited: Vec<NodeId> = self.grow_random_dijkstras_ball(center, ball_size);
    
        //shuffle the result and extract source nodes
        visited.shuffle(&mut thread_rng());
        let result_sources = visited[0..sources].to_vec();

        //shuffle the result and extract target nodes
        visited.shuffle(&mut thread_rng());
        let result_targets = visited[0..targets].to_vec();

        (result_sources, result_targets)
    }

    /// runs the dijkstra algorithm starting from a random node until ball_size nodes are settled
    /// then picks sample_size, unique, random nodes from the given settled nodes as well as one additonal random node from the discovered ball
    /// the ball size has to be larger or equal to the sample size
    pub fn get_random_nodes_and_another_single_random_node(&mut self, ball_size: usize, sample_size: usize) -> (NodeId, Vec<NodeId>) {
        assert!(ball_size >= sample_size);

        let mut settled = self.grow_random_dijkstras_ball(self.get_random_node(), ball_size);
    
        // shuffle the result and extract sample_size nodes
        let mut rng = rand::thread_rng();
        settled.shuffle(&mut rng);
        let result_nodes = settled[0..sample_size].to_vec();
        let random_single_node = settled[rng.gen_range(0..settled.len()) as usize] as NodeId;

        (random_single_node, result_nodes)
    }

    /// runs the dijkstra algorithm starting from a random node until ball_size nodes are settled
    /// then picks sample_size, unique, random nodes from the given settled nodes as well as one additonal random node
    /// the ball size has to be larger or equal to the sample size
    pub fn get_random_nodes_and_separate_single_random_node(&mut self, ball_size: usize, sample_size: usize) -> (NodeId, Vec<NodeId>) {
        assert!(ball_size >= sample_size);

        let mut settled = self.grow_random_dijkstras_ball(self.get_random_node(), ball_size);
    
        // shuffle the result and extract sample_size nodes
        let mut rng = rand::thread_rng();
        settled.shuffle(&mut rng);
        let result_nodes = settled[0..sample_size].to_vec();
        let random_single_node = self.get_random_node();

        (random_single_node, result_nodes)
    }

    /// picks a single random node from the given graph
    pub fn get_random_node(&self) -> NodeId {
        let mut rng = rand::thread_rng();
        
        rng.gen_range(0..(self.first_edge.len() -1)) as NodeId
    }

    pub fn create_many_to_many_sources_equals_targets(&mut self, iterations: usize, path: &Path) {
        let combinations: Vec<(String, usize, Vec<usize>)> = vec![
            (String::from("many_to_many_2^10"), 2_i32.pow(10) as usize, (10..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()),
            (String::from("many_to_many_2^12"), 2_i32.pow(12) as usize, (12..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()),
            (String::from("many_to_many_2^14"), 2_i32.pow(14) as usize, (14..25).map(|exponent| 2_i32.pow(exponent) as usize).collect())
        ];

        for (source_target_size, sample_size, ball_sizes) in combinations {
            for ball_size in ball_sizes {
                for iteration in 0..iterations {
                    let node_ids = self.get_random_nodes(self.get_random_node(), ball_size, sample_size);
                    let file_name_sources = format!("{}/sources/iteration_{}_ball_size_{}", source_target_size, iteration, ball_size);
                    NodePicker::save_file(&node_ids, path, file_name_sources);

                    let file_name_targets = format!("{}/targets/iteration_{}_ball_size_{}", source_target_size, iteration, ball_size);
                    NodePicker::save_file(&node_ids, path, file_name_targets);
                }
            }
        }
    }

    pub fn create_many_to_many_same_ball(&mut self, iterations: usize, path: &Path) {
        let combinations: Vec<(String, usize, Vec<usize>)> = vec![
            (String::from("many_to_many_2^10"), 2_i32.pow(10) as usize, (10..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()),
            (String::from("many_to_many_2^12"), 2_i32.pow(12) as usize, (12..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()),
            (String::from("many_to_many_2^14"), 2_i32.pow(14) as usize, (14..25).map(|exponent| 2_i32.pow(exponent) as usize).collect())
        ];

        for (source_target_size, sample_size, ball_sizes) in combinations {
            for ball_size in ball_sizes {
                for iteration in 0..iterations {
                    let (sources, targets) = self.get_two_random_sets(self.get_random_node(), ball_size, sample_size, sample_size);
                    let file_name_sources = format!("{}/sources/iteration_{}_ball_size_{}", source_target_size, iteration, ball_size);
                    NodePicker::save_file(&sources, path, file_name_sources);

                    let file_name_targets = format!("{}/targets/iteration_{}_ball_size_{}", source_target_size, iteration, ball_size);
                    NodePicker::save_file(&targets, path, file_name_targets);
                }
            }
        }
    }

    pub fn create_many_to_many_separate_balls(&mut self, iterations: usize, path: &Path) {
        let combinations: Vec<(String, (usize, Vec<usize>), (usize, Vec<usize>))> = vec![
            (
                String::from("many_to_many_2^10"),
                (2_i32.pow(10) as usize, (10..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()), // sources
                (2_i32.pow(10) as usize, (10..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()) // targets
            ),
            (
                String::from("many_to_many_2^12"),
                (2_i32.pow(12) as usize, (12..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()), // sources
                (2_i32.pow(12) as usize, (12..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()) // targets
            ),
            (
                String::from("many_to_many_2^14"),
                (2_i32.pow(14) as usize, (14..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()), // sources
                (2_i32.pow(14) as usize, (14..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()) // targets
            ),
            (
                String::from("many_to_many_2^10_to_2^14"),
                (2_i32.pow(10) as usize, (10..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()), // sources
                (2_i32.pow(14) as usize, (14..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()) // targets
            ),
            (
                String::from("many_to_many_2^14_to_2^10"),
                (2_i32.pow(14) as usize, (14..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()), // sources
                (2_i32.pow(10) as usize, (10..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()) // targets
            )
        ];

        for (source_target_size, (source_size, source_balls), (target_size, target_balls)) in combinations {
            
            // create source files
            for source_ball_size in source_balls {
                for iteration in 0..iterations {
                    let node_ids = self.get_random_nodes(self.get_random_node(), source_ball_size, source_size);
                    let file_name = format!("{}/sources/iteration_{}_ball_size_{}", source_target_size, iteration, source_ball_size);

                    NodePicker::save_file(&node_ids, path, file_name);
                }
            }

            // create target files
            for target_ball_size in target_balls {
                for iteration in 0..iterations {
                    let node_ids = self.get_random_nodes(self.get_random_node(), target_ball_size, target_size);
                    let file_name = format!("{}/targets/iteration_{}_ball_size_{}", source_target_size, iteration, target_ball_size);

                    NodePicker::save_file(&node_ids, path, file_name);
                }
            }
        }
    }

    pub fn create_one_to_many_same_balls(&mut self, iterations: usize, path: &Path) {
        let combinations: Vec<(String, String, usize, Vec<usize>)> = vec![
            (String::from("one_to_many_2^10"), String::from("many_to_one_2^10"), 2_i32.pow(10) as usize, (10..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()),
            (String::from("one_to_many_2^12"), String::from("many_to_one_2^12"), 2_i32.pow(12) as usize, (12..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()),
            (String::from("one_to_many_2^14"), String::from("many_to_one_2^14"), 2_i32.pow(14) as usize, (14..25).map(|exponent| 2_i32.pow(exponent) as usize).collect())
        ];

        for (source_target_size_one_to_many, source_target_size_many_to_one, sample_size, ball_sizes) in combinations {
            
            // generate one to many data
            for ball_size in &ball_sizes {
                for iteration in 0..iterations {
                    let (single_random, random_ids) = self.get_random_nodes_and_another_single_random_node(*ball_size as usize, sample_size);
    
                    let file_name_sources_one_to_many = format!("{}/sources/iteration_{}_ball_size_{}", source_target_size_one_to_many, iteration, ball_size);
                    let file_name_targets_one_to_many = format!("{}/targets/iteration_{}_ball_size_{}", source_target_size_one_to_many, iteration, ball_size);

                    NodePicker::save_file(&vec![single_random], path, file_name_sources_one_to_many);
                    NodePicker::save_file(&random_ids, path, file_name_targets_one_to_many);

                    let file_name_sources_many_to_one = format!("{}/sources/iteration_{}_ball_size_{}", source_target_size_many_to_one, iteration, ball_size);
                    let file_name_targets_many_to_one = format!("{}/targets/iteration_{}_ball_size_{}", source_target_size_many_to_one, iteration, ball_size);

                    NodePicker::save_file(&vec![single_random], path, file_name_targets_many_to_one);
                    NodePicker::save_file(&random_ids, path, file_name_sources_many_to_one);
                }
            }
        }
    }

    pub fn create_one_to_many_separate_balls(&mut self, iterations: usize, path: &Path) {
        let one_to_many_combinations: Vec<(String, String, usize, Vec<usize>)> = vec![
            (String::from("one_to_many_2^10"), String::from("many_to_one_2^10"), 2_i32.pow(10) as usize, (10..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()),
            (String::from("one_to_many_2^12"), String::from("many_to_one_2^12"), 2_i32.pow(12) as usize, (12..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()),
            (String::from("one_to_many_2^14"), String::from("many_to_one_2^14"), 2_i32.pow(14) as usize, (14..25).map(|exponent| 2_i32.pow(exponent) as usize).collect())
        ];

        for (source_target_size_one_to_many, source_target_size_many_to_one, sample_size, ball_sizes) in one_to_many_combinations {
            
            // generate one to many data
            for ball_size in &ball_sizes {
                for iteration in 0..iterations {
                    let (single_random, random_ids) = self.get_random_nodes_and_separate_single_random_node(*ball_size as usize, sample_size);
    
                    let file_name_sources_one_to_many = format!("{}/sources/iteration_{}_ball_size_{}", source_target_size_one_to_many, iteration, ball_size);
                    let file_name_targets_one_to_many = format!("{}/targets/iteration_{}_ball_size_{}", source_target_size_one_to_many, iteration, ball_size);

                    NodePicker::save_file(&vec![single_random], path, file_name_sources_one_to_many);
                    NodePicker::save_file(&random_ids, path, file_name_targets_one_to_many);

                    let file_name_sources_many_to_one = format!("{}/sources/iteration_{}_ball_size_{}", source_target_size_many_to_one, iteration, ball_size);
                    let file_name_targets_many_to_one = format!("{}/targets/iteration_{}_ball_size_{}", source_target_size_many_to_one, iteration, ball_size);

                    NodePicker::save_file(&vec![single_random], path, file_name_targets_many_to_one);
                    NodePicker::save_file(&random_ids, path, file_name_sources_many_to_one);
                }
            }
        }
    }

    fn save_file(node_ids: &Vec<NodeId>, path: &Path, file_name: String) {
        match export_node_ids(&path.join(file_name.clone()), node_ids) {
            Ok(_) => println!("successfully saved: {:?}", &path.join(file_name)),
            Err(e) => println!("error saving: {:?} {}", &path.join(file_name), e)
        }
    }
} 