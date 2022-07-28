use std::path::Path;

use crate::{types::{EdgeId, NodeId, Weight, EdgeIds, NodeIds, Weights, Ranks}, many_to_many::{many_to_many_algorithm::ManyToManyAlgorithm, bucket::BucketManyToMany, advanced_bucket::AdvancedBucketManyToMany, rphast::RPHAST, rphast_multiple_trees_fwd_rank::RPHASTMultipleTreesFwdRank, rphast_multiple_trees_fwd_dijkstra::RPHASTMultipleTreesFwdDijkstra, rphast_multiple_trees_fwd_buckets::RPHASTMultipleTreesFwdBuckets, rphast_simd::SIMDRPHAST, rphast_batched_simd::BatchedSIMDRPHAST, hub_labels_partially_pruned::HubLabelsPartiallyPruned, hub_labels_top_down::HubLabelsTopDown, lazy_simd::LazySIMD, lazy_batched_buckets::LazyBatchedBuckets, lazy::LazyRPHAST, dijkstra::DijkstraManyToMany, advanced_bucket_batched::AdvancedBucketBatchedManyToMany, forward_backward_buckets_merged::ForwardBackwardBucketsMergedManyToMany, lazy_one_to_many::LazyRPHASTOneToMany}, utils::{combine_restricted_graph, io::read_graph_data, create_restricted_graph, measure_time, node_picker::CachedNodePicker}};

pub struct Experiment {
    pub name: String, // a unique identifier for this experiment
    pub function: Box<dyn Fn(&mut ExperimentData, &SourceTargetConfiguration, CustomExperimentVariables) -> ()> // method that executes the experiment
}

impl Experiment {
    
    pub fn new(name: String, function: impl Fn(&mut ExperimentData, &SourceTargetConfiguration, CustomExperimentVariables) -> () + 'static) -> Self {
        Experiment {
            name,
            function: Box::new(function)
        }
    }

    pub fn all_experiments() -> Vec<Self> {
        let mut experiments = Vec::new();

        experiments.push(Experiment::new(
            String::from("buckets_no_pruning_bucket_size"), 
            |data, source_target_configuration, variables| {
                let mut algorithm = BucketManyToMany::new(&data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight, &data.ranks);
        
                for ball_size in &source_target_configuration.ball_sizes {
                    for iteration in 0..variables.num_iterations {
                        let (sources, targets) = get_sources_and_targets(&mut data.node_picker, iteration, *ball_size, &source_target_configuration.folder);
                        algorithm.initialize(&sources, &targets);
                        algorithm.select_targets_no_stall_on_demand(&targets);
                        
                        let (total_entries, visited_buckets) = algorithm.get_average_number_of_bucket_entries();
                        println!("{},{},{}", ball_size, total_entries, visited_buckets);
                    }
                }
            }
        ));

        experiments.push(Experiment::new(
            String::from("simultaneous_buckets_bucket_size"), 
            |data, source_target_configuration, variables| {
                let mut algorithm = AdvancedBucketManyToMany::new(&data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight, &data.ranks);
        
                for ball_size in &source_target_configuration.ball_sizes {
                    for iteration in 0..variables.num_iterations {
                        let (sources, targets) = get_sources_and_targets(&mut data.node_picker, iteration, *ball_size, &source_target_configuration.folder);
                        algorithm.initialize(&sources, &targets);
                        algorithm.select(&sources, &targets);

                        let (total_entries, visited_buckets) = algorithm.get_average_number_of_bucket_entries();
                        println!("{},{},{}", ball_size, total_entries, visited_buckets);
                    }
                }
            }
        ));

        experiments.push(Experiment::new(
            String::from("buckets_bucket_size"), 
            |data, source_target_configuration, variables| {
                let mut algorithm = BucketManyToMany::new(&data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight, &data.ranks);
        
                for ball_size in &source_target_configuration.ball_sizes {
                    for iteration in 0..variables.num_iterations {
                        let (sources, targets) = get_sources_and_targets(&mut data.node_picker, iteration, *ball_size, &source_target_configuration.folder);
                        algorithm.initialize(&sources, &targets);
                        algorithm.select(&sources, &targets);

                        let (total_entries, visited_buckets) = algorithm.get_average_number_of_bucket_entries();
                        println!("{},{},{}", ball_size, total_entries, visited_buckets);
                    }
                }
            }   
        ));

        experiments.push(Experiment::new(
            String::from("rphast_bwd_dfs"),
            |data, source_target_configuration, variables| {
                let mut algorithm = RPHAST::new(&data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight, &data.ranks);
                
                for ball_size in &source_target_configuration.ball_sizes {
                    for iteration in 0..variables.num_iterations {
                        let (sources, targets) = get_sources_and_targets(&mut data.node_picker, iteration, *ball_size, &source_target_configuration.folder);
                        algorithm.initialize(&sources, &targets);
                        algorithm.select(&sources, &targets);
                        
                        let query_time = measure_time(|| {
                            algorithm.query(&sources, &targets);
                        });
                        
                        println!("{},{:?}", ball_size, query_time.as_nanos());
                    }
                }
            }
        ));

        experiments.push(Experiment::new(
            String::from("rphast_bwd_ranks"),
            |data, source_target_configuration, variables| {
                let mut algorithm = RPHAST::new(&data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight, &data.ranks);
                
                for ball_size in &source_target_configuration.ball_sizes {
                    for iteration in 0..variables.num_iterations {
                        let (sources, targets) = get_sources_and_targets(&mut data.node_picker, iteration, *ball_size, &source_target_configuration.folder);
                        algorithm.initialize(&sources, &targets);
                        algorithm.select_targets_decreasing_rank(&targets);
                        
                        let query_time = measure_time(|| {
                            algorithm.query(&sources, &targets);
                        });
                        
                        println!("{},{:?}", ball_size, query_time.as_nanos());
                    }
                }
            }
        ));

        experiments.push(Experiment::new(
            String::from("sbi_fwd_bwd_mutual_buckets"),
            |data, source_target_configuration, variables| {
                let mut algorithm = ForwardBackwardBucketsMergedManyToMany::new(&data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight, &data.ranks);
                
                for ball_size in &source_target_configuration.ball_sizes {
                    for iteration in 0..variables.num_iterations {
                        let (sources, targets) = get_sources_and_targets(&mut data.node_picker, iteration, *ball_size, &source_target_configuration.folder);
                        let query_time = measure_time(|| {
                            algorithm.initialize(&sources, &targets);
                            algorithm.populate_buckets(&sources, &targets);
                            let visited = algorithm.calculate_mutual_buckets();
                            algorithm.merge_buckets(&visited);
                        });
                        
                        let visited = algorithm.calculate_mutual_buckets();
                        let (fwd_size, bwd_size) = algorithm.get_total_bucket_size(&visited);
                        
                        println!("{},{:?},{},{},{}", ball_size, query_time.as_nanos(), visited.len(), fwd_size, bwd_size);
                    }
                }
            }
        ));

        experiments.push(Experiment::new(
            String::from("rphast_restricted_graph_size"),
            |data, source_target_configuration, variables| {
                let mut algorithm = RPHAST::new(&data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight, &data.ranks);
                
                for ball_size in &source_target_configuration.ball_sizes {
                    for iteration in 0..variables.num_iterations {
                        let (sources, targets) = get_sources_and_targets(&mut data.node_picker, iteration, *ball_size, &source_target_configuration.folder);
                        algorithm.initialize(&sources, &targets);
                        algorithm.select(&sources, &targets);
                        
                        let (restricted_vertices, restricted_edges) = algorithm.get_restricted_graph_size();
                        println!("{},{},{}", ball_size, restricted_vertices, restricted_edges);
                    }
                }
            }
        ));

        experiments
    }
}

pub struct ExperimentData<'a> {
    fwd_first_edge: Vec<EdgeId>,
    fwd_head: Vec<NodeId>,
    fwd_weight: Vec<Weight>,
    bwd_first_edge: Vec<EdgeId>,
    bwd_head: Vec<NodeId>,
    bwd_weight: Vec<Weight>,
    ranks: Vec<usize>,
    pub node_picker: CachedNodePicker<'a>
}

impl<'a> ExperimentData<'a> {

    pub fn new(graph_path: &'a Path, data_path: &'a Path) -> Self {
        let (fwd_first_edge, fwd_head, fwd_weight, bwd_first_edge, bwd_head, bwd_weight, ranks) = read_graph(graph_path);

        let node_picker: CachedNodePicker = CachedNodePicker::new(data_path);

        ExperimentData {
            fwd_first_edge,
            fwd_head,
            fwd_weight,
            bwd_first_edge,
            bwd_head,
            bwd_weight,
            ranks,
            node_picker
        }
    }
}

pub struct SourceTargetConfiguration {
    pub folder: String,
    pub source_size: usize,
    pub target_size: usize,
    pub ball_sizes: Vec<usize>
}

impl SourceTargetConfiguration {

    pub fn new(folder: String, source_size: usize, target_size: usize, ball_sizes: Vec<usize>) -> Self {
        SourceTargetConfiguration {
            folder,
            source_size,
            target_size,
            ball_sizes
        }
    }

    /// creates a list of all source target configurations used for these experiments
    pub fn all_source_target_configuration() -> Vec<Self> {
        let mut configurations: Vec<Self> = Vec::new();
        
        configurations.push(Self::new(
            String::from("many_to_many_2^10"), 
            2_i32.pow(10) as usize, 
            2_i32.pow(10) as usize,
            (10..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()
        ));

        configurations.push(Self::new(
            String::from("many_to_many_2^12"),
            2_i32.pow(12) as usize,
            2_i32.pow(12) as usize,
            (12..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()
        ));
        
        configurations.push(Self::new(
            String::from("many_to_many_2^14"),
            2_i32.pow(14) as usize,
            2_i32.pow(14) as usize, 
            (14..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()
        ));

        configurations.push(Self::new(
            String::from("many_to_many_2^10_to_2^14"),
            2_i32.pow(10) as usize,
            2_i32.pow(14) as usize, 
            (14..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()
        ));

        configurations.push(Self::new(
            String::from("many_to_many_2^14_to_2^10"),
            2_i32.pow(14) as usize, 
            2_i32.pow(10) as usize, 
            (14..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()
        ));

        configurations.push(Self::new(
            String::from("one_to_many_2^10"),
            1, 
            2_i32.pow(10) as usize, 
            (10..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()
        ));

        configurations.push(Self::new(
            String::from("one_to_many_2^12"),
            1, 
            2_i32.pow(12) as usize, 
            (12..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()
        ));

        configurations.push(Self::new(
            String::from("one_to_many_2^14"),
            1, 
            2_i32.pow(14) as usize, 
            (14..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()
        ));

        configurations.push(Self::new(
            String::from("many_to_one_2^10"),
            2_i32.pow(10) as usize, 
            1, 
            (10..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()
        ));

        configurations.push(Self::new(
            String::from("many_to_one_2^12"),
            2_i32.pow(12) as usize, 
            1, 
            (12..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()
        ));

        configurations.push(Self::new(
            String::from("many_to_one_2^14"),
            2_i32.pow(14) as usize, 
            1, 
            (14..25).map(|exponent| 2_i32.pow(exponent) as usize).collect()
        ));

        configurations
    }
}

// additional variables to modify for each experiment
pub struct ExperimentVariables {

    pub batch_size: usize
}

pub struct CustomExperimentVariables {
    pub num_iterations: usize,
    pub batch_size: usize
}

pub struct ManyToManyExperimentAlgorithm {
    pub name: String,
    pub algorithm: Box<dyn Fn(&mut ExperimentData, ExperimentVariables) -> Box<dyn ManyToManyAlgorithm>>  
}

impl ManyToManyExperimentAlgorithm {
    
    fn new(name: String, function: impl Fn(&mut ExperimentData, ExperimentVariables) -> Box<dyn ManyToManyAlgorithm + 'static> + 'static) -> Self {
        Self {
            name,
            algorithm: Box::new(function)
        }
    }

    /// creates a list of all many to many algorithms used for these experiments
    pub fn all_algorithms() -> Vec<Self> {
        let mut algorithms = Vec::new();

        algorithms.push(Self::new(
            String::from("forward_backward_sbi"),
            |data, _configuration| {
                Box::new(ForwardBackwardBucketsMergedManyToMany::new(&data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight, &data.ranks))
            }
        ));

        algorithms.push(Self::new(
            String::from("buckets"),
            |data, _configuration| {
                Box::new(BucketManyToMany::new(&data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight, &data.ranks))
            }
        ));
    
        algorithms.push(Self::new(
            String::from("simultaneous_buckets"),
            |data, _configuration| {
                Box::new(AdvancedBucketManyToMany::new(&data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight, &data.ranks))
            }
        ));

        algorithms.push(Self::new(
            String::from("simultaneous_buckets_batched"),
            |data, configuration| {
                Box::new(AdvancedBucketBatchedManyToMany::new(configuration.batch_size, &data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight, &data.ranks))
            }
        ));
    
        algorithms.push(Self::new(
            String::from("rphast"),
            |data, _configuration| {
                Box::new(RPHAST::new(&data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight, &data.ranks))
            }
        ));
    
        algorithms.push(Self::new(
            String::from("rphast_multiple_trees_fwd_rank"),
            |data, configuration| {
                Box::new(RPHASTMultipleTreesFwdRank::new(configuration.batch_size, &data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight, &data.ranks))
            }
        ));
    
        algorithms.push(Self::new(
            String::from("rphast_multiple_trees_fwd_dijkstra"),
            |data, configuration| {
                Box::new(RPHASTMultipleTreesFwdDijkstra::new(configuration.batch_size, &data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight))
            }
        ));
    
        algorithms.push(Self::new(
            String::from("rphast_multiple_trees_fwd_buckets"),
            |data, configuration| {
                Box::new(RPHASTMultipleTreesFwdBuckets::new(configuration.batch_size, &data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight, &data.ranks))
            }
        ));
    
        algorithms.push(Self::new(
            String::from("simd_rphast"),
            |data, _configuration| {
                Box::new(SIMDRPHAST::new(&data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight, &data.ranks))
            }
        ));
    
        algorithms.push(Self::new(
            String::from("simd_rphast_batched"),
            |data, configuration| {
                Box::new(BatchedSIMDRPHAST::new(configuration.batch_size, &data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight, &data.ranks))
            }
        ));
    
        algorithms.push(Self::new(
            String::from("hub_labels_partially_pruned"),
            |data, _configuration| {
                Box::new(HubLabelsPartiallyPruned::new(&data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight))
            }
        ));
    
        algorithms.push(Self::new(
            String::from("hub_labels_top_down"),
            |data, _configuration| {
                Box::new(HubLabelsTopDown::new(&data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight, &data.ranks))
            }
        ));
    
        algorithms.push(Self::new(
            String::from("lazy_simd"),
            |data, _configuration| {
                Box::new(LazySIMD::new(&data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight, &data.ranks))
            }
        ));
    
        algorithms.push(Self::new(
            String::from("lazy_batched_advanced_buckets"),
            |data, configuration| {
                Box::new(LazyBatchedBuckets::new(configuration.batch_size, &data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight, &data.ranks))
            }
        ));
    
        algorithms.push(Self::new(
            String::from("lazy"),
            |data, _configuration| {
                Box::new(LazyRPHAST::new(&data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight))
            }
        ));

        algorithms.push(Self::new(
            String::from("lazy_one_to_many"),
            |data, _configuration| {
                Box::new(LazyRPHASTOneToMany::new(&data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight))
            }
        ));
    
        algorithms.push(Self::new(
            String::from("dijkstra"),
            |data, _configuration| {
                let num_vertices = data.ranks.len();
                let (first_out, head, weight) = combine_restricted_graph(num_vertices, &data.fwd_first_edge, &data.fwd_head, &data.fwd_weight, &data.bwd_first_edge, &data.bwd_head, &data.bwd_weight);
                
                Box::new(DijkstraManyToMany::new(&first_out, &head, &weight))
            }        
        ));

        algorithms
    }
}

pub fn read_graph(path: &Path) -> (EdgeIds, NodeIds, Weights, EdgeIds, NodeIds, Weights, Ranks) {
    // read the graph to compute querys on
    let (first_edge, target_node, weight, rank) = read_graph_data(&path).unwrap();

    let (fwd_first_edge, fwd_target_node, fwd_weight, bwd_first_edge, bwd_source_node, bwd_weight) = create_restricted_graph(&first_edge, &target_node, &weight, &rank);
    (fwd_first_edge, fwd_target_node, fwd_weight, bwd_first_edge, bwd_source_node, bwd_weight, rank)
}

pub fn get_sources_and_targets(node_picker: &mut CachedNodePicker, iteration: usize, ball_size: usize, source_target_size: &str) -> (Vec<NodeId>, Vec<NodeId>) {
    let sources = node_picker.get_sources(iteration, ball_size, source_target_size);
    let targets = node_picker.get_targets(iteration, ball_size, source_target_size);
    
    (sources, targets)
}