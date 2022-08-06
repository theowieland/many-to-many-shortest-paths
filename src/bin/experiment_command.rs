use std::path::Path;
use clap::Parser;
use shortest_path_algorithms::{experiments::experiment::{SourceTargetConfiguration, ManyToManyExperimentAlgorithm, ExperimentData, Experiment, get_sources_and_targets, ExperimentVariables, read_graph, CustomExperimentVariables}, utils::{combine_restricted_graph, node_picker::NodePicker, measure_time}};

#[derive(Parser)]
struct Opts {

    #[clap(subcommand)]
    subcmd: SubCommand
}

#[derive(Parser)]
enum SubCommand {
    CustomExperiment(CustomExperimentCommand),
    Experiment(ExperimentCommand),
    CreateData(CreateDataCommand)
}

#[derive(Parser)]
struct ExperimentCommand {

    /// name of the algorithm to perfom the experiment with
    #[clap(short, long)]
    algorithm: String,

    /// path to the folder that stores used source and target nodes
    #[clap(short, long)]
    data_path: String,

    /// path to the folder that stores the used graph
    #[clap(short, long)]
    graph_path: String,

    /// source-target size configuration
    #[clap(short, long)]
    source_target_size: String,

    /// number of iterations for the used experiment
    #[clap(short, long, default_value="25")]
    num_iterations: usize,

    /// batch size for applicable algorithms
    #[clap(short, long, default_value="16")]
    batch_size: usize
}

#[derive(Parser)]
struct CustomExperimentCommand {

    /// name of the experiment to perfom
    #[clap(short, long)]
    name: String,

    /// path to the folder that stores used source and target nodes
    #[clap(short, long)]
    data_path: String,

    /// path to the folder that stores the used graph
    #[clap(short, long)]
    graph_path: String,

    /// source-target size configuration
    #[clap(short, long)]
    source_target_size: String,

    /// number of iterations for the used experiment
    #[clap(short, long, default_value="25")]
    num_iterations: usize,

    /// batch size for applicable algorithms
    #[clap(short, long, default_value="16")]
    batch_size: usize
}

#[derive(Parser)]
struct CreateDataCommand {

    /// path to folder that contains the graph used to create experiment data
    #[clap(short, long)]
    graph_path: String,

    /// path to output folder
    #[clap(short, long)]
    output_path: String, 

    /// number of files per size/ball_size
    #[clap(short, long)]
    num_iterations: usize,

    /// special configuration for data creation, e.g. single_node, ball_size_14_to_25
    #[clap(short, long)]
    configuration: String
}

fn main() {
    let opts: Opts = Opts::parse();

    // stores all available source target size configurations
    let source_target_configurations = SourceTargetConfiguration::all_source_target_configuration();
    let available_algorithms: Vec<ManyToManyExperimentAlgorithm> = ManyToManyExperimentAlgorithm::all_algorithms();
    let custom_experiments: Vec<Experiment> = Experiment::all_experiments();

    match opts.subcmd {
        SubCommand::CustomExperiment(experiment) => {
            let graph_path = Path::new(&experiment.graph_path);
            let data_path = Path::new(&experiment.data_path);

            let mut data = ExperimentData::new(graph_path, data_path);
            for custom_experiment in &custom_experiments {
                if custom_experiment.name == experiment.name {
                    for configuration in &source_target_configurations {
                        if *configuration.folder == experiment.source_target_size {
                            (custom_experiment.function)(&mut data, configuration, CustomExperimentVariables {num_iterations: experiment.num_iterations, batch_size: experiment.batch_size});
                        }
                    }
                    
                }
            }
        },
        SubCommand::Experiment(experiment) => {
            let graph_path = Path::new(&experiment.graph_path);
            let data_path = Path::new(&experiment.data_path);

            let mut data = ExperimentData::new(graph_path, data_path);
            
            for many_to_many_algorithm in available_algorithms {
                if many_to_many_algorithm.name == experiment.algorithm {
                    let mut algorithm = (many_to_many_algorithm.algorithm)(&mut data, ExperimentVariables {batch_size: experiment.batch_size});
                    
                    for configuration in &source_target_configurations {
                        if *configuration.folder == experiment.source_target_size {
                            for ball_size in &configuration.ball_sizes {
                                for iteration in 0..experiment.num_iterations {
                                    let (sources, targets) = get_sources_and_targets(&mut data.node_picker, iteration, *ball_size, &configuration.folder);

                                    let (select_time, _result) = measure_time(|| {
                                        algorithm.initialize(&sources, &targets);
                                        algorithm.select(&sources, &targets);
                                    });

                                    let (query_time, _result) = measure_time(|| algorithm.query(&sources, &targets));
        
                                    println!("{},{:?},{:?}", ball_size, select_time.as_nanos(), query_time.as_nanos());
                                }
                            } 
                        }
                    }
                }
            }
        },
        SubCommand::CreateData(create_data) => {
            let graph_path = Path::new(&create_data.graph_path);
            let output_path = Path::new(&create_data.output_path);

            let (fwd_first_edge, fwd_head, fwd_weight, bwd_first_edge, bwd_head, bwd_weight, _ranks) = read_graph(graph_path);
            let (first_edge, head, weight) = combine_restricted_graph(fwd_first_edge.len() - 1, &fwd_first_edge, &fwd_head, &fwd_weight, &bwd_first_edge, &bwd_head, &bwd_weight);
            let mut node_picker = NodePicker::new(first_edge.len() -1, first_edge, head, weight);
        
            match create_data.configuration.as_str() {
                "many_to_many_sources_equals_targets" => {
                    node_picker.create_many_to_many_sources_equals_targets(create_data.num_iterations, &output_path.join("sources_equals_targets"));
                },
                "many_to_many_same_ball" => {
                    node_picker.create_many_to_many_same_ball(create_data.num_iterations, &output_path.join("same_balls"));
                },
                "many_to_many_different_balls" => {
                    node_picker.create_many_to_many_separate_balls(create_data.num_iterations, &output_path.join("different_balls"));
                },
                "one_to_many_same_balls" => {
                    node_picker.create_one_to_many_same_balls(create_data.num_iterations, &output_path.join("same_balls"));
                }
                "one_to_many_different_balls" => {
                    node_picker.create_one_to_many_separate_balls(create_data.num_iterations, &output_path.join("different_balls"));
                }
                _ => panic!("unknown create data configuration")
            };    
        }
    }
}