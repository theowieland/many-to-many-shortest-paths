use shortest_path_algorithms::contraction::contraction_hierarchies::contract_graph_with_order;
use shortest_path_algorithms::graph_representation::GraphArray;
use shortest_path_algorithms::many_to_many::forward_backward_buckets_merged::ForwardBackwardBucketsMergedManyToMany;
use shortest_path_algorithms::many_to_many::hub_labels_partially_pruned::HubLabelsPartiallyPruned;
use shortest_path_algorithms::many_to_many::hub_labels_top_down::HubLabelsTopDown;
use shortest_path_algorithms::many_to_many::lazy_batched_buckets::LazyBatchedBuckets;
use shortest_path_algorithms::many_to_many::lazy_one_to_many::LazyRPHASTOneToMany;
use shortest_path_algorithms::many_to_many::lazy_simd::LazySIMD;
use shortest_path_algorithms::many_to_many::rphast_bwd_rank::RPHASTBwdRank;
use shortest_path_algorithms::many_to_many::rphast_multiple_trees_fwd_buckets::RPHASTMultipleTreesFwdBuckets;
use shortest_path_algorithms::many_to_many::rphast_multiple_trees_fwd_dijkstra::RPHASTMultipleTreesFwdDijkstra;
use shortest_path_algorithms::many_to_many::rphast_multiple_trees_fwd_rank::RPHASTMultipleTreesFwdRank;
use shortest_path_algorithms::many_to_many::rphast_simd::SIMDRPHAST;
use shortest_path_algorithms::types::*;
use shortest_path_algorithms::graph_permutations::convert_ordering_to_ranks;
use shortest_path_algorithms::many_to_many::many_to_many_algorithm::ManyToManyAlgorithm;
use shortest_path_algorithms::many_to_many::lazy::LazyRPHAST;
use shortest_path_algorithms::many_to_many::bucket::BucketManyToMany;
use shortest_path_algorithms::many_to_many::advanced_bucket::AdvancedBucketManyToMany;
use shortest_path_algorithms::many_to_many::rphast::RPHAST;
use shortest_path_algorithms::utils::create_restricted_graph;

//       ┌─┐
//  ┌────┤3├────┐
//  │    └┬┘    │
//  │     │1    │
//  │     │     │
// 1│    ┌┴┐    │1
//  │    │1│    │
//  │  4 ├─┤ 4  │
//  │ ┌──┘ └──┐ │
//  │ │       │ │
//  ├─┤       ├─┤
//  │0│       │2│
//  └─┘       └─┘
fn get_test_graph() -> (Vec<EdgeId>, Vec<NodeId>, Vec<Weight>, Vec<u32>, Vec<Weight>) {
    let first_out = vec![0, 2, 5, 7, 10];
    let head = vec![1, 3, 0, 2, 3, 1, 3, 0, 1, 2];
    let weight = vec![4, 1, 4, 4, 1, 4, 1, 1, 1, 1];
    let order = vec![0, 2, 1, 3];
    
    let expected_many_to_many_result = vec![0, 2, 2, 1, 2, 0, 2, 1, 2, 2, 0, 1, 1, 1, 1, 0]; 

    (first_out, head, weight, order, expected_many_to_many_result)
}

fn get_bwd_stall_on_demand_test_graph() -> (Vec<EdgeId>, Vec<NodeId>, Vec<Weight>, Vec<EdgeId>, Vec<NodeId>, Vec<Weight>, Vec<usize>, Vec<Weight>) {
    let first_out = vec![0, 1, 1, 3, 4, 6, 7];
    let head = vec![4, 1, 3, 1, 2, 5, 2];
    let weight = vec![4, 8, 3, 3, 1, 2, 3];
    let order = vec![0, 1, 2, 3, 4, 5];
    let graph: GraphArray = GraphArray::new(first_out, head, weight);
    let (first_edge, target_node, weight, rank) = contract_graph_with_order(&graph, &order);

    let expected_many_to_many_result = vec![
        0, 11, 5, 8, 4, 6, 
        INFINITY, 0, INFINITY, INFINITY, INFINITY, INFINITY,
        INFINITY, 6, 0, 3, INFINITY, INFINITY,
        INFINITY, 3, INFINITY, 0, INFINITY, INFINITY,
        INFINITY, 7, 1, 4, 0, 2,
        INFINITY, 9, 3, 6, INFINITY, 0
    ];

    let (fwd_first_edge, fwd_target_node, fwd_weight, bwd_first_edge, bwd_source_node, bwd_weight) = create_restricted_graph(&first_edge, &target_node, &weight, &rank);
    (fwd_first_edge, fwd_target_node, fwd_weight, bwd_first_edge, bwd_source_node, bwd_weight, rank, expected_many_to_many_result)
}

fn get_double_shortcut_graph() -> (Vec<EdgeId>, Vec<NodeId>, Vec<Weight>, Vec<EdgeId>, Vec<NodeId>, Vec<Weight>, Vec<usize>, Vec<Weight>){
    let first_out = vec![0, 2, 2, 4, 5, 6];
    let head = vec![2, 3, 1, 4, 2, 1];
    let weight = vec![10, 2, 10, 2, 2, 2];
    let order = vec![0, 1, 2, 3, 4];
    let graph: GraphArray = GraphArray::new(first_out, head, weight);
    let (first_edge, target_node, weight, rank) = contract_graph_with_order(&graph, &order);

    let expected_many_to_many_result = vec![
        0, 8, 4, 2, 6,
        INFINITY, 0, INFINITY, INFINITY, INFINITY,
        INFINITY, 4, 0, INFINITY, 2,
        INFINITY, 6, 2, 0, 4,
        INFINITY, 2, INFINITY, INFINITY, 0
    ];

    let (fwd_first_edge, fwd_target_node, fwd_weight, bwd_first_edge, bwd_source_node, bwd_weight) = create_restricted_graph(&first_edge, &target_node, &weight, &rank);
    (fwd_first_edge, fwd_target_node, fwd_weight, bwd_first_edge, bwd_source_node, bwd_weight, rank, expected_many_to_many_result)
}

/// tests the many to many algorithms on the get_test_graph
#[test]
fn test_many_to_many() {
    let (first_out, head, weight, ordering, expected_many_to_many_result) = get_test_graph();
    let ranks = convert_ordering_to_ranks(&ordering);
    let (fwd_first_edge, fwd_head, fwd_weight, bwd_first_edge, bwd_head, bwd_weight) = create_restricted_graph(&first_out, &head, &weight, &ranks);

    let many_to_many_algorithms = get_all_many_to_many_algorithms(&fwd_first_edge, &fwd_head, &fwd_weight, &bwd_first_edge, &bwd_head, &bwd_weight, &ranks);

    let sources = vec![0, 1, 2, 3];
    let targets = vec![0, 1, 2, 3];

    for mut algorithm in many_to_many_algorithms {
        algorithm.calculate(&sources, &targets);

        let distance_table = algorithm.get_distance_table();

        println!("{:?}", distance_table.data);

        for source_index in 0..sources.len() {
            for target_index in 0..targets.len() {
                assert_eq!(distance_table.get(source_index, target_index), expected_many_to_many_result[source_index * sources.len() + target_index]);
            }
        }
    }
}

#[test]
fn test_many_to_many_stall_on_demand() {
    let (fwd_first_edge, fwd_head, fwd_weight, bwd_first_edge, bwd_head, bwd_weight, ranks, expected_many_to_many_result) = get_bwd_stall_on_demand_test_graph();

    let many_to_many_algorithms = get_all_many_to_many_algorithms(&fwd_first_edge, &fwd_head, &fwd_weight, &bwd_first_edge, &bwd_head, &bwd_weight, &ranks);

    let sources = vec![0, 1, 2, 3, 4, 5];
    let targets = vec![0, 1, 2, 3, 4, 5];

    for mut algorithm in many_to_many_algorithms {
        algorithm.calculate(&sources, &targets);

        let distance_table = algorithm.get_distance_table();

        println!("{:?}", distance_table.data);

        for source_index in 0..sources.len() {
            for target_index in 0..targets.len() {
                assert_eq!(distance_table.get(source_index, target_index), expected_many_to_many_result[source_index * sources.len() + target_index]);
            }
        }
    }
}

#[test]
fn test_double_shortcut_graph() {
    let (fwd_first_edge, fwd_head, fwd_weight, bwd_first_edge, bwd_head, bwd_weight, ranks, expected_many_to_many_result) = get_double_shortcut_graph();

    let many_to_many_algorithms = get_all_many_to_many_algorithms(&fwd_first_edge, &fwd_head, &fwd_weight, &bwd_first_edge, &bwd_head, &bwd_weight, &ranks);

    let sources = vec![0, 1, 2, 3, 4];
    let targets = vec![0, 1, 2, 3, 4];

    for mut algorithm in many_to_many_algorithms {
        algorithm.calculate(&sources, &targets);

        let distance_table = algorithm.get_distance_table();

        println!("{:?}", distance_table.data);

        for source_index in 0..sources.len() {
            for target_index in 0..targets.len() {
                assert_eq!(distance_table.get(source_index, target_index), expected_many_to_many_result[source_index * sources.len() + target_index]);
            }
        }
    }
}

fn get_all_many_to_many_algorithms(fwd_first_edge: &Vec<EdgeId>, fwd_head: &Vec<NodeId>, fwd_weight: &Vec<Weight>, bwd_first_edge: &Vec<EdgeId>, bwd_head: &Vec<NodeId>, bwd_weight: &Vec<Weight>, ranks: &Vec<usize>) -> Vec<Box<dyn ManyToManyAlgorithm>> {
    let many_to_many_algorithms: Vec<Box<dyn ManyToManyAlgorithm>> = vec![
        Box::new(AdvancedBucketManyToMany::new(fwd_first_edge, fwd_head, fwd_weight, bwd_first_edge, bwd_head, bwd_weight, ranks)),   
        Box::new(BucketManyToMany::new(fwd_first_edge, fwd_head, fwd_weight, bwd_first_edge, bwd_head, bwd_weight, ranks)),
        Box::new(ForwardBackwardBucketsMergedManyToMany::new(fwd_first_edge, fwd_head, fwd_weight, bwd_first_edge, bwd_head, bwd_weight, ranks)),
        Box::new(HubLabelsPartiallyPruned::new(fwd_first_edge, fwd_head, fwd_weight, bwd_first_edge, bwd_head, bwd_weight)),
        Box::new(HubLabelsTopDown::new(fwd_first_edge, fwd_head, fwd_weight, bwd_first_edge, bwd_head, bwd_weight, ranks)),
        Box::new(LazyBatchedBuckets::new(16, fwd_first_edge, fwd_head, fwd_weight, bwd_first_edge, bwd_head, bwd_weight, ranks)),
        Box::new(LazyRPHASTOneToMany::new(fwd_first_edge, fwd_head, fwd_weight, bwd_first_edge, bwd_head, bwd_weight)),
        Box::new(LazyRPHAST::new(fwd_first_edge, fwd_head, fwd_weight, bwd_first_edge, bwd_head, bwd_weight)),
        Box::new(LazySIMD::new(fwd_first_edge, fwd_head, fwd_weight, bwd_first_edge, bwd_head, bwd_weight, ranks)),
        Box::new(RPHASTBwdRank::new(16, fwd_first_edge, fwd_head, fwd_weight, bwd_first_edge, bwd_head, bwd_weight, ranks)),
        Box::new(RPHASTMultipleTreesFwdBuckets::new(16, fwd_first_edge, fwd_head, fwd_weight, bwd_first_edge, bwd_head, bwd_weight, ranks)),
        Box::new(RPHASTMultipleTreesFwdDijkstra::new(16, fwd_first_edge, fwd_head, fwd_weight, bwd_first_edge, bwd_head, bwd_weight)),
        Box::new(RPHASTMultipleTreesFwdRank::new(16, fwd_first_edge, fwd_head, fwd_weight, bwd_first_edge, bwd_head, bwd_weight, ranks)),
        Box::new(SIMDRPHAST::new(fwd_first_edge, fwd_head, fwd_weight, bwd_first_edge, bwd_head, bwd_weight, ranks)),
        Box::new(RPHAST::new(fwd_first_edge, fwd_head, fwd_weight, bwd_first_edge, bwd_head, bwd_weight, ranks)),
    ];

    many_to_many_algorithms
}