#![feature(portable_simd)]
#![feature(slice_split_at_unchecked)]

extern crate rand;
extern crate core_simd;
extern crate rudac;
extern crate clap;

pub mod types;

pub mod graph_representation;
pub mod witness_search;
pub mod contraction_hierarchies;
pub mod graph_algorithms;
pub mod utils;
pub mod data_structures;
pub mod bucket_containers;
pub mod graph_permutations;

pub mod many_to_many;
pub mod experiments;