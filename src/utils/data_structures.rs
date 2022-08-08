/// this module contains different data structures used by some of the algorithms

use std::ops::{Index, Add};
use std::clone::Clone;
use std::cmp::{Eq, min};
use std::fmt::Debug;

use crate::types::{INFINITY, Weight};

pub trait ArrayStructure<T: Clone + Copy + Eq>: Index<usize> {

    fn reset(&mut self);
    fn set(&mut self, index: usize, value: T);
}

/// ValidFlags can be used to efficiently store a large amount of values that need to be invalidated very quickly
/// the whole array can be invalidated with a single instruction without iterating over the whole array
pub struct ValidFlags<T: Clone + Copy + Eq> {
    valid_flags: Vec<usize>,
    valid_flag: usize,
    default_value: T,
    data: Vec<T>
}

impl<T: Clone + Copy + Eq> Index<usize> for ValidFlags<T> {

    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        if self.valid_flags[index] == self.valid_flag {
            return &self.data[index];
        }
        
        return &self.default_value;
    }
}

impl<T: Clone + Copy + Eq> ArrayStructure<T> for ValidFlags<T> {

    fn reset(&mut self) {
        if self.valid_flag == usize::MAX {
            self.valid_flag = 1;

            let current_len = self.data.len();
            self.data = vec![self.default_value; current_len];
            self.valid_flags = vec![0; current_len];
        }
        else {
            self.valid_flag += 1;
        }
    }

    fn set(&mut self, index: usize, value: T) {
        self.data[index] = value;
        self.valid_flags[index] = self.valid_flag;
    }
}

impl<T: Clone + Copy + Eq> ValidFlags<T> {

    pub fn new(size: usize, default_value: T) -> Self {
        ValidFlags {
            valid_flags: vec![0; size],
            valid_flag: 1,
            default_value,
            data: vec![default_value; size]
        }
    }

    pub fn is_valid(&self, index: usize) -> bool {
        self.valid_flags[index] == self.valid_flag
    }

    pub fn get_unchecked(&self, index: usize) -> T {
        self.data[index]
    }
}

/// array that keeps track of all modified indices and offers a reset method to reset all those modified values
pub struct ResettableArray<T: Clone + Copy + Eq> {
    default_value: T,
    modified_indicies: Vec<usize>, // stores the indicies of all modidified entries
    data: Vec<T>
}

impl<T: Clone + Copy + Eq> Index<usize> for ResettableArray<T> {

    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl<T: Clone + Copy + Eq> ArrayStructure<T> for ResettableArray<T> {

    fn reset(&mut self) {
        for index in &self.modified_indicies {
            self.data[*index] = self.default_value;
        }

        self.modified_indicies.clear();
    }

    fn set(&mut self, index: usize, value: T) {
        if self.data[index] == self.default_value {
            self.modified_indicies.push(index);
        }

        self.data[index] = value;
    }
}

impl<T: Clone + Copy + Eq> ResettableArray<T> {

    pub fn new(size: usize, default_value: T) -> Self {
        ResettableArray {
            default_value: default_value,
            modified_indicies: Vec::new(),
            data: vec![default_value; size]
        }
    }

    pub unsafe fn access_unsafe(&self, index: usize) -> T {
        *self.data.get_unchecked(index)
    }

    pub unsafe fn get_unchecked(&mut self, index: usize) -> &T {
        self.data.get_unchecked(index)
    }
}

/// stores a 2d matrix inside a single 1d array
pub struct Matrix<T: Clone + Copy + Ord + Add<Output = T>> {
    pub data: Vec<T>,

    rows: usize,
    cols: usize,
}

impl<T: Clone + Copy + Ord + Add<Output = T>> Matrix<T> {
    
    pub fn new(rows: usize, cols: usize, initial_value: T) -> Self {
        Matrix {
            data: vec![initial_value; rows * cols],
            rows,
            cols
        }
    }

    pub fn empty() -> Self {
        Matrix {
            data: Vec::new(),
            rows: 0,
            cols: 0
        }
    }

    pub fn get(&self, row: usize, col: usize) -> T {
        self.data[row * self.cols + col]
    }

    pub fn set(&mut self, row: usize, col: usize, value: T) {
        self.data[row * self.cols + col] = value;
    }

    pub unsafe fn set_unsafe(&mut self, row: usize, col: usize, value: T) {
        *self.data.get_unchecked_mut(row * self.cols + col) = value;
    }

    /// returns true in case the value has been reduced, false otherwise
    pub fn reduce_value(&mut self, row: usize, col: usize, value: T) -> bool {
        let index = row * self.cols + col;

        let current_value = &mut self.data[index];

        if value < *current_value {
            *current_value = value;
            return true;
        }

        false
    }

    pub unsafe fn reduce_buckets(&mut self, row: usize, buckets: &[(usize, T)], value: T) {
        let current_values = self.data.get_unchecked_mut((row * self.cols)..((row + 1) * self.cols));

        buckets.iter().for_each(|(index, bucket_distance)| {
            let current_value = current_values.get_unchecked_mut(*index);

            if *bucket_distance + value < *current_value {
                *current_value = *bucket_distance + value;
            }
        });
    }

    pub unsafe fn reduce_buckets_batch_offset(&mut self, row: usize, buckets: &[(usize, T)], value: T, target_batch_index_offset: usize) {
        let current_values = self.data.get_unchecked_mut((row * self.cols)..((row + 1) * self.cols));

        buckets.iter().for_each(|(index, bucket_distance)| {
            let current_value = current_values.get_unchecked_mut(*index + target_batch_index_offset);

            if *bucket_distance + value < *current_value {
                *current_value = *bucket_distance + value;
            }
        });
    }

    /// resizes the vec to the given dimension and sets all values to the given initial value
    pub fn resize(&mut self, rows: usize, cols: usize, initial_value: T) {
        self.data.resize(rows * cols, initial_value);
        
        self.data.iter_mut().for_each(|entry| *entry = initial_value);
        
        self.rows = rows;
        self.cols = cols;
    }
}

impl<T: Clone + Copy + Ord + Add<Output = T>> Clone for Matrix<T> {

    fn clone(&self) -> Self {
        Self { data: self.data.clone(), rows: self.rows, cols: self.cols }
    }

    fn clone_from(&mut self, source: &Self) {
        self.data = source.data.clone();
        self.rows = source.rows;
        self.cols = source.cols;
    }
}

impl<T: Clone + Copy + Ord + Add<Output = T>> PartialEq for Matrix<T> {
    
    fn ne(&self, other: &Self) -> bool {
        self.data != other.data
    }

    fn eq(&self, other: &Self) -> bool {
        self.data == other.data && self.rows == other.rows && self.cols == other.cols
    }
}

const TREE_ARITY: usize = 4;

//doest support element removal, only decrease and maintains a max node
pub struct LazyMaxHeap<T: Copy + Clone + Eq + PartialOrd + Ord + Debug> {
    size: usize,
    sorted_indices: Vec<usize>, //saves the indices in a heap structure
    current_values: Vec<T>, //the current values based on which the heap is sorted
    updated_values: Vec<T> //updates the current values on demand
}

impl<T: Copy + Clone + Eq + PartialOrd + Ord + Debug> LazyMaxHeap<T> {
    
    pub fn new(size: usize, initial_value: T) -> Self {
        let mut sorted_indices = vec![0; size];
        for index in 0..size {
            sorted_indices[index] = index;
        }

        LazyMaxHeap {
            size,
            sorted_indices,
            current_values: vec![initial_value; size], 
            updated_values: vec![initial_value; size]
        }
    }

    pub fn reset(&mut self, size: usize, initial_value: T) {
        self.size = size;

        let mut sorted_indices = vec![0; self.size];
        for index in 0..self.size {
            sorted_indices[index] = index;
        }

        self.current_values = vec![initial_value; self.size];
        self.updated_values = vec![initial_value; self.size];
    }

    pub fn get_current_max(&self) -> (usize, T) {
        let max_index = self.sorted_indices[0];

        (max_index, self.current_values[max_index])
    }

    pub fn decrease_value(&mut self, index: usize, value: T) {
        self.updated_values[index] = value;
    }

    fn move_down_in_tree(&mut self, index: usize) {
        let mut current_index = index;
        
        loop {
            if let Some(biggest_child_index) = self.children_index_range(current_index).max_by_key(|&child_index| self.current_values[self.sorted_indices[child_index]]) {
                if self.current_values[self.sorted_indices[biggest_child_index]] > self.current_values[self.sorted_indices[current_index]] {
                    self.sorted_indices.swap(biggest_child_index, current_index);
                    current_index = biggest_child_index;
                }
                else {
                    return;
                }
            }
            else {
                return;
            }
        }
    }

    fn children_index_range(&self, parent_index: usize) -> std::ops::Range<usize> {
        let first_child = TREE_ARITY * parent_index + 1;
        let last_child = min(TREE_ARITY * parent_index + TREE_ARITY + 1, self.size);
        first_child..last_child
    }

    pub fn recompute_max_element(&mut self) {
        let mut current_max_index = self.sorted_indices[0];

        while self.current_values[current_max_index] != self.updated_values[current_max_index] {
            self.current_values[current_max_index] = self.updated_values[current_max_index];
            self.move_down_in_tree(0);

            current_max_index = self.sorted_indices[0];
        }
    }

    pub fn print(&self) {
        let mut sorted_values = Vec::new();
        println!("indices: {:?}", self.sorted_indices);

        for index in &self.sorted_indices {
            sorted_values.push(self.current_values[*index]);
        }
        println!("values: {:?}", sorted_values);
    }
}


/// valid flags but without associated data
pub struct EmptyValidFlags {
    size: usize,
    valid_flags: Vec<usize>,
    valid_flag: usize
}

impl EmptyValidFlags {

    pub fn new(size: usize) -> Self {
        EmptyValidFlags {
            size,
            valid_flags: vec![0; size],
            valid_flag: 1
        }
    }

    //resets and resizes the valid flags array
    pub fn resize(&mut self, size: usize) {
        if self.size != size {
            self.size = size;
            self.valid_flags = vec![0; self.size];
            self.valid_flag = 1;
        }
    }

    pub fn is_valid(&self, index: usize) -> bool {
        self.valid_flags[index] == self.valid_flag
    }

    pub fn set_valid(&mut self, index: usize) {
        self.valid_flags[index] = self.valid_flag;
    }

    pub fn reset(&mut self) {
        if self.valid_flag == usize::MAX {
            self.valid_flag = 1;

            self.valid_flags = vec![0; self.size];
        }
        else {
            self.valid_flag += 1;
        }
    }
}

//memory efficient batched valid flags that reduces memory space by not initializing all batches
pub struct LazyBatchedValidFlags {
    size: usize,
    batch_size: usize,
    valid_flag: usize,
    valid_batches: Vec<usize>,
    data: Vec<Option<Box<[Weight]>>>,
    initialized_batches: Vec<usize>,
}

impl LazyBatchedValidFlags {

    pub fn new(size: usize, batch_size: usize) -> Self {
        LazyBatchedValidFlags {
            size,
            batch_size,
            valid_flag: 1,
            valid_batches: vec![0; size],
            data: vec![None; size],
            initialized_batches: Vec::new()
        }
    }

    pub fn is_valid(&self, index: usize) -> bool {
        self.valid_batches[index] == self.valid_flag
    }

    pub fn initialize_batch(&mut self, index: usize) {
        self.valid_batches[index] = self.valid_flag;
        self.data[index] = Some((0..self.batch_size).map(|_| INFINITY as Weight).collect());
        self.initialized_batches.push(index);
    }

    pub fn reset(&mut self) {
        if self.valid_flag == usize::MAX {
            self.valid_flag = 1;

            self.data = vec![None; self.size];
            self.valid_batches = vec![0; self.size];
        }
        else {
            self.valid_flag += 1;

            for modified_batch in &self.initialized_batches {
                self.data[*modified_batch] = None;
            }
            self.initialized_batches.clear();
        }
    }

    pub fn set(&mut self, index: usize, batch_index: usize, value: Weight) {
        self.data[index].as_mut().unwrap()[batch_index] = value;
    }

    pub fn get(&self, index: usize, batch_index: usize) -> Weight {
        self.data[index].as_ref().unwrap()[batch_index]
    }

    pub fn get_row(&mut self, index: usize) -> &mut Option<Box<[Weight]>> {
        &mut self.data[index]
    }

    pub fn reduce_batch(&mut self, index: usize, offset: Weight, reduce_index: usize) {
        for batch_index in 0..self.batch_size {

            self.data[index].as_mut().unwrap()[batch_index] = min(self.data[index].as_ref().unwrap()[batch_index], offset + self.data[reduce_index].as_ref().unwrap()[batch_index]);
        }
    }
}