use std::{ops::Range, cmp};

pub struct MinBinaryHeap<T: HeapElement + Ord> {
    data: Vec<T>,
    indices: Vec<usize>
}

pub trait HeapElement {

    fn unique_index(&self) -> usize;
}

const INVALID_POSITION: usize = usize::MAX;

impl<T: HeapElement + Ord> MinBinaryHeap<T> {

    pub fn new(size: usize) -> Self {
        MinBinaryHeap { 
            data: Vec::new(), 
            indices: vec![INVALID_POSITION; size]
        }
    } 

    pub fn len(&self) -> usize {
        self.data.len()    
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn insert(&mut self, entry: T) {
        let initial_index = self.data.len();
        let entry_unique_id = entry.unique_index();
        self.data.push(entry);
        self.indices[entry_unique_id] = initial_index;

        self.sift_up(initial_index);
    }

    pub fn insert_or_decrease(&mut self, entry: T) {
        if self.indices[entry.unique_index()] != INVALID_POSITION {
            self.decrease(entry);
        }
        else {
            self.insert(entry);
        }
    }

    pub fn decrease(&mut self, entry: T) {
        let unique_id = entry.unique_index();

        self.data[self.indices[unique_id]] = entry;
        self.sift_up(self.indices[unique_id]);
    }

    pub fn get_min(&mut self) -> Option<&T> {
        self.data.first()
    }

    pub fn clear(&mut self) {
        for element in &self.data {
            self.indices[element.unique_index()] = INVALID_POSITION;
        }

        self.data.clear();
    }

    pub fn pop(&mut self) -> Option<T> {
        if self.len() == 0 {
            return None;
        }

        self.swap(0, self.len() - 1); // move first element to the back
        let min = self.data.pop();
        let min_element = min.unwrap();
        let min_unique_id = min_element.unique_index();
        self.indices[min_unique_id] = INVALID_POSITION;

        if self.len() > 0 {
            self.sift_down(0);
        }

        Some(min_element)
    } 

    pub fn get_by_unique_index(&mut self, unique_index: usize) -> Option<&mut T> {
        let position = self.indices[unique_index];

        Some(&mut self.data[position])
    }

    pub fn remove_by_unique_index(&mut self, unique_index: usize) -> Option<T> {
        let position = self.indices[unique_index];

        self.swap(position, self.len() -1); // move last element to current postion
        let element = self.data.pop();
        self.indices[unique_index] = INVALID_POSITION;

        self.sift_down(position);

        element
    }

    pub fn contains_unique_index(&mut self, unique_index: usize) -> bool {
        self.indices[unique_index] != INVALID_POSITION
    }

    fn sift_up(&mut self, index: usize) {
        let mut current_index = index;

        while 0 < current_index {
            let parent_index = self.parent_index(current_index);

            if self.data[parent_index] <= self.data[current_index] {
                break;
            }
            
            // switch elements
            self.swap(parent_index, current_index);
            current_index = parent_index;
        }
    }

    fn sift_down(&mut self, index: usize) {
        let mut current_index = index;

        while current_index < self.len() {
            let mut min_entry_index: usize = current_index;

            for child_index in self.children_indices(current_index) {
                if child_index < self.len() && self.data[child_index] < self.data[min_entry_index] {
                    min_entry_index = child_index;
                }
            }

            if min_entry_index != current_index {
                self.swap(min_entry_index, current_index);

                current_index = min_entry_index;
            }
            else {
                break;
            }
        }
    }

    pub fn print(&self) {
        println!("min heap: {}", self.data.len());
        for index in 0..self.len() {
            println!("{} {} ", index, self.data[index].unique_index());
        }
    }

    fn parent_index(&mut self, index: usize) -> usize {
        (index - 1) / 2
    }

    fn children_indices(&mut self, index: usize) -> Range<usize> {
        let min_index = cmp::min(self.len(), 2 * index + 1);
        let max_index = cmp::min(self.len(), 2 * index + 2);

        min_index..(max_index + 1)
    }

    fn swap(&mut self, first: usize, second: usize) {
        self.indices[self.data[first].unique_index()] = second;
        self.indices[self.data[second].unique_index()] = first;

        self.data.swap(first, second);
    }

}