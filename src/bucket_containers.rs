use crate::data_structures::EmptyValidFlags;

pub struct BucketContainer<T: Copy + Clone> {
    buckets: Vec<Vec<T>>,
    modified_buckets: Vec<usize>,
    modified_valid_flags: EmptyValidFlags
}

impl<T: Copy + Clone> BucketContainer<T> {

    pub fn new(num_buckets: usize) -> Self {
        let mut buckets: Vec<Vec<T>> = Vec::new();

        for _ in 0..num_buckets {
            buckets.push(Vec::new());
        }

        BucketContainer {
            buckets,
            modified_buckets: Vec::new(),
            modified_valid_flags: EmptyValidFlags::new(num_buckets)
        }
    }

    pub fn get_bucket(&self, index: usize) -> &Vec<T> {
        &self.buckets[index]
    }

    pub fn add_bucket_entry(&mut self, index: usize, value: T) {
        self.buckets[index].push(value);

        if !self.modified_valid_flags.is_valid(index) {
            self.modified_buckets.push(index);
            self.modified_valid_flags.set_valid(index);
        }
    }

    pub fn clear_bucket(&mut self, index: usize) {
        self.buckets[index].clear();
    }

    pub fn bucket_retain<F>(&mut self, index: usize, retain: F) where F: FnMut(&T) -> bool {
        self.buckets[index].retain(retain);
    }

    pub fn get_modified_buckets(&self) -> &Vec<usize> {
        &self.modified_buckets
    }

    pub fn reset(&mut self) {
        for index in &self.modified_buckets {
            self.buckets[*index].clear();
        }

        self.modified_valid_flags.reset();
        self.modified_buckets.clear();
    }

    pub fn sort_all_buckets<F>(&mut self, extractor: F) where F: Fn(&T) -> usize {
        for bucket in &mut self.buckets {
            bucket.sort_by_key(|t| extractor(t));
        }
    }
}

pub struct PreparationAccessBucketContainer<T: Copy + Clone> {
    data: Vec<T>,
    first_entry: Vec<usize>,
    preparation_entries: Vec<Vec<T>>
}

impl<T: Copy + Clone> PreparationAccessBucketContainer<T> {

    pub fn new(size: usize) -> Self {
        PreparationAccessBucketContainer {
            data: Vec::new(),
            first_entry: Vec::new(),
            preparation_entries: vec![Vec::new(); size],
        }
    }

    pub fn add_preparation_entry(&mut self, index: usize, entry: T) {
        self.preparation_entries[index].push(entry);
    }

    pub fn get_preparation_bucket(&mut self, index: usize) -> &[T] {
        &self.preparation_entries[index]
    }

    pub fn prepare(&mut self) {
        let mut first_index = 0;

        for index in 0..self.preparation_entries.len() {
            self.first_entry.push(first_index);

            let to_append = &mut self.preparation_entries[index];
            first_index += to_append.len();

            self.data.append(to_append);
            to_append.clear();
        }

        self.first_entry.push(first_index);
    }

    pub fn get_bucket(&self, index: usize) -> &[T] {
        &self.data[self.first_entry[index]..self.first_entry[index + 1]]
    }

    pub fn reset(&mut self) {
        self.data.clear();
        self.first_entry.clear();
    }
}