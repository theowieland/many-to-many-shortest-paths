use shortest_path_algorithms::utils::binary_heap::{HeapElement, MinBinaryHeap};

#[derive(Copy, Clone, Eq, PartialEq, Debug, Ord, PartialOrd)]
pub struct MinElement {
    pub value: usize,
    pub unique_id: usize
}

impl HeapElement for MinElement {

    fn unique_index(&self) -> usize {
        self.unique_id as usize
    }
}

#[test]
fn insert_pop_test() {
    let mut min_heap = MinBinaryHeap::new(10);

    min_heap.insert(MinElement {value: 10, unique_id: 0});
    min_heap.insert(MinElement {value: 8, unique_id: 1});
    min_heap.insert(MinElement {value: 12, unique_id: 2});
    min_heap.insert(MinElement {value: 7, unique_id: 3});

    assert_eq!(min_heap.len(), 4);

    assert_eq!(min_heap.pop(), Some(MinElement {value: 7, unique_id: 3}));
    assert_eq!(min_heap.pop(), Some(MinElement {value: 8, unique_id: 1}));
    assert_eq!(min_heap.pop(), Some(MinElement {value: 10, unique_id: 0}));
    assert_eq!(min_heap.pop(), Some(MinElement {value: 12, unique_id: 2}));
}

#[test]
fn decrease_test() {
    let mut min_heap = MinBinaryHeap::new(10);
    min_heap.insert(MinElement {value: 10, unique_id: 0});
    min_heap.print();

    min_heap.insert(MinElement {value: 8, unique_id: 1});
    min_heap.print();

    min_heap.insert(MinElement {value: 12, unique_id: 2});
    min_heap.print();

    min_heap.insert(MinElement {value: 7, unique_id: 3});
    min_heap.print();

    min_heap.decrease(MinElement {value: 4, unique_id: 2});
    min_heap.print();

    assert_eq!(min_heap.len(), 4);
    assert_eq!(min_heap.pop(), Some(MinElement {value: 4, unique_id: 2}));
    assert_eq!(min_heap.pop(), Some(MinElement {value: 7, unique_id: 3}));
    assert_eq!(min_heap.pop(), Some(MinElement {value: 8, unique_id: 1}));
    assert_eq!(min_heap.pop(), Some(MinElement {value: 10, unique_id: 0}));
    
}