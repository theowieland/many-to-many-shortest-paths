use shortest_path_algorithms::utils::data_structures::{Matrix, LazyMaxHeap};


#[test]
fn test_matrix() {
    let mut matrix = Matrix::empty();
    let initial_value = 0;

    matrix.resize(10, 10, initial_value);
    assert_eq!(matrix.data.len(), 10 * 10);

    matrix.resize(10, 20, initial_value);
    assert_eq!(matrix.data.len(), 10 * 20);

    matrix.resize(10, 10, initial_value);
    assert_eq!(matrix.data.len(), 10 * 10);
}

#[test]
fn test_lazy_heap() {
    let mut lazy_heap = LazyMaxHeap::<usize>::new(10, 100);
    lazy_heap.reset(10, 100);
    for index in 0..5 {
        lazy_heap.decrease_value(index, index);
    }
    lazy_heap.decrease_value(9, 1);
    lazy_heap.decrease_value(8, 88);
    lazy_heap.decrease_value(7, 77);
    lazy_heap.decrease_value(6, 66);
    lazy_heap.decrease_value(5, 55);
    lazy_heap.recompute_max_element();
    lazy_heap.print();
}