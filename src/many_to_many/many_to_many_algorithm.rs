use crate::{types::*, utils::data_structures::Matrix};

/// a many to many algorith is able to compute shortest path distances between a given set of sources and targets
pub trait ManyToManyAlgorithm {

    fn initialize(&mut self, sources: &[NodeId], targets: &[NodeId]);

    /// selection phase used by some algorihtms - some select targets others may select sources or nothing at all
    fn select(&mut self, sources: &[NodeId], targets: &[NodeId]);

    /// caluclate all shortest path distance pairs and store result in distance matrix
    fn query(&mut self, sources: &[NodeId], targets: &[NodeId]);

    /// caluclate the shortest path between the given source and target nodes
    fn calculate(&mut self, sources: &[NodeId], targets: &[NodeId]) {
        self.initialize(sources, targets);
        self.select(sources, targets);
        self.query(sources, targets);
    }

    /// return the previously calculated distance table
    fn get_distance_table(&mut self) -> &Matrix<Weight>;
}