use std;

pub type NodeId = u32;
pub type EdgeId = u32;
pub type Weight = u32;
pub type Rank = usize;

pub const INFINITY: Weight = std::u32::MAX / 2;

pub type NodeIds = Vec<NodeId>;
pub type EdgeIds = Vec<EdgeId>;
pub type Weights = Vec<Weight>;
pub type Arclist = [(NodeId, Weight)];
pub type Ranks = Vec<Rank>;