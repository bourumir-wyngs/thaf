use std::fmt;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum Strand {
    Plus,
    Minus,
}

impl Strand {
    pub fn from_char(x: char) -> Strand {
        match x {
            '+' => Strand::Plus,
            '-' => Strand::Minus,
            _ => panic!("Invalid strand [{x}]"),
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Strand::Plus => write!(f, "+"),
            Strand::Minus => write!(f, "-"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Region {
    pub id: String,
    pub start: usize,
    pub end: usize,
    pub strand: Strand,
}

#[derive(Debug, Clone)]
pub struct Transcript {
    pub id: String,
    pub chromosome: String, // added chromosome field
    pub regions: Vec<Region>,
}

impl Transcript {
    /// Computes the total length of all regions in the transcript.
    pub fn size(&self) -> usize {
        self.regions.iter().map(|r| r.end - r.start + 1).sum()
    }
}

#[derive(Debug, Clone)]
pub struct TranscriptRegion {
    pub chromosome: String,
    pub start: usize,
    pub end: usize,
    pub strand: Strand,
    pub transcript_id: String,
    pub region_id: String,
    pub gene_id: Option<String>,
}
