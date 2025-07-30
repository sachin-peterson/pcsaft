// PC-SAFT Model
#[derive(Debug, Clone)]
pub struct PCSAFTModel {
    pub params: PCSAFTParams<f64>,
    pub assoc: AssocOptions
}

// SAFT Component
#[derive(Debug, Clone)]
pub struct PCSAFTComponent {
    pub mw: f64,                                // Molecular weight      [g/mol]
    pub m: f64,                                 // No. of segments       []
    pub sigma: f64,                             // Segment diameter      [Angstrom]
    pub epsilon: f64,                           // Dispersion energy     [K]
    pub bondvol: f64,                           // Assodication energy   [K]
    pub epsilon_assoc: f64,                     // Association volume    [m^3]
    pub sites: Vec<String>                      // Association sites     []
}

// PC-SAFT Parameters
#[derive(Debug, Clone)]
pub struct PCSAFTParams<T> {
    pub mw: Vec<T>,                             // Molecular weight      [g/mol]
    pub m: Vec<T>,                              // No. of segments       []
    pub sigma: Vec<Vec<T>>,                     // Segment diameter      [Angstrom]
    pub epsilon: Vec<Vec<T>>,                   // Dispersion energy     [K]
    pub bondvol: AssocMatrix<T>,                // Assodication energy   [K]
    pub epsilon_assoc: AssocMatrix<T>,          // Association volume    [m^3]
    pub kij: Option<Vec<Vec<T>>>,               // Interaction parameter []
    pub lij: Option<Vec<Vec<T>>>                // Interaction parameter []
}

// Association matrix
#[derive(Debug, Clone)]
pub struct AssocMatrix<T> {
    pub values: Vec<T>,                         // Parameter values
    pub outer_indices: Vec<(usize, usize)>,     // Index of components
    pub inner_indices: Vec<(usize, usize)>,     // Index of sites
    pub outer_size: (usize, usize),             // Size of component matrix
    pub inner_size: (usize, usize)              // Size of site matrices
}

// Association parameters
#[derive(Debug, Clone)]
pub struct AssocOptions {
    pub sites: SiteParam,                       // Association sites
    pub alpha: f64,                             // Damping factor
    pub atol: f64,                              // Absolute tolerance
    pub rtol: f64,                              // Relative tolerance
    pub max_iters: usize,                       // For Newton method
    pub combining: String                       // Mixing rules
}

// Association site scheme
#[derive(Debug, Clone)]
pub struct SiteParam {
    pub components: Vec<String>,                // Component names
    pub sites: Vec<Vec<String>>,                // Site labels for each component
    pub n_sites: PackedVector<usize>            // Flat site count per site index
}

// Packed vector
#[derive(Debug, Clone)]
pub struct PackedVector<T> {
    pub v: Vec<T>,
    pub p: Vec<usize> 
}

impl<T: Clone> PackedVector<T> {
    // Returns the i-th row
    pub fn get_row(&self, i: usize) -> &[T] {
        let start = self.p[i];
        let end = self.p[i+1];
        &self.v[start..end]
    }

    // Returns a single element
    pub fn get(&self, i: usize, j: usize) -> T {
        self.get_row(i)[j].clone()
    }
}

// Packed vector of vectors
#[derive(Debug, Clone)]
pub struct PackedVecVec<T> {
    pub v: Vec<Vec<T>>,
    pub p: Vec<usize>,
}

impl<T: Clone> PackedVecVec<T> {
    // Creates a new PackedVecVec
    pub fn new(p: Vec<usize>, flat: Vec<T>) -> Self {
        let mut v = Vec::with_capacity(p.len());
        for i in 0..p.len() {
            let start = p[i];
            let end = if i + 1 < p.len() { p[i + 1] } else { flat.len() };
            v.push(flat[start..end].to_vec());
        }
        Self {p, v}
    }
}

// Methods for AssocMatrix
impl AssocMatrix<f64> {
    // Returns the value at [i][j][a][b]
    pub fn get(&self, i: usize, j: usize, a: usize, b: usize) -> f64 {
        for (idx, &(ii, jj)) in self.outer_indices.iter().enumerate() {
            if ii == i && jj == j {
                let (aa, bb) = self.inner_indices[idx];
                if aa == a && bb == b {
                    return self.values[idx];
                }
            }
        }
        0.0
    }

    // Creates a new AssocMatrix
    pub fn from_parts(values: Vec<f64>, outer_indices: Vec<(usize, usize)>, inner_indices: Vec<(usize, usize)>) -> Self {
        let outer_size = outer_indices.iter()
        .fold((0, 0), |(max_i, max_j), &(i, j)| {
            (max_i.max(i+1), max_j.max(j+1))
        });

        let inner_size = inner_indices.iter()
        .fold((0, 0), |(max_a, max_b), &(a, b)| {
            (max_a.max(a+1), max_b.max(b+1))
        });

        Self {values, outer_indices, inner_indices, outer_size, inner_size}
    }
}
