# filter_pLDDT

**Removes residues from AlphaFold PDB structures where pLDDT (stored in B-factor column)
  is below the specified cutoff.**

**Install:**
1- Download filter_plddt.cpp
2- run :  *g++ -O3 -march=native -flto -fopenmp filter_plddt.cpp -o filter_plddt*
   
Requires OPENMP

 **Usage:**
  ./filter_plddt <input_dir> <output_dir> <cutoff> [threads]

**Example:**
  ./filter_plddt_fast af_models filtered_models 50 112

 **Assumptions**
	pLDDT values are stored in columns 61–66 of ATOM/HETATM records
	Input files follow standard PDB formatting
	AlphaFold models store identical pLDDT values for all atoms of a residue

**Arguments:**
  input_dir   Directory containing PDB files
  output_dir  Directory where filtered PDB files will be written
  cutoff      pLDDT threshold, e.g. 50
  threads     Optional number of threads to use. Defaults to the number of hardware cores. 

