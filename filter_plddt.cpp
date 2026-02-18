#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <omp.h>
#include <cstdlib>
#include <atomic>

namespace fs = std::filesystem;

inline double fast_bfactor(const char* line) {
    return std::atof(line + 60);
}

void print_help() {
    std::cout <<
    "filter_plddt\n\n"
    "Description:\n"
    "  Removes residues from AlphaFold PDB structures where pLDDT (stored in B-factor column)\n"
    "  is below the specified cutoff.\n\n"
    "Usage:\n"
    "  ./filter_plddt <input_dir> <output_dir> <cutoff> [threads]\n\n"
    "Arguments:\n"
    "  input_dir   Directory containing PDB files\n"
    "  output_dir  Directory where filtered PDB files will be written\n"
    "  cutoff      pLDDT threshold, e.g. 50\n"
    "  threads     Optional number of threads to use\n\n"
    "Example:\n"
    "  ./filter_plddt_fast af_models filtered_models 50 112\n\n"
    "Notes:\n"
    "  pLDDT values are read from columns 61–66 of ATOM/HETATM records.\n"
    "  Files are processed in parallel using OpenMP.\n"
    "  Ridvan A. Ayaz - razizayaz@gmail.com\n "
    "  ! This program comes with ZERO WARRANTY. Always do your QC and take backups !\n" ;
}

int main(int argc, char* argv[]) {

    if (argc == 2) {
        std::string arg = argv[1];
        if (arg == "-h" || arg == "--help") {
            print_help();
            return 0;
        }
    }

    if (argc < 4 || argc > 5) {
        std::cerr << "Invalid arguments.\n"
        "Usage: ./filter_plddt <input_dir> <output_dir> <cutoff> [threads]\n"
        " Use -h for help.\n";
        return 1;
    }

    std::string input_dir = argv[1];
    std::string output_dir = argv[2];
    double cutoff = std::atof(argv[3]);

    if (argc == 5) {
        int threads = std::atoi(argv[4]);
        if (threads > 0) {
            omp_set_num_threads(threads);
        }
    }

    fs::create_directories(output_dir);

    std::vector<fs::path> pdb_files;
    pdb_files.reserve(25000);

    for (auto& p : fs::recursive_directory_iterator(input_dir)) {
        if (p.path().extension() == ".pdb") {
            pdb_files.push_back(p.path());
        }
    }

    size_t total = pdb_files.size();
    std::cout << "Found " << total << " PDB files\n";

    std::atomic<size_t> progress(0);
    const size_t report_interval = 100;

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < total; ++i) {

        const fs::path& infile = pdb_files[i];
        fs::path outfile = fs::path(output_dir) / infile.filename();

        std::ifstream fin(infile, std::ios::binary);
        if (!fin) continue;

        fin.seekg(0, std::ios::end);
        size_t size = fin.tellg();
        fin.seekg(0);

        std::string buffer;
        buffer.resize(size);
        fin.read(buffer.data(), size);
        fin.close();

        std::string output;
        output.reserve(size);

        const char* data = buffer.c_str();
        size_t pos = 0;

        while (pos < size) {

            size_t line_end = buffer.find('\n', pos);
            if (line_end == std::string::npos) line_end = size;

            const char* line = data + pos;
            size_t line_len = line_end - pos;

            if (line_len > 66 &&
                ((line[0]=='A' && line[1]=='T' && line[2]=='O' && line[3]=='M') ||
                 (line[0]=='H' && line[1]=='E' && line[2]=='T' && line[3]=='A'))) {

                double b = fast_bfactor(line);

                if (b >= cutoff) {
                    output.append(line, line_len);
                    output.push_back('\n');
                }

            } else {
                output.append(line, line_len);
                output.push_back('\n');
            }

            pos = line_end + 1;
        }

        std::ofstream fout(outfile, std::ios::binary);
        fout.write(output.data(), output.size());
        fout.close();

        size_t done = ++progress;

        if (done % report_interval == 0) {
            #pragma omp critical
            {
                double percent = (100.0 * done) / total;
                std::cout << "Processed " << done << "/" << total
                          << " (" << percent << "%)\n";
            }
        }
    }

    std::cout << "Completed processing " << total << " files\n";
    return 0;
}
