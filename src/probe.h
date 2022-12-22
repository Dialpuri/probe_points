//
// Created by jordan on 12/14/22.
//

#include <iostream>
#include <string>
#include <fstream>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>
#include <clipper/core/container_map.h>

#include <stdexcept>
#include <array>


#ifndef PROBE_POINTS_PROBE_H
#define PROBE_POINTS_PROBE_H


class LibraryItem {

public:
    LibraryItem(std::string pdb_code);
    const std::string pdb_base_dir = "./data/test_library/";
    clipper::MiniMol load_pdb();

    clipper::Cell return_cell(clipper::MMonomer& monomer);

    clipper::Coord_orth calculate_center_point(std::vector<clipper::MAtom> &atoms);
    clipper::RTop_orth align_fragment(clipper::MMonomer& monomer);

    void convert_map_to_array(clipper::Xmap<float> &xmap);
    void calculate_electron_density(clipper::MiniMol& test_mol);

    void dump_minimol(clipper::MiniMol& output_model, std::string file_path, std::string file_name);
    void dump_electron_density(std::string path);
    std::vector<std::pair<clipper::MMonomer, clipper::Xmap<float>>> m_density;

private:
    std::string m_pdb_code;
    std::string m_pdb_file_path;
};

class Library {

public:
    Library();
    Library(std::string library_file_path);
    std::vector <LibraryItem> read_library_item();

    void combine_density();

private:
    std::vector <LibraryItem> m_library;
    std::string m_library_path;
};

class MapReader {
public:
    MapReader();

    void _test();
};


#endif //PROBE_POINTS_PROBE_H
