//
// Created by jordan on 12/14/22.
//

#include "hash.h"


void Hasher::load(std::string file_path) {
    clipper::CCP4MTZfile mtz;
    mtz.set_column_label_mode( clipper::CCP4MTZfile::Legacy );

    mtz.open_read(file_path);

    clipper::HKL_info hkls;
    hkls.init( mtz.spacegroup(), mtz.cell(), mtz.resolution(), true );

    clipper::HKL_data<clipper::data32::F_sigF>  wrk_f ( hkls );
    clipper::HKL_data<clipper::data32::F_phi>   fphi( hkls );

    mtz.import_hkl_data( wrk_f ,"FP,SIGFP");
    mtz.import_hkl_data( fphi,  "FWT,PHWT" );

    mtz.close_read();

    clipper::Spacegroup cspg = hkls.spacegroup();
    clipper::Cell       cxtl = hkls.cell();
    clipper::Grid_sampling grid( cspg, cxtl, hkls.resolution() );
    clipper::Xmap<float>   xwrk( cspg, cxtl, grid );
    xwrk.fft_from( fphi );

    std::cout << std::endl;
    std::cout << " Spgr " << hkls.spacegroup().symbol_xhm() << std::endl;
    std::cout << hkls.cell().format() << std::endl;
    std::cout << " Nref " << hkls.num_reflections() << " " << fphi.num_obs() << std::endl;

    xmap = xwrk;

}


void Hasher::slice(float slice_index) {

    std::cout << "Slicing map with " << slice_index << std::endl;

    clipper::Xmap_base::Map_reference_coord iu, iv, iw;
    clipper::Cell cell = xmap.cell();


    clipper::Grid grid = clipper::Grid(cell.a(),cell.b(),cell.c());

    clipper::Coord_orth base_coord_orth = clipper::Coord_orth(0,0,0);
    clipper::Coord_grid base_coord = base_coord_orth.coord_frac(xmap.cell()).coord_grid(grid);
    clipper::Xmap_base::Map_reference_coord i0 = clipper::Xmap_base::Map_reference_coord(xmap, base_coord);

    std::cout  << cell.format() << std::endl;
    std::cout << grid.format() << std::endl;

//    clipper::Coord_orth end_coord_orth = clipper::Coord_orth(cell.a(),cell.b(),cell.c());
//    clipper::Coord_grid end_coord = end_coord_orth.coord_frac(xmap.cell()).coord_grid(grid);

//    std::cout << "END COORD - " << end_coord_orth.format() << std::endl;
    clipper::Coord_grid end_coord = clipper::Coord_grid(cell.a(), cell.b(), cell.c());
//    clipper::Xmap_base::Map_reference_coord test = clipper::Xmap_base::Map_reference_coord(xmap, x);

    clipper::Xmap_base::Map_reference_coord iend = clipper::Xmap_base::Map_reference_coord(xmap, end_coord);

    int nu = grid.nu();
    int nv = grid.nv();

//    std::cout << "Nu " << nu << " Nv " << nv << "\n";

    for (int i = 0; i <= nu; i++) {
        std::vector<float> i_list = {};
        for (int j = 0; j <= nv; j++) {
            i_list.push_back(0.0f);
        }
        m_grid_values.push_back(i_list);
    }

//    std::cout << "m_grid " << m_grid_values.size() << " " << m_grid_values[0].size() << std::endl;

    int i_index = 0;
    for (iu = i0; iu.coord().u() < iend.coord().u(); iu.next_u()) {
        int j_index = 0;
        for (iv = iu; iv.coord().v() < iend.coord().v(); iv.next_v()) {

//            std::cout << iv.coord().u() << "/" << iend.coord().u() << " " << iv.coord().v() << "/" << iend.coord().v() << std::endl;

            for (iw = iv; iw.coord().w() <= iend.coord().w(); iw.next_w()) {
                if (iw.coord().w() == slice_index) {
                    SliceData data = SliceData(xmap[iw], iw.coord().u(), iw.coord().v(), iw.coord().w());
                    m_slice.push_back(data);
//                    std::cout << i_index << "/" << m_grid_values.size() << " " << j_index << "/" << m_grid_values[0].size() << std::endl;
                    m_grid_values[i_index][j_index] = xmap[iw];
                    break;
                }
            }
            j_index++;
        }
        i_index++;
    }

    std::cout << "Finished slicing..." << std::endl;
}

void Hasher::dump_slice(std::string file_name) {
    std::ofstream output_csv ;
    output_csv.open(file_name);

    output_csv << "u,v,data\n";
    for (SliceData data: m_slice) {
        output_csv << data.u() << "," << data.v() << "," << data.data() << "\n";
    }

    output_csv.close();
}


void Hasher::dump_slice(std::string file_name, std::vector<SliceData> data) {

    std::cout << "Dumping " << file_name << "\n";

    std::ofstream output_csv ;
    output_csv.open(file_name);

    output_csv << "u,v,data\n";
    for (SliceData data: data) {
        output_csv << data.u() << "," << data.v() << "," << data.data() << "\n";
    }

    output_csv.close();
}

float Hasher::gaussian_1d(float x,  float sigma) {
    float kernel = (1 / (sigma * sqrt(2 * M_PI))) * (exp((-(pow(x,2)))/(2*(pow(sigma,2)))));
    return kernel;
}

float Hasher::gaussian_2d(float x, float y,  float sigma) {
    float kernel = (1 / (sigma * sqrt(2 * M_PI))) * (exp((-((pow(x,2)+(pow(y,2)))))/(2*(pow(sigma,2)))));
    return kernel;
}


float Hasher::gaussian_3d(float x, float y, float z, float sigma) {
    float kernel = (1 / (sigma * 2 * sqrt(2 * M_PI))) * (exp((-((pow(x,2)+(pow(y,2))+(pow(z,2)))))/(2*(pow(sigma,2)))));
    return kernel;}



Matrix Hasher::generate_gaussian_kernel(int sigma) {

    std::cout << "Generating gaussian kernel with " << sigma << std::endl;
    int matrix_dimension = 1;
    Matrix kernel_matrix;

    int index_i = 0;
    for (int i = -matrix_dimension; i <= matrix_dimension; i++) {
        int index_j = 0;
        for (int j = -matrix_dimension; j <= matrix_dimension; j++) {
            kernel_matrix.m_matrix[index_i][index_j] = gaussian_2d(i, j, sigma);
            index_j++;
        }
        index_i++;
    }

    kernel_matrix.print();

    return kernel_matrix;
}

float Hasher::convolute(Matrix& kernel, Matrix& base) {

    float sum = 0.0f;

    for (int i = 0; i < base.m_matrix.size(); i++) {
        for (int j = 0; j < base.m_matrix[i].size(); j++) {
            sum = sum + (base.m_matrix[i][j] * kernel.m_matrix[i][j]);
        }
    }

    return sum;
}

std::vector<SliceData> Hasher::apply_gaussian(Matrix kernel) {

    auto m_grid_values_out = m_grid_values;
    std::vector<SliceData> m_grid_out;

    for (int i = 0; i < m_grid_values.size(); i++) {
        std::vector<float> grid = m_grid_values[i];
        for (int j = 0; j < grid.size(); j++) {

            int i_index = i;
            int j_index = j;

            Matrix grid;

            int i_iter_pos = 1;
            int i_iter_neg = -1;
            int j_iter_pos = 1;
            int j_iter_neg = -1;

            if (i == m_grid_values.size() - 1) {
                i_iter_pos = -(m_grid_values.size()-1);
            }

            if (j == m_grid_values[0].size() - 1) {
                j_iter_pos = -(m_grid_values[0].size()-1);
            }

            if (i == 0) {
                i_iter_neg = m_grid_values.size()-1;
            }

            if (j == 0) {
                j_iter_neg = m_grid_values[0].size()-1;
            }



            grid.m_matrix[0][0] = m_grid_values[i_index + i_iter_neg][j_index + j_iter_neg];
            grid.m_matrix[0][1] = m_grid_values[i_index + i_iter_neg][j_index];
            grid.m_matrix[0][2] = m_grid_values[i_index + i_iter_neg][j_index + j_iter_pos];

            grid.m_matrix[1][0] = m_grid_values[i_index][j_index + j_iter_neg];
            grid.m_matrix[1][1] = m_grid_values[i_index][j_index];
            grid.m_matrix[1][2] = m_grid_values[i_index][j_index + j_iter_pos];

            grid.m_matrix[2][0] = m_grid_values[i_index + i_iter_pos][j_index + j_iter_neg];
            grid.m_matrix[2][1] = m_grid_values[i_index + i_iter_pos][j_index];
            grid.m_matrix[2][2] = m_grid_values[i_index + i_iter_pos][j_index + j_iter_pos];

            float convolute_sum = convolute(kernel, grid);
            m_grid_values_out[i][j] = convolute_sum;
            m_grid_out.push_back(SliceData(convolute_sum, i, j, 0));

        }
    }

    return m_grid_out;
}

std::vector<SliceData> Hasher::difference_of_gaussian(std::vector<SliceData>& top, std::vector<SliceData>& bottom) {

    std::vector<SliceData> return_list;

    for (int i = 0; i < top.size(); i++) {

//        std::cout << i << std::endl;

        for (int j = 0; j < bottom.size(); j++) {

            SliceData top_slice = top[i];
            SliceData bottom_slice = bottom[i];

            if (top_slice.u() == bottom_slice.u() && top_slice.v() == bottom_slice.v()) {
                SliceData difference = SliceData(top_slice.data()-bottom_slice.data(), top_slice.u(), top_slice.v(), 1);
                return_list.push_back(difference);
                break;
            }
        }
    }

    return return_list;
}

int main() {

    std::cout << "Running" << std::endl;
    Hasher hash;

    auto s1 = std::chrono::high_resolution_clock::now();
    hash.load("./data/1hr2_final.mtz");
    auto e1 = std::chrono::high_resolution_clock::now();

    std::vector<SliceData> data;

    auto s2 = std::chrono::high_resolution_clock::now();
    hash.slice(3);
    auto e2 = std::chrono::high_resolution_clock::now();

    hash.dump_slice("./debug/slice_data.csv");

    auto s3 = std::chrono::high_resolution_clock::now();
    Matrix kernel = hash.generate_gaussian_kernel(3);
    auto e3 = std::chrono::high_resolution_clock::now();

    auto s4 = std::chrono::high_resolution_clock::now();
    std::vector<SliceData> blur_1 = hash.apply_gaussian(kernel);
    auto e4 = std::chrono::high_resolution_clock::now();

    hash.dump_slice("./debug/kernel1.csv", blur_1);

    Matrix kernel_2 = hash.generate_gaussian_kernel(5);
    std::vector<SliceData> blur_2 = hash.apply_gaussian(kernel_2);

    hash.dump_slice("./debug/kernel2.csv", blur_2);

    std::vector<SliceData> difference = hash.difference_of_gaussian(blur_1, blur_2);

    hash.dump_slice("./debug/kerneldifference.csv", difference);

    auto d1 = std::chrono::duration_cast<std::chrono::milliseconds>(e1-s1);
    auto d2 = std::chrono::duration_cast<std::chrono::milliseconds>(e2-s2);
    auto d3 = std::chrono::duration_cast<std::chrono::milliseconds>(e3-s3);
    auto d4 = std::chrono::duration_cast<std::chrono::milliseconds>(e4-s4);

    std::cout << "==Timings==\n" << "Loading " << d1.count() << "ms\n"
            << "Slicing " << d2.count() << "ms\n"
            << "Generating Kernel " << d3.count() << "ms\n"
            << "Applying Kernel " << d4.count() << "ms" << std::endl;


    return 0;
}