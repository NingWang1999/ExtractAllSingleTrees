#include <iomanip>
 // PCL
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/console/parse.h>
// MakeDataset
#include "stem_map.h"
#include <filesystem>

using namespace pcl::console;

int
main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <raw_data_folder> <rough_data_folder>\n";
        return -1;
    }
#if defined(_OPENMP)
    print_info("[PARALLEL PROCESSING USING ");
    print_value("%d", omp_get_max_threads());
    print_info(" THREADS] \n\n");
#else
    print_info("[NON-PARALLEL PROCESSING] \n\n");
#endif

    std::string raw_data_directory = argv[1];   // e.g., "raw_data"
    std::string rough_data_directory = argv[2]; // e.g., "rough_data"

    double tic, toc, time_val = 0.0;

    for (const auto& entry : std::filesystem::directory_iterator(raw_data_directory)) {
        //============================= Load Data =============================
        if (entry.path().extension() == ".pcd") {  // only process .pcd file
            std::string filename = entry.path().filename().string();   // e.g., "View_1.pcd"
            std::string base_name = filename.substr(0, filename.find_last_of(".")); // e.g., "View_1"

            std::cout << "Processing file: " << filename << std::endl;

            Cloud3D::Ptr cloud_src(new Cloud3D);
            tic = omp_get_wtime();
            if (pcl::io::loadPCDFile<Point3D>(entry.path().string(), *cloud_src) == -1) {
                std::cerr << "Failed to load file: " << filename << std::endl;
                continue;
            }
            toc = omp_get_wtime();

            print_info("  Loaded ");
            print_value("%d", cloud_src->size());
            print_info(" points from ");
            print_value("%s", filename.c_str());
            print_info(" in ");
            print_value("%f", toc - tic);
            print_info(" s.\n");

            //======================= Extract All Single Trees =======================
            std::cout << "Extracting trees from " << filename << "..." << std::endl;
            tic = omp_get_wtime();

            std::string view_output_dir = rough_data_directory + "/" + base_name; // e.g., "rough_data/View_1"
           
            {
                Mapping mapping;
                mapping.setInputCloud(cloud_src->makeShared());

                /*mapping.extract(view_output_dir, base_name);*/
                try {
                    mapping.extract(view_output_dir, base_name);
                }
                catch (const std::exception& e) {
                    std::cerr << "Error in mapping.extract(): " << e.what() << std::endl;
                }
            }
            toc = omp_get_wtime();
            time_val += (toc - tic);
            print_info("Extracted trees from ");
            print_value("%s", filename.c_str());
            print_info(" in ");
            print_value("%f", toc - tic);
            print_info(" s.\n\n");
        }
    }
    return 0;
}