#include <sys/stat.h>
#include <stdio.h>
#include <stdarg.h>
#include <getopt.h>
#include <errno.h>

#include <map>


#include <vtkSmartPointer.h>
#include <vtkExtractEdges.h>
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"

#include "H5hut.h"

#if !defined (PARALLEL_IO)
#define MPI_Init(argc, argv)
#define MPI_Comm_size(comm, nprocs) { *nprocs = 1; }
#define MPI_Comm_rank(comm, myproc) { *myproc = 0; }
#define MPI_Finalize()
#define MPI_COMM_WORLD (0)
#endif

const char* version = "0.1.0";
int convert_boundary = 1;
int convert_volume = 0;
double x_shift = 0.0;
double y_shift = 0.0;
double z_shift = 0.0;

const struct option longopts[] = {
        {"version",     no_argument, 0,                 'v'},
        {"help",        no_argument, 0,                 'h'},
        {"boundary",    no_argument, &convert_boundary, 1},
        {"volume",      no_argument, &convert_volume,   1},
        {"no-boundary", no_argument, &convert_boundary, 0},
        {"no-volume",   no_argument, &convert_volume,   0},
        {"shift",       required_argument,  0,          's'},
        {0,0,0,0},
};

typedef std::map<int, size_t> IdMap;

static void
usage (
        char* cmd
        ) {
        std::cout << std::endl << "Usage: " << std::endl;
        std::cout << "  " << cmd << " [OPTIONS] FILE..." << std::endl << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "  --(no-)boundary   do (not) convert boundary triangles." << std::endl;
        std::cout << "                    default is yes" << std::endl;
        std::cout << "  --(no-)volume     do (not) convert volume mesh." << std::endl;
        std::cout << "                    default is no" << std::endl;
        std::cout << std::endl;
        exit (1);
}

static void
print_version (
        char* cmd
        ) {
        std::cout << "Version: " << cmd << " " << version << std::endl;
        exit (1);
}

static int
init (
        int argc,
        char* argv[]
        ) {
        int index;
        int iarg = 0;

        //turn off getopt error message
        opterr = 1; 

        while(iarg != -1) {
                iarg = getopt_long(argc, argv, "vh?", longopts, &index);

                switch (iarg) {
                case 'h':
                case '?':
                        usage (argv[0]);
                        break;

                case 'v':
                        print_version (argv[0]);
                        break;
                case 's':
                        sscanf (optarg, "%lf,%lf,%lf", &x_shift, &y_shift, &z_shift);
                        cout << "shift = (" << x_shift << ", " << y_shift << ", " << z_shift << ")" << endl;
                        break;
                }
        }
        return argc - optind;
}

void
convert_vtk2h5grid (
        vtkUnstructuredGrid* vtk_grid,
        h5t_mesh_t* h5_grid,
        int cell_type
        ) {
        IdMap idmap;
        vtkIdType num_vtk_cells = vtk_grid->GetNumberOfCells ();
        vtkIdType num_vtk_pts = vtk_grid->GetNumberOfPoints ();
        h5_loc_idx_t (*cells)[4] = new h5_loc_idx_t[num_vtk_cells][4];
        int h5_cell_idx = 0;
        int h5_vertex_idx = 0;
        H5FedBeginStoreVertices (h5_grid, num_vtk_pts);
        for (vtkIdType vtk_cell_id = 0; vtk_cell_id < num_vtk_cells; vtk_cell_id++) {
                if (vtk_grid->GetCellType (vtk_cell_id) != cell_type) 
                        continue;
                
                vtkIdType num_pts;
                vtkIdType* pts;
                vtk_grid->GetCellPoints (vtk_cell_id, num_pts, pts);
                for (vtkIdType i = 0; i < num_pts; i++) {
                        IdMap::iterator it = idmap.find (pts[i]);
                        if (it == idmap.end ()) {
                                // add point to H5hut mesh
                                double pt[3];
                                vtk_grid->GetPoint (pts[i], pt);
                                pt[0] += x_shift;
                                pt[1] += y_shift;
                                pt[2] += z_shift;
                                H5FedStoreVertex (h5_grid, -1, pt);
                                // map pt index in vtk file to pt index in H5hut file
                                idmap.insert (IdMap::value_type (pts[i], h5_vertex_idx));
                                cells[h5_cell_idx][i] = h5_vertex_idx++;
                        } else {
                                cells[h5_cell_idx][i] = it->second;
                        }
                }
                h5_cell_idx++;
        }
        H5FedEndStoreVertices (h5_grid);
        int num_h5_vertices = h5_vertex_idx;
        int num_h5_cells = h5_cell_idx;
        cout << "  number of points in mesh: " << num_h5_vertices << endl;
        cout << "  number of cells in mesh:  " << num_h5_cells << endl;
        // add cells to H5hut file
        H5FedBeginStoreElements (h5_grid, num_h5_cells);
        for (int i = 0; i < num_h5_cells; i++) {
                H5FedStoreElement (h5_grid, cells[i]);
        }
        H5FedEndStoreElements (h5_grid);
        delete cells;
}

h5_err_t
H5Errorhandler (
	const char* fmt,
	va_list ap
	) {
	vfprintf (stderr, fmt, ap);
        exit (42);
	return -2;
}

int
main (
        int argc,
        char* argv[]
        ) {
        MPI_Init (&argc, &argv);
        MPI_Comm comm = MPI_COMM_WORLD;
        int comm_size;
        MPI_Comm_size (comm,&comm_size);

        argc = init (argc, argv);
        if (argc == 0) {
                std::cout << "No file(s) to convert?!" << std::endl;
                usage (argv[0]);
        }
        argv += optind;
        H5SetErrorHandler (H5Errorhandler);
        MPI_Init (&argc, &argv);
        for (int i = 0; i < argc; i++) {
                std::string vtk_filename = argv[i];
                std::string h5grid_filename = vtk_filename.substr (0, vtk_filename.find_last_of ('.')) + ".h5";
                struct stat st;
                if (stat (vtk_filename.c_str(), &st) < 0) {
                        perror (vtk_filename.c_str());
                        exit (1);
                }
                if (!(st.st_mode & S_IRUSR)) {
                        std::cerr << vtk_filename << ": not readable" << std::endl;
                        exit (1);
                }
                if (stat (h5grid_filename.c_str(), &st) == 0) {
                        std::cerr << h5grid_filename << ": already exits" << std::endl;
                        exit (1);
                }
                if (errno != ENOENT) {
                        perror (h5grid_filename.c_str());
                        exit (1);
                }
                std::cout << vtk_filename << ": converting ...." << std::endl;

                // read vtk file
                vtkSmartPointer<vtkUnstructuredGridReader> reader =
                        vtkSmartPointer<vtkUnstructuredGridReader>::New ();
                reader->SetFileName (vtk_filename.c_str ());
                reader->Update ();
                vtkUnstructuredGrid* vtk_grid = reader->GetOutput ();

                // open new H5hut file
                h5_file_t f = H5OpenFile (h5grid_filename.c_str(), H5_O_WRONLY, 0);

                h5t_mesh_t* h5_grid;
                if (convert_boundary) {
                        // add boundary mesh
                        H5FedAddTriangleMesh (f, "0", &h5_grid);
                        convert_vtk2h5grid (vtk_grid, h5_grid, VTK_TRIANGLE);
                        H5FedCloseMesh (h5_grid);
                }
                if (convert_volume) {
                        // add volume mesh
                        H5FedAddTetrahedralMesh (f, "0", &h5_grid);
                        convert_vtk2h5grid (vtk_grid, h5_grid, VTK_TETRA);
                        H5FedCloseMesh (h5_grid);
                }
                H5CloseFile (f);
        }
        MPI_Finalize ();
        return 0;
}
