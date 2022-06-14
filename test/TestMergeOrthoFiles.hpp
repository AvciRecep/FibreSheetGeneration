#include <cxxtest/TestSuite.h>
#include "UblasIncludes.hpp"
#include "TetrahedralMesh.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "OutputFileHandler.hpp"
#include "PetscSetupAndFinalize.hpp"
using namespace std;
#include "VtkMeshWriter.hpp"
#include "VtkMeshReader.hpp"
#include <../src/Eigen/Core>
#include <../src/Eigen/Geometry>
using Eigen::Vector3d;

std::string file_name = "rat_cm_32_32_8_lm_32_32_2.1";

class TestMergeOrthoFiles : public CxxTest::TestSuite
{
public:
    void TestReadOrthoMesh()
    {
      // Read mesh files
      std::string full_path = "projects/mesh/Stomach3D/" + file_name;
      TrianglesMeshReader<3,3> mesh_reader(full_path);
      // Now declare a tetrahedral mesh with the same dimensions... //
      TetrahedralMesh<3,3> mesh;
      // ... and construct the mesh using the mesh reader. //
      mesh.ConstructFromMeshReader(mesh_reader);

      // Output file
      full_path = "test_laplace_ortho_"+file_name+"_merged";
      OutputFileHandler output_file_handler(full_path);
      full_path = file_name+".ortho";
      out_stream p_file = output_file_handler.OpenOutputFile(full_path);

      // Read ortho files
      std::string ortho_file_lm = "/hpc/ravc486/Projects/CHASTE/chaste-hpc/chaste_ravc486/testoutput/test_laplace_ortho_"+file_name+".lm/"+file_name+".ortho";
      std::ifstream ortho_lm(ortho_file_lm);

      std::string ortho_file_cm = "/hpc/ravc486/Projects/CHASTE/chaste-hpc/chaste_ravc486/testoutput/test_laplace_ortho_"+file_name+".cm/"+file_name+".ortho";
      std::ifstream ortho_cm(ortho_file_cm);

      std::string line;
      std::getline(ortho_lm, line);
      std::getline(ortho_cm, line);
      unsigned int numElem;
      stringstream numElemLine(line);
      numElemLine >> numElem;
      (*p_file) << numElem << "\n";

      // Initialise variables
      Vector3d lm_normal;
      Vector3d lm_sheet;
      Vector3d lm_fibre;
      Vector3d cm_normal;
      Vector3d cm_sheet;
      Vector3d cm_fibre;
      std::string lm_line;
      std::string cm_line;
      std::vector<c_vector<double, 3u> > lm_fibre_dir;
      std::vector<c_vector<double, 3u> > lm_sheet_dir;
      std::vector<c_vector<double, 3u> > lm_normal_dir;
      std::vector<c_vector<double, 3u> > cm_fibre_dir;
      std::vector<c_vector<double, 3u> > cm_sheet_dir;
      std::vector<c_vector<double, 3u> > cm_normal_dir;
      double f1, f2, f3, s1, s2, s3, n1, n2, n3;
      // Mesh iterator
      for (TetrahedralMesh<3,3>::ElementIterator iter = mesh.GetElementIteratorBegin();
                      iter != mesh.GetElementIteratorEnd(); ++iter)
      {
          // Read Attributes
          double attribute = iter->GetAttribute();

          //
          std::getline(ortho_lm, lm_line);
          std::getline(ortho_cm, cm_line);
          //
          if (attribute == 1)
          {
            std::stringstream lineStream(lm_line);
            lineStream >> f1 >> f2 >> f3 >> s1 >> s2 >> s3 >> n1 >> n2 >> n3;
            lm_fibre(0) = f1;
            lm_fibre(1) = f2;
            lm_fibre(2) = f3;
            lm_sheet(0) = s1;
            lm_sheet(1) = s2;
            lm_sheet(2) = s3;
            lm_normal(0) = n1;
            lm_normal(1) = n2;
            lm_normal(2) = n3;
            cm_fibre(0) = 0;
            cm_fibre(1) = 0;
            cm_fibre(2) = 0;
            cm_sheet(0) = 0;
            cm_sheet(1) = 0;
            cm_sheet(2) = 0;
            cm_normal(0) = 0;
            cm_normal(1) = 0;
            cm_normal(2) = 0;
          }
          else if (attribute == 4)
          {
            std::stringstream lineStream(cm_line);
            lineStream >> f1 >> f2 >> f3 >> s1 >> s2 >> s3 >> n1 >> n2 >> n3;
            cm_fibre(0) = s1;
            cm_fibre(1) = s2;
            cm_fibre(2) = s3;
            cm_sheet(0) = f1;
            cm_sheet(1) = f2;
            cm_sheet(2) = f3;
            cm_normal(0) = n1;
            cm_normal(1) = n2;
            cm_normal(2) = n3;
            lm_fibre(0) = 0;
            lm_fibre(1) = 0;
            lm_fibre(2) = 0;
            lm_sheet(0) = 0;
            lm_sheet(1) = 0;
            lm_sheet(2) = 0;
            lm_normal(0) = 0;
            lm_normal(1) = 0;
            lm_normal(2) = 0;
          }
          else
          {
            std::stringstream lineStream(cm_line);
            lineStream >> f1 >> f2 >> f3 >> s1 >> s2 >> s3 >> n1 >> n2 >> n3;
            cm_fibre(0) = 0;
            cm_fibre(1) = 0;
            cm_fibre(2) = 0;
            cm_sheet(0) = 0;
            cm_sheet(1) = 0;
            cm_sheet(2) = 0;
            cm_normal(0) = 0;
            cm_normal(1) = 0;
            cm_normal(2) = 0;
            lm_fibre(0) = 0;
            lm_fibre(1) = 0;
            lm_fibre(2) = 0;
            lm_sheet(0) = 0;
            lm_sheet(1) = 0;
            lm_sheet(2) = 0;
            lm_normal(0) = 0;
            lm_normal(1) = 0;
            lm_normal(2) = 0;
          }
          // Write new ortho file
          (*p_file) << f1 << " " << f2 << " " << f3 << " "
                    << s1 << " " << s2 << " " << s3 << " "
                    << n1 << " " << n2 << " " << n3 << "\n";

          //
          c_vector<double, 3u> fibre_direction;
          fibre_direction[0] = lm_fibre(0);
          fibre_direction[1] = lm_fibre(1);
          fibre_direction[2] = lm_fibre(2);
          lm_fibre_dir.push_back(fibre_direction);
          fibre_direction[0] = cm_fibre(0);
          fibre_direction[1] = cm_fibre(1);
          fibre_direction[2] = cm_fibre(2);
          cm_fibre_dir.push_back(fibre_direction);

          c_vector<double, 3u> sheet_direction;
          sheet_direction[0] = lm_sheet(0);
          sheet_direction[1] = lm_sheet(1);
          sheet_direction[2] = lm_sheet(2);
          lm_sheet_dir.push_back(sheet_direction);
          sheet_direction[0] = cm_sheet(0);
          sheet_direction[1] = cm_sheet(1);
          sheet_direction[2] = cm_sheet(2);
          cm_sheet_dir.push_back(sheet_direction);

          c_vector<double, 3u> normal_direction;
          normal_direction[0] = lm_normal(0);
          normal_direction[1] = lm_normal(1);
          normal_direction[2] = lm_normal(2);
          lm_normal_dir.push_back(normal_direction);
          normal_direction[0] = cm_normal(0);
          normal_direction[1] = cm_normal(1);
          normal_direction[2] = cm_normal(2);
          cm_normal_dir.push_back(normal_direction);
      }
      p_file->close();

      // Write mesh with fibres as a vtk file
      full_path = "test_laplace_ortho_"+file_name+"_merged";
      VtkMeshWriter<3u, 3u> mesh_writer(full_path, "mesh", false);
      mesh_writer.AddCellData("LM Fibre", lm_fibre_dir);
      mesh_writer.AddCellData("LM Sheet", lm_sheet_dir);
      mesh_writer.AddCellData("LM Normal", lm_normal_dir);
      mesh_writer.AddCellData("CM Fibre", cm_fibre_dir);
      mesh_writer.AddCellData("CM Sheet", cm_sheet_dir);
      mesh_writer.AddCellData("CM Normal", cm_normal_dir);
      mesh_writer.WriteFilesUsingMesh(mesh);
    }
};
