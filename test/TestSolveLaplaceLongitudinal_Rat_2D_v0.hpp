#ifndef TESTSOLVELAPLACEFORFIBRE_HPP_
#define TESTSOLVELAPLACEFORFIBRE_HPP_

#include <cxxtest/TestSuite.h>

#include "UblasIncludes.hpp"
/* This is the class that is needed to solve a linear elliptic PDE. */
#include "SimpleLinearEllipticSolver.hpp"
/* This is needed to read mesh datafiles of the 'Triangles' format. */
//#include "TrianglesMeshReader.hpp"
/* This class represents the mesh internally. */
#include "TetrahedralMesh.hpp"
/* These are used to specify boundary conditions for the PDEs. */
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
/* This class helps us deal with output files. */
#include "OutputFileHandler.hpp"
/* The following header must be included in every test that uses PETSc. Note that it
 * cannot be included in the source code. */
#include "PetscSetupAndFinalize.hpp"

#include "VtkMeshWriter.hpp"

// My includes
#include <math.h>
#include <fstream>
#include <vector>
#include <string>
#include "Debug.hpp"

#define PI 3.14159265

using namespace std;

struct nodeXYZ_st
{
    double x;
    double y;
    double z;
};
struct nodeInfo_st
{
    unsigned int index;
    double x;
    double y;
    double z;
    unsigned int cmEle;
    double Xi1;
    double Xi2;
    double Xi3;
};
struct nodeBoun_st
{
    unsigned int index;
    double x;
    double y;
    double z;
    double long_boun;
    double circ_boun;
};
string file_name = "rat_scaffold_16_16_2_2D";
string full_path = "";

class MyPde : public AbstractLinearEllipticPde<2,3>
{
private:

public:

    MyPde()
    {
    }

    double ComputeConstantInUSourceTerm(const ChastePoint<3>& rX, Element<2,3>* pElement)
    {
        return 0.0;
    }

    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<3>& rX, Element<2,3>* pElement)
    {
        return 0.0;
    }

    c_matrix<double,3,3> ComputeDiffusionTerm(const ChastePoint<3>& rX)
    {
        return identity_matrix<double>(3);
    }
};


class TestSolveLaplaceLongitudinal : public CxxTest::TestSuite
{

private:
    std::vector<nodeInfo_st> tetNodeInfo;
    std::vector<nodeXYZ_st> all_nodes;
    std::vector<nodeXYZ_st> dir_bound_0;
    std::vector<nodeXYZ_st> dir_bound_1;
    std::set<unsigned int> face_node;
    std::vector<nodeBoun_st> all_boun_Info;

private:
  void ReadFilesIntoMap_Rat()
  {
        std::cout << "Read Files Into Map\n";
        // Read face file
        full_path = "projects/mesh/Stomach3D/" + file_name + ".face";
        std::ifstream inFace(full_path);
        if (!inFace)
        {
            cout << "There was a problem opening faces for reading " << endl;
        }
        std::string line;
        if(!std::getline(inFace, line))
        {
            cout << "Error reading file line" << endl;
        }
        unsigned int numFaces, dummy;
        stringstream numFaceLine(line);
        numFaceLine >> numFaces >> dummy;
        while (numFaces > 0)
        {
            unsigned int temp1, temp2, temp3;
            std::getline(inFace, line);
            stringstream faceInfo(line);
            faceInfo >> dummy >> temp1 >> temp2 >> temp3;
            face_node.insert(temp1);
            face_node.insert(temp2);
            face_node.insert(temp3);
            numFaces -- ;
        }

        cout << "Number of nodes in face: " << face_node.size() << endl;

        // Read node file
        full_path = "projects/mesh/Stomach3D/" + file_name + ".node";
        std::ifstream inNode(full_path);
        if (!inNode)
        {
            cout << "There was a problem opening nodes for reading " << endl;
        }
        //std::string line;
        if(!std::getline(inNode, line))
        {
            cout << "Error reading file line" << endl;
        }
        unsigned int numNodes, numNodesTmp, dummy1, dummy2, dummy3;
        stringstream numNodeLine(line);
        numNodeLine >> numNodes >> dummy1 >> dummy2 >> dummy3;
        nodeXYZ_st nodeStruct;
        numNodesTmp = numNodes;
        while (numNodesTmp > 0)
        {
            std::getline(inNode, line);
            stringstream nodeInfo(line);
            nodeInfo >> dummy1 >> nodeStruct.x >> nodeStruct.y >> nodeStruct.z;
            all_nodes.push_back(nodeStruct);
            numNodesTmp -- ;
        }

        cout << "Number of nodes in mesh: " << all_nodes.size() << endl;

        // Read boundary condition file
        full_path = "projects/mesh/Stomach3D/" + file_name + ".lm.boun";
        ifstream inBoun(full_path);
        if (!inBoun)
        {
            cout << "There was a problem opening boundary file for reading " << endl;
        }
        nodeBoun_st bounStruct;
        numNodesTmp = numNodes;
        while (numNodesTmp > 0)
        {
            std::getline(inBoun, line);
            stringstream bounInfo(line);
            bounInfo >> bounStruct.x >> bounStruct.y >> bounStruct.z >> bounStruct.long_boun >> bounStruct.circ_boun;
            all_boun_Info.push_back(bounStruct);
            numNodesTmp--;
        }
        cout << "Vector size -- " << all_boun_Info.size() << endl;
    }

    void sortDirchletAndNeumann_Rat() //throw(Exception)
    {

        for(std::vector<nodeBoun_st>::iterator itr = all_boun_Info.begin(); itr != all_boun_Info.end(); itr++)
        {
            nodeBoun_st myNodeInfo = *itr;
            if (myNodeInfo.long_boun == 1)
            {
              nodeXYZ_st nodeSt;
              nodeSt.x= myNodeInfo.x;
              nodeSt.y= myNodeInfo.y;
              nodeSt.z= myNodeInfo.z;
              dir_bound_1.push_back(nodeSt);
            }
            else if (myNodeInfo.long_boun == 0)
            {
              nodeXYZ_st nodeSt;
              nodeSt.x= myNodeInfo.x;
              nodeSt.y= myNodeInfo.y;
              nodeSt.z= myNodeInfo.z;
              dir_bound_0.push_back(nodeSt);
            }
        }
        cout << "0 -- " << dir_bound_0.size() << endl;
        cout << "1 -- " << dir_bound_1.size() << endl;
    }

public:
    void TestSolvingFibre() //throw(Exception)
    {
        ///// Read mesh files
        full_path = "projects/mesh/Stomach2D/" + file_name;
        TrianglesMeshReader<2,3> mesh_reader(full_path);
        TetrahedralMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        ///// Next we instantiate an instance of our PDE we wish to solve.
        MyPde pde;

        ///// Define boundary conditions
        BoundaryConditionsContainer<2,3,1> bcc;
        ConstBoundaryCondition<3>* zero_cond = new ConstBoundaryCondition<3>(0.0);
        ConstBoundaryCondition<3>* hundred_cond = new ConstBoundaryCondition<3>(100);
        TetrahedralMesh<2,3>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        int n = 0;
        int n0 = 0;
        int n100 = 0;
        while (iter < mesh.GetBoundaryNodeIteratorEnd())
        {
            double x = (*iter)->GetPoint()[0];
            double y = (*iter)->GetPoint()[1];
            double z = (*iter)->GetPoint()[2];
            // If x=0 or y=0... //
            if (x*x+(y-1)*(y-1)+(z+4)*(z+4) < 0.1)
            {
                bcc.AddDirichletBoundaryCondition(*iter, hundred_cond);
                n100++;
            }
            if ((x+0.02)*(x+0.02)+(y-1.35)*(y-1.35)+(z+1.79)*(z+1.79) < 0.1)
            {
                bcc.AddDirichletBoundaryCondition(*iter, zero_cond);
                n0++;
            }
            n++;
            iter++;
        }

        ///// Solve
        SimpleLinearEllipticSolver<2,3> solver(&mesh, &pde, &bcc);
        Vec result = solver.Solve();
        ReplicatableVector result_repl(result);

        ///// Write the solution
        full_path = "test_laplace_longi_" + file_name;
        OutputFileHandler output_file_handler(full_path);
        full_path = file_name+"_laplace_longi.txt";
        out_stream p_file = output_file_handler.OpenOutputFile(full_path);

        for (unsigned i=0; i<result_repl.GetSize(); i++) //result_repl.GetSize()
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            double z = mesh.GetNode(i)->rGetLocation()[2];
            double u = result_repl[i];

            (*p_file) << x << " " << y << " " << z << " " << u << "\n";
            cout << i << " " << x << " " << y << " " << z << " " << u << "\n";
        }

        TRACE("Completed writing the linear solve values");
        full_path = file_name+"_laplace_longi_grad.txt";
        out_stream p_file_grad = output_file_handler.OpenOutputFile(full_path);
        full_path = file_name+"_laplace_longi_grad_mag.txt";
        out_stream p_file_grad_mag = output_file_handler.OpenOutputFile(full_path);
        std::vector<c_vector<double,3u> > fibre_directions;
        c_vector<double,3u> Node1, Node2, Node3;
        c_vector<double,3> potVec, gradVec;
        c_matrix<double,3,2> element_jacobian;
        c_matrix<double,2,3> inverse_jacobian;
        double determinant;
        for(unsigned i = 0; i < mesh.GetNumElements(); i++)
        {
            double L1 = result_repl[mesh.GetElement(i)->GetNodeGlobalIndex(0)];
            double L2 = result_repl[mesh.GetElement(i)->GetNodeGlobalIndex(1)];
            double L3 = result_repl[mesh.GetElement(i)->GetNodeGlobalIndex(2)];
            mesh.GetElement(i)->CalculateInverseJacobian(element_jacobian,
                                  determinant,inverse_jacobian);

            potVec[0] = L2-L1;
            potVec[1] = L3-L1;
            gradVec = prod(trans(inverse_jacobian), potVec);
            double magnitude = sqrt(gradVec[0]* gradVec[0] + gradVec[1] * gradVec[1] + gradVec[2] * gradVec[2]);
            c_vector<double,3u> fibre_direction;
            if (magnitude < 0.5)
            {
                fibre_direction = fibre_directions[i-1];
                gradVec[0] = fibre_direction[0];
                gradVec[1] = fibre_direction[1];
                gradVec[2] = fibre_direction[2];
            }

            (*p_file_grad) << gradVec[0] << " " << gradVec[1] << " " << gradVec[2] << "\n";
            (*p_file_grad_mag) << magnitude << "\n";

            fibre_direction[0] = gradVec[0];
            fibre_direction[1] = gradVec[1];
            fibre_direction[2] = gradVec[2];
            fibre_directions.push_back(fibre_direction);
/*
            if(i==10000000000)
            {
              TRACE("BEGIN\n");
              cout << "Element ID: " << i << " L1: " << L1 << " L2: " << L2 << " L3: " << L3 << " L4: " << L4 << "\n";
              cout << "Node 0 (ID, x, y, z): " << mesh.GetElement(i)->GetNodeGlobalIndex(0) << " " << mesh.GetNode(mesh.GetElement(i)->GetNodeGlobalIndex(0))->rGetLocation()[0]<< " " << mesh.GetNode(mesh.GetElement(i)->GetNodeGlobalIndex(0))->rGetLocation()[1]<< " " << mesh.GetNode(mesh.GetElement(i)->GetNodeGlobalIndex(0))->rGetLocation()[2]<< "\n";
              cout << "Node 1 (ID, x, y, z): " << mesh.GetElement(i)->GetNodeGlobalIndex(1) << " " << mesh.GetNode(mesh.GetElement(i)->GetNodeGlobalIndex(1))->rGetLocation()[0]<< " " << mesh.GetNode(mesh.GetElement(i)->GetNodeGlobalIndex(1))->rGetLocation()[1]<< " " << mesh.GetNode(mesh.GetElement(i)->GetNodeGlobalIndex(1))->rGetLocation()[2]<< "\n";
              cout << "Node 2 (ID, x, y, z): " << mesh.GetElement(i)->GetNodeGlobalIndex(2) << " " << mesh.GetNode(mesh.GetElement(i)->GetNodeGlobalIndex(2))->rGetLocation()[0]<< " " << mesh.GetNode(mesh.GetElement(i)->GetNodeGlobalIndex(2))->rGetLocation()[1]<< " " << mesh.GetNode(mesh.GetElement(i)->GetNodeGlobalIndex(2))->rGetLocation()[2]<< "\n";
              cout << "Node 2 (ID, x, y, z): " << mesh.GetElement(i)->GetNodeGlobalIndex(3) << " " << mesh.GetNode(mesh.GetElement(i)->GetNodeGlobalIndex(3))->rGetLocation()[0]<< " " << mesh.GetNode(mesh.GetElement(i)->GetNodeGlobalIndex(3))->rGetLocation()[1]<< " " << mesh.GetNode(mesh.GetElement(i)->GetNodeGlobalIndex(3))->rGetLocation()[2]<< "\n";
              cout << "GradVec: " << gradVec[0] << " " << gradVec[1] << " " << gradVec[2] << "\n";
              TRACE("END");
            }*/
        }
        full_path = "test_laplace_longi_" + file_name;
        VtkMeshWriter<2u, 3u> mesh_writer(full_path, "mesh", false);
        mesh_writer.AddCellData("Fibre Direction", fibre_directions);
        mesh_writer.WriteFilesUsingMesh(mesh);

        PetscTools::Destroy(result);
    }
};
#endif
