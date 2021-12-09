
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

class MyPde : public AbstractLinearEllipticPde<3,3>
{
private:

public:

    MyPde()
    {
    }

    double ComputeConstantInUSourceTerm(const ChastePoint<3>& rX, Element<3,3>* pElement)
    {
        return 0.0;
    }

    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<3>& rX, Element<3,3>* pElement)
    {
        return 0.0;
    }

    c_matrix<double,3,3> ComputeDiffusionTerm(const ChastePoint<3>& rX)
    {
        return identity_matrix<double>(3);
    }
};


class TestSolvingLaplaceForFibre : public CxxTest::TestSuite
{

private:
    std::vector<nodeInfo_st> tetNodeInfo;
    std::vector<nodeXYZ_st> all_nodes;
    std::vector<nodeXYZ_st> dir_bound_0;
    std::vector<nodeXYZ_st> dir_bound_1;
    std::set<unsigned int> face_node;
    std::vector<nodeBoun_st> all_boun_Info;

private:

  void ReadFilesIntoMap_rat() //throw(Exception)
  {
        std::cout << "Read Files Into Map\n";
        // Read face file
        std::ifstream inFace("projects/mesh/FibreSheetGeneration/rat_16_16_1.1.face");
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
        std::ifstream inNode("projects/mesh/FibreSheetGeneration/rat_16_16_1.1.node");
        if (!inNode)
        {
            cout << "There was a problem opening nodes for reading " << endl;
        }
        //std::string line;
        if(!std::getline(inNode, line))
        {
            cout << "Error reading file line" << endl;
        }
        unsigned int numNodes, numNodesIter, dummy1, dummy2, dummy3;
        stringstream numNodeLine(line);
        numNodeLine >> numNodes >> dummy1 >> dummy2 >> dummy3;
        nodeXYZ_st nodeStruct;
        numNodesIter = numNodes;
        while (numNodesIter > 0)
        {
            std::getline(inNode, line);
            stringstream nodeInfo(line);
            nodeInfo >> dummy1 >> nodeStruct.x >> nodeStruct.y >> nodeStruct.z;
            all_nodes.push_back(nodeStruct);
            numNodesIter -- ;
        }

        cout << "Number of nodes in mesh: " << all_nodes.size() << endl;

        ifstream inBoun("projects/mesh/FibreSheetGeneration/rat_16_16_1.1.boun");
        if (!inBoun)
        {
            cout << "There was a problem opening boundary file for reading " << endl;
        }

        nodeBoun_st bounStruct;
        numNodesIter = numNodes;
        while (numNodesIter > 0)
        {
            std::getline(inBoun, line);
            stringstream bounInfo(line);
            bounInfo >> bounStruct.x >> bounStruct.y >> bounStruct.z >> bounStruct.long_boun >> bounStruct.circ_boun;
            all_boun_Info.push_back(bounStruct);
            numNodesIter--;
        }
        cout << "Vector size -- " << all_boun_Info.size() << endl;
    }

    void sortDirchletAndNeumann_rat()
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

        TrianglesMeshReader<3,3> mesh_reader("projects/mesh/FibreSheetGeneration/rat_16_16_1.1");
        // Now declare a tetrahedral mesh with the same dimensions... //
        TetrahedralMesh<3,3> mesh;
        // ... and construct the mesh using the mesh reader. //
        mesh.ConstructFromMeshReader(mesh_reader);

        // Next we instantiate an instance of our PDE we wish to solve. //
        MyPde pde;

        TRACE("Read files into map");
        ReadFilesIntoMap_rat();

        TRACE("Sort dirichilet boundaries");
        sortDirchletAndNeumann_rat();

        TRACE("Begin Fibre solve process");
        BoundaryConditionsContainer<3,3,1> bcc;

        ConstBoundaryCondition<3>* p_zero_boundary_condition = new ConstBoundaryCondition<3>(0.0);
        ConstBoundaryCondition<3>* p_in_boundary_condition = new ConstBoundaryCondition<3>(100);
        // We then get a boundary node iterator from the mesh... //
        TetrahedralMesh<3,3>::BoundaryNodeIterator iter = mesh.GetBoundaryNodeIteratorBegin();
        // ...and loop over the boundary nodes, getting the x and y values. //
        unsigned int inCount = 0;
        unsigned int outCount = 0;
        while (iter < mesh.GetBoundaryNodeIteratorEnd())
        {
            double x = (*iter)->GetPoint()[0];
            double y = (*iter)->GetPoint()[1];
            double z = (*iter)->GetPoint()[2];
            // If x=0 or y=0... //
            for(std::vector<nodeXYZ_st>::iterator itr = dir_bound_1.begin(); itr != dir_bound_1.end(); itr++)
            {
                nodeXYZ_st myNode = *itr;
                if (x < (myNode.x + 0.000001) && x > (myNode.x - 0.000001) && y < (myNode.y + 0.000001) && y > (myNode.y - 0.000001) && z < (myNode.z + 0.000001) && z > (myNode.z - 0.000001))
                {
                    bcc.AddDirichletBoundaryCondition(*iter, p_in_boundary_condition);
                    inCount++;
                }
            }

            for(std::vector<nodeXYZ_st>::iterator itr = dir_bound_0.begin(); itr != dir_bound_0.end(); itr++)
            {
                nodeXYZ_st myNode = *itr;
                if (x < (myNode.x + 0.000001) && x > (myNode.x - 0.000001) && y < (myNode.y + 0.000001) && y > (myNode.y - 0.000001) && z < (myNode.z + 0.000001) && z > (myNode.z - 0.000001))
                {
                    bcc.AddDirichletBoundaryCondition(*iter, p_zero_boundary_condition);
                    outCount++;
                }
            }
            iter++;
        }
        cout << "Compared and found IN: " << inCount << endl;
        cout << "Compared and found OUT: " << outCount << endl;
        SimpleLinearEllipticSolver<3,3> solver(&mesh, &pde, &bcc);

        // To solve, just call {{{Solve()}}}. A PETSc vector is returned. //
        Vec result = solver.Solve();

        ReplicatableVector result_repl(result);

        OutputFileHandler output_file_handler("TestLaplace_rat_16_16_1_longi_Dec2021");

        out_stream p_file = output_file_handler.OpenOutputFile("rat_16_16_1_linear_sol_longi.txt");

        PRINT_VARIABLE(result_repl.GetSize());


        // Loop over the entries of the solution. //
        for (unsigned i=0; i<result_repl.GetSize(); i++) //result_repl.GetSize()
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            double z = mesh.GetNode(i)->rGetLocation()[2];
            double u = result_repl[i];

            (*p_file) << x << " " << y << " " << z << " " << u << "\n";
            //cout << i << " " << x << " " << y << " " << z << " " << u << "\n";
        }

        TRACE("Completed writing the linear solve values");

        out_stream p_file_grad = output_file_handler.OpenOutputFile("rat_16_16_1_grad_longi.txt");
        out_stream p_file_grad_mag = output_file_handler.OpenOutputFile("rat_16_16_1_mag_grad_longi.txt");
        std::vector<c_vector<double,3u> > fibre_directions;
        c_vector<double,3u> Node1, Node2, Node3, Node4;
        c_vector<double,3> potVec, gradVec;
        c_matrix<double,3,3> element_jacobian, inverse_jacobian;
        double dummy;
        for(unsigned i = 0; i < mesh.GetNumElements(); i++)
        {
            double L1 = result_repl[mesh.GetElement(i)->GetNodeGlobalIndex(0)];
            double L2 = result_repl[mesh.GetElement(i)->GetNodeGlobalIndex(1)];
            double L3 = result_repl[mesh.GetElement(i)->GetNodeGlobalIndex(2)];
            double L4 = result_repl[mesh.GetElement(i)->GetNodeGlobalIndex(3)];

            mesh.GetElement(i)->CalculateInverseJacobian(element_jacobian,
                                  dummy,inverse_jacobian);

            potVec[0] = L2-L1;
            potVec[1] = L3-L1;
            potVec[2] = L4-L1;

            gradVec = prod(trans(inverse_jacobian), potVec);
            double magnitude = sqrt(gradVec[0]* gradVec[0] + gradVec[1] * gradVec[1] + gradVec[2] * gradVec[2]);
            c_vector<double,3u> fibre_direction;
/*
            if (magnitude < 5)
            {
                fibre_direction = fibre_directions[i-1];
                gradVec[0] = fibre_direction[0];
                gradVec[1] = fibre_direction[1];
                gradVec[2] = fibre_direction[2];
            }
*/
            (*p_file_grad) << gradVec[0] << " " << gradVec[1] << " " << gradVec[2] << "\n";
            (*p_file_grad_mag) << magnitude  <<  "\n";

            fibre_direction[0] = gradVec[0];
            fibre_direction[1] = gradVec[1];
            fibre_direction[2] = gradVec[2];
            fibre_directions.push_back(fibre_direction);
        }

        VtkMeshWriter<3u, 3u> mesh_writer("TestLaplace_rat_16_16_1_longi_Dec2021", "mesh", false);
        mesh_writer.AddCellData("Fibre Direction", fibre_directions);
        mesh_writer.WriteFilesUsingMesh(mesh);

        PetscTools::Destroy(result);
    }

};
#endif
