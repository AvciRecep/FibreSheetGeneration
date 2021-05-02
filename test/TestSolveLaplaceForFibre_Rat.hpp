
#ifndef TESTSOLVELAPLACEFORFIBRE_HPP_
#define TESTSOLVELAPLACEFORFIBRE_HPP_

#include <cxxtest/TestSuite.h>

#include "UblasIncludes.hpp"
/* This is the class that is needed to solve a linear elliptic PDE. */
#include "SimpleLinearEllipticSolver.hpp"
/* This is needed to read mesh datafiles of the 'Triangles' format. */
#include "TrianglesMeshReader.hpp"
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

private:

    void ReadFilesIntoMap() //throw(Exception)
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
            cout << "There was a problem opening faces for reading " << endl;
        }
        //std::string line;
        if(!std::getline(inNode, line))
        {
            cout << "Error reading file line" << endl;
        }
        unsigned int numNodes, dummy1, dummy2, dummy3;
        stringstream numNodeLine(line);
        numNodeLine >> numNodes >> dummy1 >> dummy2 >> dummy3;
        nodeXYZ_st nodeStruct;
        while (numNodes > 0)
        {
            std::getline(inNode, line);
            stringstream nodeInfo(line);
            nodeInfo >> dummy1 >> nodeStruct.x >> nodeStruct.y >> nodeStruct.z;
            all_nodes.push_back(nodeStruct);
            numNodes -- ;
        }

        cout << "Number of nodes in mesh: " << all_nodes.size() << endl;
        /*
        nodeInfo_st nodeStructure;

        std::ifstream inCoordinate("projects/mesh/FibreSheetGeneration/rat_2_2_1.1.node");
        if (!inCoordinate)
        {
            cout << "There was a problem opening coordinates for reading " << endl;
        }

        // read ipxi file
        ifstream inElemDetails("projects/mesh/Stomach3D/hsb16_tesel_4.1.ipxi");
        if (!inElemDetails)
        {
            cout << "There was a problem opening element details for reading " << endl;
        }
        std::string lineEle;
        if(!std::getline(inCoordinate, line))
        {
            cout << "Error reading file line" << endl;
        }

        std::vector<int> interestedElem;

        // Top elements //
        interestedElem.push_back(1);
        interestedElem.push_back(2);
        interestedElem.push_back(3);
        interestedElem.push_back(4);
        interestedElem.push_back(5);
        interestedElem.push_back(6);
        interestedElem.push_back(7);
        interestedElem.push_back(8);

        // Bottom elements //
        interestedElem.push_back(57);
        interestedElem.push_back(58);
        interestedElem.push_back(59);
        interestedElem.push_back(60);
        interestedElem.push_back(61);
        interestedElem.push_back(62);
        interestedElem.push_back(63);
        interestedElem.push_back(64);

        unsigned int numNodes;
        stringstream numNodeLine(line);
        numNodeLine >> numNodes;
        while (numNodes > 0)
        {
            std::getline(inCoordinate, line);
            stringstream nodeCoor(line);
            std::getline(inElemDetails, lineEle);
            stringstream eleInfo(lineEle);
            nodeCoor >>nodeStructure.index >> nodeStructure.x >> nodeStructure.y >> nodeStructure.z;
            eleInfo >> dummy >> nodeStructure.cmEle >> nodeStructure.Xi1 >> nodeStructure.Xi2 >> nodeStructure.Xi3;
            if (find(interestedElem.begin(), interestedElem.end(),nodeStructure.cmEle)!= interestedElem.end())
                tetNodeInfo.push_back(nodeStructure);
            numNodes -- ;
        }
        cout << "Vector size -- " << tetNodeInfo.size() << endl;
        */
    }


    void sortDirchletAndNeumann() //throw(Exception)
    {

        for(std::vector<nodeInfo_st>::iterator itr = tetNodeInfo.begin(); itr != tetNodeInfo.end(); itr++)
        {
            nodeInfo_st myNodeInfo = *itr;
            unsigned int nodeIdx = myNodeInfo.index;
            if(face_node.find(nodeIdx) != face_node.end())
            {
                if(myNodeInfo.cmEle >= 1 && myNodeInfo.cmEle < 9 )
                {
                    if(myNodeInfo.Xi2 < 0.001)
                    {
                        nodeXYZ_st nodeSt;
                        nodeSt.x= myNodeInfo.x;
                        nodeSt.y= myNodeInfo.y;
                        nodeSt.z= myNodeInfo.z;
                        dir_bound_1.push_back(nodeSt);
                    }
                }
                else if (myNodeInfo.cmEle >= 57 && myNodeInfo.cmEle < 64 )
                {
                    if(myNodeInfo.Xi2 > 0.99)
                    {
                        nodeXYZ_st nodeSt;
                        nodeSt.x= myNodeInfo.x;
                        nodeSt.y= myNodeInfo.y;
                        nodeSt.z= myNodeInfo.z;
                        dir_bound_0.push_back(nodeSt);
                    }
                }
                else
                {
                	//do nothing
                }
            }
        }
        cout << "0 -- " << dir_bound_0.size() << endl;
        cout << "1 -- " << dir_bound_1.size() << endl;
    }
    void sortDirchletAndNeumann_rat() //throw(Exception)
    {

        for(std::vector<nodeXYZ_st>::iterator itr = all_nodes.begin(); itr != all_nodes.end(); itr++)
        {
            nodeXYZ_st myNodeInfo = *itr;
            //unsigned int nodeIdx = myNodeInfo.index;
            //if(face_node.find(nodeIdx) != face_node.end())
            //{
                double xdiff1, ydiff1, zdiff1, dist1;
                xdiff1 = myNodeInfo.x-0.0135;
                ydiff1 = myNodeInfo.y-0.9818;
                zdiff1 = myNodeInfo.z-(-4.08);
                dist1 = sqrt(xdiff1*xdiff1+ydiff1*ydiff1+zdiff1*zdiff1);
                //cout << dist1 << ',' << myNodeInfo.x << ',' << myNodeInfo.y << ',' << myNodeInfo.z << endl;
                if(dist1<0.2505)
                {
                        nodeXYZ_st nodeSt2;
                        nodeSt2.x= myNodeInfo.x;
                        nodeSt2.y= myNodeInfo.y;
                        nodeSt2.z= myNodeInfo.z;
                        dir_bound_0.push_back(nodeSt2);
                }

                double xdiff2, ydiff2, zdiff2, dist2;
                xdiff2 = myNodeInfo.x-(-0.0201);
                ydiff2 = myNodeInfo.y-1.348;
                zdiff2 = myNodeInfo.z-(-1.793);
                dist2 = sqrt(xdiff2*xdiff2+ydiff2*ydiff2+zdiff2*zdiff2);

                if(dist2<0.2505)
                {
                        nodeXYZ_st nodeSt1;
                        nodeSt1.x= myNodeInfo.x;
                        nodeSt1.y= myNodeInfo.y;
                        nodeSt1.z= myNodeInfo.z;
                        dir_bound_1.push_back(nodeSt1);
                }
            //}
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
        ReadFilesIntoMap();

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
                if (x < (myNode.x + 0.0001) && x > (myNode.x - 0.0001) && y < (myNode.y + 0.0001) && y > (myNode.y - 0.0001) && z < (myNode.z + 0.0001) && z > (myNode.z - 0.0001))
                {
                    bcc.AddDirichletBoundaryCondition(*iter, p_in_boundary_condition);
                    inCount++;
                }
            }

            for(std::vector<nodeXYZ_st>::iterator itr = dir_bound_0.begin(); itr != dir_bound_0.end(); itr++)
            {
                nodeXYZ_st myNode = *itr;
                if (x < (myNode.x + 0.0001) && x > (myNode.x - 0.0001) && y < (myNode.y + 0.0001) && y > (myNode.y - 0.0001) && z < (myNode.z + 0.0001) && z > (myNode.z - 0.0001))
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


        OutputFileHandler output_file_handler("TestLaplace_rat_16_16_1_longi");

        out_stream p_file = output_file_handler.OpenOutputFile("rat_16_16_1_linear_sol_longi.txt");

        PRINT_VARIABLE(result_repl.GetSize());


        // Loop over the entries of the solution. //
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            double z = mesh.GetNode(i)->rGetLocation()[2];
            double u = result_repl[i];

            (*p_file) << x << " " << y << " " << z << " " << u << "\n";
        }

        TRACE("Completed writing the linear solve values");

        out_stream p_file_grad = output_file_handler.OpenOutputFile("rat_16_16_1_grad_longi.txt");
        out_stream p_file_grad_mag = output_file_handler.OpenOutputFile("rat_16_16_1_mag_grad_longi.txt");
        std::vector<c_vector<double, 3u> > fibre_directions;
        c_vector<double, 3u> Node1, Node2, Node3, Node4;

        c_vector<double, 3> potVec, gradVec;
        c_matrix<double, 3, 3> element_jacobian, inverse_jacobian;
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

            gradVec =  prod(trans(inverse_jacobian), potVec);
            double magnitude = sqrt(gradVec[0]* gradVec[0]+ gradVec[1] * gradVec[1] + gradVec[2] * gradVec[2]);
            c_vector<double, 3u> fibre_direction;

            if (magnitude < 5 )
            {
                fibre_direction = fibre_directions[i-1];
                gradVec[0] = fibre_direction[0];
                gradVec[1] = fibre_direction[1];
                gradVec[2] = fibre_direction[2];
            }

            (*p_file_grad) << gradVec[0] << " " << gradVec[1] << " " << gradVec[2] << "\n";
            (*p_file_grad_mag) << magnitude  <<  "\n";

            fibre_direction[0] = gradVec[0];
            fibre_direction[1] = gradVec[1];
            fibre_direction[2] = gradVec[2];
            fibre_directions.push_back(fibre_direction);
        }

        VtkMeshWriter<3u, 3u> mesh_writer("TestLaplace_rat_16_16_1_longi", "mesh", false);
        mesh_writer.AddCellData("Fibre Direction", fibre_directions);
        mesh_writer.WriteFilesUsingMesh(mesh);

        PetscTools::Destroy(result);
    }

};

#endif
