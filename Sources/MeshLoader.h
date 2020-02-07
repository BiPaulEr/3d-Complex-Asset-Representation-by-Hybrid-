#ifndef MESH_LOADER_H
#define MESH_LOADER_H

#include <string>
#include <memory>

#include "Mesh.h"

namespace MeshLoader {

	//Explaination of the function of each attribute is in Mesh.h

// Loads an OFF mesh file. See https://en.wikipedia.org/wiki/OFF_(file_format)
void loadOFF (const std::string & filename, std::shared_ptr<Mesh> meshPtr);

//Write List of indices of adjacent triangles in the file
void writeEdges(const std::string & filename, std::shared_ptr<Mesh> meshPtr);
// Load attribut triangle_Adj_Edges from mesh object
void loadEdges(const std::string & filename, std::shared_ptr<Mesh> meshPtr);

//Write in files the ListofTriangles and Vertex from DataInterpolation, size_radius used to analyse the mesh, the transfer data between thow sizes, the new position of vertex after the decimation
void writeScaleDataInterpolation(const std::string & filename, std::shared_ptr<Mesh> meshPtr, float &taille_radius);
//Load data from file to fill the attributes from mesh object : DataInterpolationVertexLoad/DataInterpolationTrianglesLoad/vertex_transferGroupe_load/vertexPositions_NEWLoad/size_radius
void loadScaleDataInterpolation(const std::string & filename, std::shared_ptr<Mesh> meshPtr);
//Show Data from m_vertexPositions_NEWLoad/m_DataInterpolationTrianglesLoad/m_DataInterpolationVertexLoad/m_vertex_transferGroupe_load
void showScaleDataInterpolation(std::shared_ptr<Mesh> meshPtr);

//Write  Potentiel and real macro surfaces current groupe object : Need to have component_Triangle_Potentiel_MacroSurface and component_Triangle_Real_MacroSurface initialise 
void writePotentielAndRealMacrosurface(const std::string & filename, std::shared_ptr<Mesh> meshPtr);
//Load  component_Triangle_Potentiel_MacroSurface_Load and component_Triangle_Real_MacroSurface_Load which is currentGroupe of vertex
void loadPotentielAndRealMacrosurface(const std::string & filename, std::shared_ptr<Mesh> meshPtr);


void writeDecimationAndVoxel(const std::string & filename, std::shared_ptr<Mesh> meshPtr);
//Load attributes vertexDecimationPositionsLoad/ triangleDecimationIndicesLoad/vertexVoxelPositionsLoad/triangleVoxelIndicesLoad from files to Mesh input
void loadDecimationAndVoxel(const std::string & filename, std::shared_ptr<Mesh> meshPtr);

//Write all attribute needed for interpolation/decmation/voxellisation in the files 
void WriteAllDataToFiles(std::string Name, std::shared_ptr<Mesh> meshPtr, float taille_radius);
//Load all attribute needed for interpolation/decmation/voxellisation from the files 
void LoadAllDataFromFiles(std::string Name, std::shared_ptr<Mesh> & meshPtr);
//Show to the console all the attributes load to check
void ShowAllDataLoaded(std::shared_ptr<Mesh> MeshPtr);
}


#endif // MESH_LOADER_H