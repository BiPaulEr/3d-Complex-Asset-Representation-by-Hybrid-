// ----------------------------------------------
// Base code for practical computer graphics
// assignments.
//
// Copyright (C) 2018 Tamy Boubekeur
// All rights reserved.
// ----------------------------------------------

#define _USE_MATH_DEFINES
#define VOXELIZER_IMPLEMENTATION

#include <glad/glad.h>

#include <cstdlib>
#include <cstdio>
#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>
#include <exception>


#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <glm/gtc/quaternion.hpp>

#include "Error.h"
#include "ShaderProgram.h"
#include "Camera.h"
#include "Mesh.h"
#include "MeshLoader.h"
#include "Material.h"
#include "LightSource.h"
#include "mdmeshDecimator.h"
#include "voxelizer.h"

#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

static const std::string SHADER_PATH ("Resources/Shaders/");

static const std::string DEFAULT_MESH_FILENAME("Resources/mesh_collection/robot.off");//Models/sphere_2.off");//mesh_collection/max_50K.off");
std::string MESH_NAME("robot");

using namespace std;

// Window parameters
static GLFWwindow * windowPtr = nullptr;
ImGuiIO io;
ImVec2 test;
ImVec2 test2;
static char * glwfversion = "#version 450";
// Pointer to the current camera model
static std::shared_ptr<Camera> cameraPtr;

static std::vector<std::string> ListNameMesh = { "robot","homer","pegaso","dancer2","bozbezbozzel" };
static std::vector<std::string> ListNameObjectToLoad = { "_edges.off","_PotentielAndRealSurfaces.off","_DecimationAndVoxel.off","_DataInterpolation.off" };
static int oldMEshLoaded = 0;
static int newMEshLoaded = 0;
static bool my_tool_active = false;
// Pointer to the displayed mesh
static std::shared_ptr<Mesh> meshPtr;

// Pointer to GPU shader pipeline i.e., set of shaders structured in a GPU program
static std::shared_ptr<ShaderProgram> shaderProgramPtr; // A GPU program contains at least a vertex shader and a fragment shader

// Camera control variables
static float meshScale = 1.0; // To update based on the mesh size, so that navigation runs at scale
glm::vec3 center;
static bool isRotating (false);
static bool isPanning (false);
static bool isZooming (false);
static double baseX (0.0), baseY (0.0);
static glm::vec3 baseTrans (0.0);
static glm::vec3 baseRot (0.0);
float taille_radius = 0.07f;
static float Scale=0.f;
static float interpolate = 0;

static float color[4] = { 0.5f,0.5f,0.5f,0.5f };
static int Mode_ = 0;
// Decimation need
pair<int, std::vector<string>> m_pairEdgesToDecime;
std::vector<glm::vec3> m_vertexDecimationPositions;
std::vector<glm::uvec3> m_triangleDecimationIndices;
pair< std::vector<glm::vec3>, std::vector<glm::uvec3>> mDecimationPairPositionIndices;
pair< std::vector<glm::vec3>, std::vector<glm::uvec3>> mDecimationPairPositionIndices2;
// Voxelisation need
std::vector<glm::uvec3> mVoxelTrianglesIndices;
std::vector<glm::vec3> mVoxelVertexPositions;
pair< std::vector<glm::vec3>, std::vector<glm::uvec3>> mVoxelPairPositionIndices;


//Tout needed For Mix Begin

pair< pair<int, std::vector<glm::vec3>>, pair<int, std::vector<glm::uvec3>>> DataStoreAfterVoxDec;


//Mix interpolation needed

std::vector<float> tranferGruppe;
pair< pair<int, std::vector<glm::vec3>>, pair<int, std::vector<glm::uvec3>>> CompleteDataInterpolationTranfer;
std::vector<glm::vec3> vPosition2;
//Rendering mode (0 : PBR, 1 : toon shading, 2 : x-toon shading)
static float Mode = 0.f;

void clear ();

#pragma region newlist
void Real_To_Two_List_Triangle_Position_Indices(std::vector<glm::vec3> VertexOriginPositions, std::vector<glm::uvec3>TrianglesIndicesOriginIndices, std::shared_ptr<Mesh> meshPtr) {
	/*TemporaireStockage.first.first.clear();
	TemporaireStockage.first.second.clear();
	TemporaireStockage.second.first.clear();
	TemporaireStockage.second.second.clear();*/
	std::cout << "[Real to Two List Split Loading]" << std::endl;
	std::cout << "<                              >" << std::endl;

	std::vector<glm::uvec3> listMacroTrianglesIndices;
	std::vector<glm::vec3> listMacroVertexPositions;
	std::vector<glm::uvec3> listMicroTrianglesIndices;
	std::vector<glm::vec3> listMicroVertexPositions;
	int nbrDecimation = 0, nbrDecimationT = 0, nbrIndices = 0, nbrIndicesT = 0;

	int indiceTriangle;
	glm::vec3 vertexPosition_1, vertexPosition_2, vertexPosition_3;
	glm::uvec3 indice_Triangle;
	for (size_t i = 0; i < meshPtr->component_Triangle_Real_MacroSurface.first.size(); i++)
	{
		indiceTriangle = meshPtr->component_Triangle_Real_MacroSurface.first[i];
		vertexPosition_1 = VertexOriginPositions[TrianglesIndicesOriginIndices[indiceTriangle][0]];
		vertexPosition_2 = VertexOriginPositions[TrianglesIndicesOriginIndices[indiceTriangle][1]];
		vertexPosition_3 = VertexOriginPositions[TrianglesIndicesOriginIndices[indiceTriangle][2]];

		if (meshPtr->component_Triangle_Real_MacroSurface.second[i] == 1) {

			pair<bool[3], int[3]> rsult(meshPtr->containsVector(listMacroVertexPositions, vertexPosition_1, vertexPosition_2, vertexPosition_3));
			
			nbrDecimationT++;
			if (rsult.first[0]) {
				indice_Triangle[0] = rsult.second[0];
			}
			else {
				listMacroVertexPositions.push_back(vertexPosition_1);
				indice_Triangle[0] = nbrDecimation;
				nbrDecimation++;
			}
			if (rsult.first[1]) {
				indice_Triangle[1] = rsult.second[1];
			}
			else {
				listMacroVertexPositions.push_back(vertexPosition_2);
				indice_Triangle[1] = nbrDecimation;
				nbrDecimation++;
			}if (rsult.first[2]) {
				indice_Triangle[2] = rsult.second[2];
			}
			else {
				listMacroVertexPositions.push_back(vertexPosition_3);
				indice_Triangle[2] = nbrDecimation;
				nbrDecimation++;
			}
			listMacroTrianglesIndices.push_back(indice_Triangle);
			if ((indice_Triangle[0] > 30000) || (indice_Triangle[1] > 30000) || (indice_Triangle[2] > 30000)) {
				std::cout << "Macro[" << indice_Triangle[0] << "," << indice_Triangle[1] << "," << indice_Triangle[2] << "] " << std::flush;
			}
			if ((indice_Triangle[0]<0) || (indice_Triangle[1]<0) || (indice_Triangle[2] <0)) {
				std::cout << "Macro_inf_[" << indice_Triangle[0] << "," << indice_Triangle[1] << "," << indice_Triangle[2] << "] " << std::flush;
			}
		}
		else {
			pair<bool[3], int[3]> rsult(meshPtr->containsVector(listMicroVertexPositions, vertexPosition_1, vertexPosition_2, vertexPosition_3));
			nbrIndicesT++;
		
			if (rsult.first[0]) {
				indice_Triangle[0] = rsult.second[0];
			}
			else {
				listMicroVertexPositions.push_back(vertexPosition_1);
				indice_Triangle[0] = nbrIndices;
				nbrIndices++;
			}
			if (rsult.first[1]) {
				indice_Triangle[1] = rsult.second[1];
			}
			else {
				listMicroVertexPositions.push_back(vertexPosition_2);
				indice_Triangle[1] = nbrIndices;
				nbrIndices++;
			}if (rsult.first[2]) {
				indice_Triangle[2] = rsult.second[2];
			}
			else {
				listMicroVertexPositions.push_back(vertexPosition_3);
				indice_Triangle[2] = nbrIndices;
				nbrIndices++;
			}
			listMicroTrianglesIndices.push_back(indice_Triangle);
			
			if ((indice_Triangle[0] > 30000) || (indice_Triangle[1] > 30000) || (indice_Triangle[2] > 30000)) {
				std::cout << "Micro[" << indice_Triangle[0] << "," << indice_Triangle[1] << "," << indice_Triangle[2] << "] " << std::flush;
			}
		}
	}
	
	std::cout << " [MacroPositions.size()/MacroTriangles.size()]=[" << listMacroVertexPositions.size() << " / " << listMacroTrianglesIndices.size() <<"]"<< std::endl;
	std::cout << "[MicroPositions.size() / MicroTriangles.size()] = [" << listMicroVertexPositions.size() << " / " << listMicroTrianglesIndices.size() << "]" << std::endl;
	
	/*TemporaireStockage.first.first = listMacroVertexPositions;
	TemporaireStockage.first.second = listMacroTrianglesIndices;
	TemporaireStockage.second.first = listMicroVertexPositions;
	TemporaireStockage.second.second = listMicroTrianglesIndices;*/
	meshPtr->real_Macro_List_Vertex=listMacroVertexPositions;
	meshPtr->real_Macro_List_Triangles = listMacroTrianglesIndices;
	meshPtr->real_Micro_List_Vertex = listMicroVertexPositions;
	meshPtr->real_Micro_List_Triangles = listMicroTrianglesIndices;
	std::cout << "<                               >" << std::endl;
	std::cout << "[Real to Two List Split Loaded]" << std::endl;
	std::cout << "" << std::endl;
	std::cout << "" << std::endl;
	
}

void Two_List_Triangles_Positions_Indices_To_New_List_Triangles_Positions_Indices(float taille_radius, std::shared_ptr<Mesh> meshPtr) {
	std::cout << "[Decimation And Voxelisation on the Two List Split Loading]" << std::endl;
	std::cout << "<                                                        >" << std::endl;


	//decimation
	std::cout << "" << std::endl;
	std::cout << "[---------------------------------------]" << std::endl;
	std::cout << "Begin  Decimation : " << "[" << meshPtr->real_Macro_List_Vertex.size() << "/" << meshPtr->real_Macro_List_Triangles.size() << "]" << std::endl;

	mDecimationPairPositionIndices.first.clear();
	mDecimationPairPositionIndices.second.clear();
	mDecimationPairPositionIndices = MeshDecimation::QEMDecimateForASpecifiedListOfEdges(meshPtr->real_Macro_List_Vertex, meshPtr->real_Macro_List_Triangles, taille_radius, meshPtr);
	

	std::cout << "End Decimation  : " << "[" << mDecimationPairPositionIndices.first.size() << "/" << mDecimationPairPositionIndices.second.size() << "]" << std::endl;
	std::cout << "[---------------------------------------]" << std::endl;

	std::cout << "" << std::endl;
	for (size_t i = 0; i < mDecimationPairPositionIndices.second.size(); i++)
	{
		if ((mDecimationPairPositionIndices.second[i][0] > 30000) || (mDecimationPairPositionIndices.second[i][0] > 30000) || (mDecimationPairPositionIndices.second[i][0] > 30000)) {
			std::cout << "mDEcim[" << mDecimationPairPositionIndices.second[i][0] << "," << mDecimationPairPositionIndices.second[i][1] << "," << mDecimationPairPositionIndices.second[i][2] << "] " << std::flush;
		}
	}
/*	for (size_t i = 0; i < meshPtr->vertex_transferGroupe.size(); i++)
	{
		std::cout << meshPtr->vertex_transferGroupe[i] << "/" << std::flush;
	}*/

	//std::cout <<" " << "/" << std::endl;
	/*for (size_t i = 0; i < meshPtr->vertex_transferGroupe.size(); i++)
	{
		std::cout << meshPtr->vertex_transferGroupe[i] << "/" << std::flush;
	}*/

	//voxelisation
	std::cout << "[---------------------------------------]" << std::endl;
	mVoxelPairPositionIndices.first.clear();
	mVoxelPairPositionIndices.second.clear();
	std::cout << "Begin Voxelization : " << "[" << meshPtr->real_Micro_List_Vertex.size() << "/" << meshPtr->real_Micro_List_Triangles.size() << "]" << std::endl;
	mVoxelPairPositionIndices = init_Voxelisation(meshPtr->real_Micro_List_Vertex, meshPtr->real_Micro_List_Triangles, taille_radius/2);
	std::cout << "End Voxelization : " << "[" << mVoxelPairPositionIndices.first.size() << "/" << mVoxelPairPositionIndices.second.size() << "]" << std::endl;
	for (size_t i = 0; i < mVoxelPairPositionIndices.second.size(); i++)
	{
		if ((mVoxelPairPositionIndices.second[i][0] > 30000000) || (mVoxelPairPositionIndices.second[i][0] > 30000000) || (mVoxelPairPositionIndices.second[i][0] > 30000000)) {
			std::cout << "mVox[" << mVoxelPairPositionIndices.second[i][0] << "," << mVoxelPairPositionIndices.second[i][1] << "," << mVoxelPairPositionIndices.second[i][2] << "] " << std::flush;
		}

	}
	std::cout << "[---------------------------------------]" << std::endl;
	std::cout << "<                                                        >" << std::endl;
	std::cout << "[Decimation And Voxelisation on the Two List Split Loaded]" << std::endl;
	std::cout << "" << std::endl;
	std::cout << "" << std::endl;
}

void New_Two_List_To_One( std::shared_ptr<Mesh> meshPtr) {
	std::cout << "[New List Concatenated Loading]" << std::endl;
	std::cout << "<                             >" << std::endl;
	std::cout << " OldList VerTex nbr et Triangles [Macro/Micro] : " << "{[" << meshPtr->real_Macro_List_Vertex.size() << "/" << meshPtr->real_Micro_List_Vertex.size() << "],[" << meshPtr->real_Macro_List_Triangles.size() << "/" << meshPtr->real_Micro_List_Triangles.size() << "]}" << std::endl;
	DataStoreAfterVoxDec.first.first = 0;
	DataStoreAfterVoxDec.second.first = 0;
	DataStoreAfterVoxDec.first.second.clear();
	DataStoreAfterVoxDec.second.second.clear();

	/*int sizePos = listPosTriIndSplit.first.first.size();
	int sizeInd = listPosTriIndSplit.first.second.size();
	DataStoreAfterVoxDec.first.first = sizePos;
	DataStoreAfterVoxDec.second.first = sizeInd;
	for (size_t i = 0; i < sizePos; i++)
	{
		DataStoreAfterVoxDec.first.second.push_back(listPosTriIndSplit.first.first[i]);
	}
	for (size_t i = 0; i < listPosTriIndSplit.second.first.size(); i++)
	{
		DataStoreAfterVoxDec.first.second.push_back(listPosTriIndSplit.second.first[i]);
	}
	for (size_t i = 0; i < sizeInd; i++)
	{
		DataStoreAfterVoxDec.second.second.push_back(listPosTriIndSplit.first.second[i]);
	}
	for (size_t i = 0; i < listPosTriIndSplit.second.second.size(); i++)
	{
		DataStoreAfterVoxDec.second.second.push_back(listPosTriIndSplit.second.second[i] += sizePos);
	}*/
	int sizePos = mDecimationPairPositionIndices.first.size();
	int sizeInd = mDecimationPairPositionIndices.second.size();
	DataStoreAfterVoxDec.first.first = sizePos;
	DataStoreAfterVoxDec.second.first = sizeInd;
	for (size_t i = 0; i < sizePos; i++)
	{

		DataStoreAfterVoxDec.first.second.push_back(mDecimationPairPositionIndices.first[i]);
	}
	for (size_t i = 0; i < mVoxelPairPositionIndices.first.size(); i++)
	{
		DataStoreAfterVoxDec.first.second.push_back(mVoxelPairPositionIndices.first[i]);

	}
	for (size_t i = 0; i < sizeInd; i++)
	{
		DataStoreAfterVoxDec.second.second.push_back(mDecimationPairPositionIndices.second[i]);
		//std::cout << "[" << DataStoreAfterVoxDec.second.second[i][0] << "," << DataStoreAfterVoxDec.second.second[i][1] << "," << DataStoreAfterVoxDec.second.second[i][2] << "] " << std::flush;
	}
	for (size_t i = 0; i < mVoxelPairPositionIndices.second.size(); i++)
	{
		DataStoreAfterVoxDec.second.second.push_back(mVoxelPairPositionIndices.second[i] += sizePos);
	}
	mDecimationPairPositionIndices.first.clear();
	mDecimationPairPositionIndices.second.clear();
	mVoxelPairPositionIndices.first.clear();
	mVoxelPairPositionIndices.second.clear();
	std::cout << " NewList : " << "{[" << DataStoreAfterVoxDec.first.first << "->" << DataStoreAfterVoxDec.first.second.size() << "],[" << DataStoreAfterVoxDec.second.first << "->" << DataStoreAfterVoxDec.second.second.size() << "]}" << std::endl;
	std::cout << "<                             >" << std::endl;
	std::cout << "[New List Concatenated Loaded]" << std::endl;
	for (size_t i = 0; i < DataStoreAfterVoxDec.second.second.size(); i++)
	{
		//std::cout <<"["<< DataStoreAfterVoxDec.second.second[i][0] << "," << DataStoreAfterVoxDec.second.second[i][1] << "," << DataStoreAfterVoxDec.second.second[i][2] << "] " << std::flush;
	}
}

void Make_New_List(float taille_radius, std::vector<int> TrianglesIndices, std::shared_ptr<Mesh> meshPtr) {

	std::cout << "[Potentiel Macro Surface Loading]" << std::endl;
	std::cout << "<                               >" << std::endl;
	meshPtr->init_component_Triangle_Potentiel_MacroSurfaces(TrianglesIndices, taille_radius);
	std::cout << "<                               >" << std::endl;
	std::cout << "[Potentiel Macro Surface Loaded]" << std::endl;
	std::cout << "" << std::endl;
	std::cout << "" << std::endl;


	std::cout << "[Real Macro Surface Loading] " << std::endl;
	std::cout << "<                          >" << std::endl;
	meshPtr->Potentiel_To_Real_Macro_Surfaces(taille_radius);
	std::cout << "<                          >" << std::endl;
	std::cout << "[Real Macro Surface Loaded]" << std::endl;
	std::cout << "" << std::endl;
	std::cout << "" << std::endl;

	std::cout << "[Real to Two List Split Loading]" << std::endl;
	std::cout << "<                              >" << std::endl;
	Real_To_Two_List_Triangle_Position_Indices(meshPtr->m_vertexPositions, meshPtr->m_triangleIndices, meshPtr);
	std::cout << "<                               >" << std::endl;
	std::cout << "[Real to Two List Split Loaded]" << std::endl;
	std::cout << "" << std::endl;
	std::cout << "" << std::endl;


	std::cout << "[Decimation And Voxelisation on the Two List Split Loading]" << std::endl;
	std::cout << "<                                                        >" << std::endl;
	Two_List_Triangles_Positions_Indices_To_New_List_Triangles_Positions_Indices( taille_radius, meshPtr);
	std::cout << "<                                                        >" << std::endl;
	std::cout << "[Decimation And Voxelisation on the Two List Split Loaded]" << std::endl;
	std::cout << "" << std::endl;
	std::cout << "" << std::endl;

	std::cout << "[New List Concatenated Loading]" << std::endl;
	std::cout << "<                             >" << std::endl;
	//std::cout << " OldList : " << "{[" << TemporaireStockage.first.first.size() << "/" << TemporaireStockage.second.first.size() << "],[" << TemporaireStockage.first.second.size() << "/" << TemporaireStockage.second.second.size() << "]}" << std::endl;
	New_Two_List_To_One( meshPtr);
	std::cout << " NewList : " << "{[" << DataStoreAfterVoxDec.first.first << "->" << DataStoreAfterVoxDec.first.second.size() << "],[" << DataStoreAfterVoxDec.second.first << "->" << DataStoreAfterVoxDec.second.second.size() << "]}" << std::endl;
	std::cout << "<                             >" << std::endl;
	std::cout << "[New List Concatenated Loaded]" << std::endl;


}
void CompleteDataInterpolation(std::vector<glm::vec3> real_Macro_List_Vertex, std::vector<glm::uvec3> real_Macro_List_Triangles, std::vector<glm::vec3> real_Micro_List_Vertex, std::vector<glm::uvec3> real_Micro_List_Triangles,float taille_radius) {
	pair< pair<std::vector<glm::vec3>, std::vector<glm::uvec3>>, pair< std::vector<glm::vec3>, std::vector<glm::uvec3>>> TemporaireStockage;
	TemporaireStockage.first.first = real_Macro_List_Vertex;
	TemporaireStockage.first.second = real_Macro_List_Triangles;
	TemporaireStockage.second.first = real_Micro_List_Vertex;
	TemporaireStockage.second.second =real_Micro_List_Triangles;

	CompleteDataInterpolationTranfer.first.first = TemporaireStockage.first.first.size();
	CompleteDataInterpolationTranfer.first.second = TemporaireStockage.first.first;
	CompleteDataInterpolationTranfer.second.first = TemporaireStockage.first.second.size();
	CompleteDataInterpolationTranfer.second.second = TemporaireStockage.first.second;
	int size_old_d = CompleteDataInterpolationTranfer.first.second.size();

	for (size_t i = 0; i < TemporaireStockage.second.first.size(); i++)
	{
		CompleteDataInterpolationTranfer.first.second.push_back(TemporaireStockage.second.first[i]);
	}
	std::cout << "OldDEcim + OldVox : " << CompleteDataInterpolationTranfer.first.second.size() << std::endl;
	int size_old_d_v = CompleteDataInterpolationTranfer.first.second.size();
	glm::uvec3 indice_triangle_tmp;
	for (size_t i = 0; i < TemporaireStockage.second.second.size(); i++)
	{
		indice_triangle_tmp = TemporaireStockage.second.second[i];
		indice_triangle_tmp += size_old_d;
		CompleteDataInterpolationTranfer.second.second.push_back(indice_triangle_tmp);
		if ((indice_triangle_tmp[0] > size_old_d_v) || (indice_triangle_tmp[1] > size_old_d_v) || (indice_triangle_tmp[2] > size_old_d_v)) {
			std::cout << "OldDEcim + OldVox traingles : Max  " << std::endl;
		}
	}
	std::cout << "OldDEcim + OldVox traingles : " << CompleteDataInterpolationTranfer.second.second.size() << std::endl;


	for (size_t i = DataStoreAfterVoxDec.first.first; i < DataStoreAfterVoxDec.first.second.size(); i++)
	{
		CompleteDataInterpolationTranfer.first.second.push_back(DataStoreAfterVoxDec.first.second[i]);
	}
	int size_old_d_v_n_v = CompleteDataInterpolationTranfer.first.second.size();
	std::cout << "OldDecim + OldVox + New Vox : " << CompleteDataInterpolationTranfer.first.second.size() << std::endl;
	int size_new = DataStoreAfterVoxDec.first.first;

	for (size_t i = DataStoreAfterVoxDec.second.first; i < DataStoreAfterVoxDec.second.second.size(); i++)
	{
		indice_triangle_tmp = (DataStoreAfterVoxDec.second.second[i] += (size_old_d_v - size_new));
		CompleteDataInterpolationTranfer.second.second.push_back(indice_triangle_tmp);
		if ((indice_triangle_tmp[0] > size_old_d_v_n_v) || (indice_triangle_tmp[2] > size_old_d_v_n_v) || (indice_triangle_tmp[2] > size_old_d_v_n_v)) {
			std::cout << "OldDEcim + OldVox traingles +Newtrianfles : Max  " << std::endl;
		}
		//std::cout << "[" << CompleteDataInterpolationTranfer.second.second[i][0] << "," << CompleteDataInterpolationTranfer.second.second[i][1] << "," << CompleteDataInterpolationTranfer.second.second[i][2] << "] " << std::flush;
	}
	std::cout << "OldDecim + OldVox + New Vox Triangle : " << CompleteDataInterpolationTranfer.second.second.size() << std::endl;
	std::cout << "\n" << CompleteDataInterpolationTranfer.first.second.size();
}
void CompleteTranferGruppe(std::vector<glm::vec3> real_Macro_List_Vertex, std::vector<glm::uvec3> real_Macro_List_Triangles, std::vector<glm::vec3> real_Micro_List_Vertex, std::vector<glm::uvec3> real_Micro_List_Triangles, float taille_radius) {
	pair< pair<std::vector<glm::vec3>, std::vector<glm::uvec3>>, pair< std::vector<glm::vec3>, std::vector<glm::uvec3>>> TemporaireStockage;
	TemporaireStockage.first.first = real_Macro_List_Vertex;
	TemporaireStockage.first.second = real_Macro_List_Triangles;
	TemporaireStockage.second.first = real_Micro_List_Vertex;
	TemporaireStockage.second.second = real_Micro_List_Triangles;
	std::cout << "vexter Gruppe : " << meshPtr->vertex_transferGroupe.size();
	int cmp = 0, cmp_ = 0;
	int final_, tmp_ = 0;
	tranferGruppe.resize(meshPtr->vertex_transferGroupe.size(), -1);
	std::cout << "[" << std::endl;
	std::vector<int> tmp_tmp;

	for (size_t i = 0; i < meshPtr->vertex_transferGroupe.size(); i++)
	{
		//std::cout << "[" << i << "->" << meshPtr->vertex_transferGroupe[i] << "] " << std::flush;
		tmp_ = meshPtr->vertex_transferGroupe[i];
		if (tmp_ != -1) {
			cmp_++;
			while (tmp_ != -1) {
				final_ = tmp_;
				tmp_ = meshPtr->vertex_transferGroupe[tmp_];
			}
			tmp_tmp.push_back(cmp_);
			tranferGruppe[i] = final_;

		}
		else {
			cmp++;
			tmp_tmp.push_back(cmp_);
		}


	}
	for (size_t i = 0; i < meshPtr->vertex_transferGroupe.size(); i++)
	{

		if (meshPtr->vertex_transferGroupe[i] != -1) {
			tranferGruppe[i] = tranferGruppe[i] - tmp_tmp[tranferGruppe[i]];
		}

	}
	std::cout << "]" << std::endl;
	std::cout << "nbr de -1 : " << cmp << "  " << std::endl;



	for (size_t i = meshPtr->vertex_transferGroupe.size(); i < CompleteDataInterpolationTranfer.first.second.size(); i++)
	{
		if (i < TemporaireStockage.first.first.size() + TemporaireStockage.second.first.size()) {
			tranferGruppe.push_back(-1.5f);
		}
		else {
			tranferGruppe.push_back(-2.f);
		}


	}


	meshPtr->m_vertexPositions_NEW.resize(CompleteDataInterpolationTranfer.first.second.size(), glm::vec3(0.f, 0.f, 0.f));
	for (size_t i = 0; i < TemporaireStockage.first.first.size(); i++)
	{
		if (tranferGruppe[i] > -0.5f) {
			meshPtr->m_vertexPositions_NEW[i] = DataStoreAfterVoxDec.first.second[tranferGruppe[i]];
		}
		else {
			//meshPtr->m_vertexPositions_NEW[i] = glm::vec3(0, 0, 0);
		}
	}
	vPosition2 = meshPtr->m_vertexPositions_NEW;
	TemporaireStockage.first.first.clear();
	TemporaireStockage.first.second.clear();
	TemporaireStockage.second.first.clear();
	TemporaireStockage.second.second.clear();
}
#pragma endregion

void LoadAllDataFromFiles(int index) {
	meshPtr->clear();
	meshPtr = std::make_shared<Mesh>();
	try {
		MeshLoader::loadOFF("Resources/mesh_collection/"+ListNameMesh[index]+".off", meshPtr);
	}
	catch (std::exception & e) {
		//exitOnCriticalError(std::string("[Error loading mesh]") + e.what());
	}
	try {
		MeshLoader::loadEdges(ListNameMesh[index]+ ListNameObjectToLoad[0], meshPtr);
	}
	catch (std::exception & e) {
		//exitOnCriticalError(std::string("[Error loading edges]") + e.what());
	}
	try {
		MeshLoader::loadPotentielAndRealMacrosurface(ListNameMesh[index] + ListNameObjectToLoad[1], meshPtr);
	}
	catch (std::exception & e) {
		//exitOnCriticalError(std::string("[Error loading edges]") + e.what());
	}
	try {
		MeshLoader::loadDecimationAndVoxel(ListNameMesh[index] + ListNameObjectToLoad[2], meshPtr);
	}
	catch (std::exception & e) {
	//	exitOnCriticalError(std::string("[Error loading edges]") + e.what());
	}
	try {
		MeshLoader::loadScaleDataInterpolation(ListNameMesh[index] + ListNameObjectToLoad[3], meshPtr);
	}
	catch (std::exception & e) {
		//exitOnCriticalError(std::string("[Error loading edges]") + e.what());
	}
	meshPtr->m_vertexOriginPositions = meshPtr->m_vertexPositions;
	meshPtr->m_triangleOriginIndices = meshPtr->m_triangleIndices;
}
/*void WriteAllDataTOFiles(int index) {
	try {
		MeshLoader::writeEdges(ListNameMesh[index] + ListNameObjectToLoad[0],meshPtr);
	}
	catch (std::exception & e) {
		//exitOnCriticalError(std::string("[Error loading edges]") + e.what());
	}
	try {
		MeshLoader::writePotentielAndRealMacrosurface(ListNameMesh[index] + ListNameObjectToLoad[1], meshPtr);
	}
	catch (std::exception & e) {
		//exitOnCriticalError(std::string("[Error loading edges]") + e.what());
	}
	try {
		MeshLoader::writeDecimationAndVoxel(ListNameMesh[index] + ListNameObjectToLoad[2], meshPtr);
	}
	catch (std::exception & e) {
		//	exitOnCriticalError(std::string("[Error loading edges]") + e.what());
	}
	try {
		MeshLoader::writeScaleDataInterpolation(ListNameMesh[index] + ListNameObjectToLoad[3],  CompleteDataInterpolationTranfer, taille_radius, tranferGruppe, vPosition2);
	}
	catch (std::exception & e) {
		//exitOnCriticalError(std::string("[Error loading edges]") + e.what());
	}
	
}*/

void printHelp () {
	std::cout << "> Help:" << std::endl
			  << "    Mouse commands:" << std::endl
			  << "    * Left button: rotate camera" << std::endl
			  << "    * Middle button: zoom" << std::endl
			  << "    * Right button: pan camera" << std::endl
			  << "    Keyboard commands:" << std::endl
   			  << "    * H: print this help" << std::endl
   			  << "    * F1: toggle wireframe rendering" << std::endl
   			  << "    * ESC: quit the program" << std::endl;
}

// Executed each time the window is resized. Adjust the aspect ratio and the rendering viewport to the current window.
void windowSizeCallback (GLFWwindow * windowPtr, int width, int height) {
	cameraPtr->setAspectRatio (static_cast<float>(width) / static_cast<float>(height));
	glViewport (0, 0, (GLint)width, (GLint)height); // Dimension of the rendering region in the window
}

/// Executed each time a key is entered.
void keyCallback (GLFWwindow * windowPtr, int key, int scancode, int action, int mods) {
	if (action == GLFW_PRESS && key == GLFW_KEY_H) {
		printHelp ();
	}
	else if (action == GLFW_PRESS && key == GLFW_KEY_F1) {
		GLint mode[2];
		glGetIntegerv (GL_POLYGON_MODE, mode);
		glPolygonMode (GL_FRONT_AND_BACK, mode[1] == GL_FILL ? GL_LINE : GL_FILL);
	}
	else if (action == GLFW_PRESS && key == GLFW_KEY_ESCAPE) {
		glfwSetWindowShouldClose (windowPtr, true); // Closes the application if the escape key is pressed
	}
	else if (action == GLFW_PRESS && key == GLFW_KEY_N) {
		//les 3 modes possibles
		Mode = Mode + 1.0f;
		if (Mode ==  5.0f) {
			Mode = 0.0f;
		}
		if (Mode == 3.0f) {
			meshPtr->vertex_CurrentGroupe = meshPtr->trianglesCurrentGroupToVertexCurrentGroup(meshPtr->component_Triangle_Potentiel_MacroSurface);
			meshPtr->init(taille_radius,false);
		}
		if (Mode == 4.0f) {
			meshPtr->vertex_CurrentGroupe = meshPtr->trianglesCurrentGroupToVertexCurrentGroup(meshPtr->component_Triangle_Real_MacroSurface);
			meshPtr->init(taille_radius,false);
		}
	
	}
	else if (action == GLFW_PRESS && key == GLFW_KEY_B) {
		//les 3 modes possibles
		Mode = Mode - 1.0f;
		if (Mode == -1.0f) {
			Mode = 0.0f;
		}
		if (Mode == 3.0f) {
			meshPtr->vertex_CurrentGroupe = meshPtr->trianglesCurrentGroupToVertexCurrentGroup(meshPtr->component_Triangle_Potentiel_MacroSurface);
			meshPtr->init(taille_radius,false);
		}
		if (Mode == 4.0f) {
			meshPtr->vertex_CurrentGroupe = meshPtr->trianglesCurrentGroupToVertexCurrentGroup(meshPtr->component_Triangle_Real_MacroSurface);
			meshPtr->init(taille_radius,false);
		}
	}
	else if (action == GLFW_PRESS && key == GLFW_KEY_1) {
		
		GLint mode[2];
		glGetIntegerv(GL_POLYGON_MODE, mode);
		glPolygonMode(GL_FRONT_AND_BACK, mode[1] == GL_FILL ? GL_LINE : GL_FILL);
		

	}	
	else if (action == GLFW_PRESS && key == GLFW_KEY_R) {
		std::cout << "R" << std::endl;
		Mode = 0.f;
		meshPtr->m_vertexPositions = meshPtr->m_vertexOriginPositions;
		meshPtr->m_triangleIndices = meshPtr->m_triangleOriginIndices;
		meshPtr->vertex_CurrentGroupe = meshPtr->trianglesCurrentGroupToVertexCurrentGroup(meshPtr->component_Triangle_Potentiel_MacroSurface);
		meshPtr->recomputePerVertexNormals();
		meshPtr->init(taille_radius,false);
		
	}
	else if (action == GLFW_PRESS && key == GLFW_KEY_D) {
		std::cout << "D" << std::endl;
		meshPtr->m_vertexPositions = m_vertexDecimationPositions;
		meshPtr->m_triangleIndices = m_triangleDecimationIndices;
		meshPtr->recomputePerVertexNormals();
		meshPtr->init(taille_radius,false);
		
	}
	else if (action == GLFW_PRESS && key == GLFW_KEY_V) {
		std::cout << "V"<<std::endl;
		meshPtr->m_vertexPositions = mVoxelVertexPositions;
		meshPtr->m_triangleIndices = mVoxelTrianglesIndices;
		meshPtr->recomputePerVertexNormals();
		meshPtr->computePlanarParameterization();
		meshPtr->init(taille_radius,false);
	}

	else if (action == GLFW_PRESS && key == GLFW_KEY_P) {
		std::cout << "P" << std::endl;
		meshPtr->m_vertexPositions = DataStoreAfterVoxDec.first.second;
		meshPtr->m_triangleIndices = DataStoreAfterVoxDec.second.second;
		std::cout << "Normals" << std::endl;
		meshPtr->recomputePerVertexNormals();
		std::cout << "Planar" << std::endl;
		meshPtr->computePlanarParameterization();
		std::cout << "Init" << std::endl;
		meshPtr->init(taille_radius,false);
	}
	else if (action == GLFW_PRESS && key == GLFW_KEY_O) {
		std::cout << "O" << std::endl;
		Mode = -3.f;
		meshPtr->m_vertexPositions = CompleteDataInterpolationTranfer.first.second;
		meshPtr->m_triangleIndices = CompleteDataInterpolationTranfer.second.second;
		meshPtr->vertex_CurrentGroupe = tranferGruppe;
		meshPtr->m_vertexPositions_NEW = vPosition2;
		meshPtr->recomputePerVertexNormals();
		meshPtr->computePlanarParameterization();
		meshPtr->init(taille_radius,true);
	}

}

/// Called each time the mouse cursor moves
void cursorPosCallback(GLFWwindow* window, double xpos, double ypos) {
	 
	if (!(((xpos > test[0]) && (xpos < test[0] + test2[0])) || ((ypos > test[1]) && (ypos < test[1] + test2[1])))){
	int width, height;
	glfwGetWindowSize(windowPtr, &width, &height);
	float normalizer = static_cast<float> ((width + height) / 2);
	float dx = static_cast<float> ((baseX - xpos) / normalizer);
	float dy = static_cast<float> ((ypos - baseY) / normalizer);
	if (isRotating) {
		glm::vec3 dRot(-dy * M_PI, dx * M_PI, 0.0);
		cameraPtr->setRotation(baseRot + dRot);
	}
	else if (isPanning) {
		cameraPtr->setTranslation(baseTrans + meshScale * glm::vec3(dx, dy, 0.0));
	}
	else if (isZooming) {
		cameraPtr->setTranslation(baseTrans + meshScale * glm::vec3(0.0, 0.0, dy));
	}
	}
}

/// Called each time a mouse button is pressed
void mouseButtonCallback (GLFWwindow * window, int button, int action, int mods) {
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
    	if (!isRotating) {
    		isRotating = true;
    		glfwGetCursorPos (window, &baseX, &baseY);
    		baseRot = cameraPtr->getRotation ();
        }
    } else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
    	isRotating = false;
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
    	if (!isPanning) {
    		isPanning = true;
    		glfwGetCursorPos (window, &baseX, &baseY);
    		baseTrans = cameraPtr->getTranslation ();
        }
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE) {
    	isPanning = false;
    } else if (button == GLFW_MOUSE_BUTTON_MIDDLE && action == GLFW_PRESS) {
    	if (!isZooming) {
    		isZooming = true;
    		glfwGetCursorPos (window, &baseX, &baseY);
    		baseTrans = cameraPtr->getTranslation ();
        }
    } else if (button == GLFW_MOUSE_BUTTON_MIDDLE && action == GLFW_RELEASE) {
    	isZooming = false;
    }
}

void initGLFW () {
	// Initialize GLFW, the library responsible for window management
	if (!glfwInit ()) {
		std::cerr << "ERROR: Failed to init GLFW" << std::endl;
		std::exit (EXIT_FAILURE);
	}

	// Before creating the window, set some option flags
	glfwWindowHint (GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint (GLFW_CONTEXT_VERSION_MINOR, 5);
	glfwWindowHint (GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint (GLFW_RESIZABLE, GL_TRUE);

	// Create the window
	windowPtr = glfwCreateWindow (1024, 768, "Computer Graphics - Practical Assignment", nullptr, nullptr);
	if (!windowPtr) {
		std::cerr << "ERROR: Failed to open window" << std::endl;
		glfwTerminate ();
		std::exit (EXIT_FAILURE);
	}

	// Load the OpenGL context in the GLFW window using GLAD OpenGL wrangler
	glfwMakeContextCurrent (windowPtr);

	/// Connect the callbacks for interactive control
	glfwSetWindowSizeCallback (windowPtr, windowSizeCallback);
	glfwSetKeyCallback (windowPtr, keyCallback);
	glfwSetCursorPosCallback(windowPtr, cursorPosCallback);
	glfwSetMouseButtonCallback (windowPtr, mouseButtonCallback);
	// Initialize OpenGL loader
#if defined(IMGUI_IMPL_OPENGL_LOADER_GL3W)
	bool err = gl3wInit() != 0;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLEW)
	bool err = glewInit() != GLEW_OK;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLAD)
	bool err = gladLoadGL() == 0;
#else
	bool err = false; // If you use IMGUI_IMPL_OPENGL_LOADER_CUSTOM, your loader is likely to requires some form of initialization.
#endif
	

	// Setup Dear ImGui context
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;
	//io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
	//io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

	// Setup Dear ImGui style
	ImGui::StyleColorsDark();
	// Setup Platform/Renderer bindings
	ImGui_ImplGlfw_InitForOpenGL(windowPtr, true);
	ImGui_ImplOpenGL3_Init(glwfversion);
}

void exitOnCriticalError (const std::string & message) {
	std::cerr << "> [Critical error]" << message << std::endl;
	std::cerr << "> [Clearing resources]" << std::endl;
	clear ();
	std::cerr << "> [Exit]" << std::endl;
	std::exit (EXIT_FAILURE);
}

void initOpenGL () {
	// Load extensions for modern OpenGL
	if (!gladLoadGLLoader ((GLADloadproc)glfwGetProcAddress))
		exitOnCriticalError ("[Failed to initialize OpenGL context]");

	glEnable (GL_DEBUG_OUTPUT); // Modern error callback functionnality
	glEnable (GL_DEBUG_OUTPUT_SYNCHRONOUS); // For recovering the line where the error occurs, set a debugger breakpoint in DebugMessageCallback
    glDebugMessageCallback (debugMessageCallback, 0); // Specifies the function to call when an error message is generated.
	glCullFace (GL_BACK);     // Specifies the faces to cull (here the ones pointing away from the camera)
	glEnable (GL_CULL_FACE); // Enables face culling (based on the orientation defined by the CW/CCW enumeration).
	glDepthFunc (GL_LESS); // Specify the depth test for the z-buffer
	glEnable (GL_DEPTH_TEST); // Enable the z-buffer test in the rasterization
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// Loads and compile the programmable shader pipeline
	try {
		shaderProgramPtr = ShaderProgram::genBasicShaderProgram (SHADER_PATH + "VertexShader.glsl",
													         	 SHADER_PATH + "FragmentShader.glsl");
	} catch (std::exception & e) {
		exitOnCriticalError (std::string ("[Error loading shader program]") + e.what ());
	}
}


LightSource lightSource1;
LightSource lightSource2;

void initScene(const std::string & meshFilename) {
	
#pragma region Camera
	// Camera
	int width, height;
	glfwGetWindowSize (windowPtr, &width, &height);
	cameraPtr = std::make_shared<Camera> ();
	cameraPtr->setAspectRatio (static_cast<float>(width) / static_cast<float>(height));
#pragma endregion

#pragma region MeshLoading

	meshPtr = std::make_shared<Mesh> ();
	try {
		MeshLoader::loadOFF (meshFilename, meshPtr);
	} catch (std::exception & e) {
		exitOnCriticalError (std::string ("[Error loading mesh]") + e.what ());
	}
	//meshPtr->computeBoundingSphere();
#pragma endregion
	pair<glm::vec3, float> tmp_cen;
	tmp_cen=meshPtr->analyseBasicGeoStat();
	center = tmp_cen.first;
	meshScale = tmp_cen.second;

#pragma region HalfEdges
	//charger les données de half edges à partir d'un fichier
	try {
	//	MeshLoader::loadEdges(MESH_NAME+"_edges.off", meshPtr);
	}
	catch (std::exception & e) {
		exitOnCriticalError(std::string("[Error loading edges]") + e.what());
	}
	//permet de construire une structure similaire qui permet d'obtenir le triangle voisins de chaque edge et l'ecrire dans un fichier
	
	meshPtr->computeTriangleAdjEdges();
	try {
		MeshLoader::writeEdges(MESH_NAME+"_edges.off", meshPtr);
	}
	catch (std::exception & e) {
		exitOnCriticalError(std::string("[Error wrtting edges 2]") + e.what());
	}
	
#pragma endregion
	
//initialisation de l'affichage des groupes des vertex pour éviter les problèmes d'affichages si il n'est pas utilisé
std::vector<float> tmp;
tmp.resize(meshPtr->m_vertexPositions.size(), -1);
meshPtr->vertex_CurrentGroupe = tmp;	
	
#pragma region PotentielMacroRegions
	std::vector<int> Triangles;
	for (int i = 0;i < meshPtr->m_triangleIndices.size();i++) {
		Triangles.push_back(i);
	}
	meshPtr->init_component_Triangle_Potentiel_MacroSurfaces(Triangles,taille_radius);

#pragma endregion
	
#pragma region RealRegions
	//************Initialisation de la liste des Real Triangles Macro Surfaces**************//
	meshPtr->Potentiel_To_Real_Macro_Surfaces(taille_radius);
	//************Initialisation de la liste des Real Triangles Macro Surfaces**************//


#pragma endregion
	try {
		//MeshLoader::writePotentielAndRealMacrosurface(MESH_NAME + "_PotentielAndRealSurfaces.off", meshPtr);
	}
	catch (std::exception & e) {
		exitOnCriticalError(std::string("[Error writting Potentiel And Real Surfaces]") + e.what());
	}

	try {
		//MeshLoader::loadPotentielAndRealMacrosurface(MESH_NAME + "_PotentielAndRealSurfaces.off", meshPtr);
	}
	catch (std::exception & e) {
		exitOnCriticalError(std::string("[Error writting Potentiel And Real Surfaces]") + e.what());
	}
#pragma region BufferInitialisation
	//************Initiation du Mesh pour préparer les buffers par défault**************//
	meshPtr->init(taille_radius,false);
	//************Initiation du Mesh pour préparer les buffers par défault**************/
#pragma endregion 
	//Pour permettre de reset le mesh : stockage des données d'origine dans un attribut de Mesh
	meshPtr->m_vertexOriginPositions = meshPtr->m_vertexPositions;
	meshPtr->m_triangleOriginIndices = meshPtr->m_triangleIndices;
	
	/*std::cout << "Nombre de Triangles : " << meshPtr->m_triangleIndices.size() << std::endl;
	std::cout << "Nombre de Vertex : " << meshPtr->m_vertexPositions.size() << std::endl;
	std::cout << "Nombre de Vertex issues de la couleur  : " << meshPtr->vertex_CurrentGroupe.size() << std::endl;*/
#pragma region Decimation
	
	
	/*cout <<" 1.2f : "<< meshPtr->NumbreOfVertexToDecimeToRespectSizeScale(1.2f, meshPtr->m_triangleIndices, meshPtr->m_vertexPositions) << std::endl;
	cout << " 0.6f : " << meshPtr->NumbreOfVertexToDecimeToRespectSizeScale(0.6f, meshPtr->m_triangleIndices, meshPtr->m_vertexPositions) << std::endl;
	cout << " 0.3f : " << meshPtr->NumbreOfVertexToDecimeToRespectSizeScale(0.3f, meshPtr->m_triangleIndices, meshPtr->m_vertexPositions) << std::endl;
	cout << " 0.15f : " << meshPtr->NumbreOfVertexToDecimeToRespectSizeScale(0.15f, meshPtr->m_triangleIndices, meshPtr->m_vertexPositions) << std::endl;
	cout << " 0.07125f : " << meshPtr->NumbreOfVertexToDecimeToRespectSizeScale(0.07125f, meshPtr->m_triangleIndices, meshPtr->m_vertexPositions) << std::endl;*/
	/*m_pairEdgesToDecime = meshPtr->WichEdgesDecimeToRespectSizeScale(taille_radius, meshPtr->m_triangleIndices, meshPtr->m_vertexPositions);
	cout << "Le nombre de Edges a Decimer pour un  float_radius : " << taille_radius << " : "<< m_pairEdgesToDecime.first<< std::endl;
	//vérifier les intersections 
	/*bool i1 = meshPtr->intersectSphere(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(1.5f, 0.0f, 0.0f), glm::vec3(0.0f, 0.0f, 0.0f), 1.0f);
	std::cout << "p1= [0.0f, 0.0f, 0.0f]  |"  <<" p2= [1.5f, 0.0f, 0.0f]   |"  << " R = 1   |" << " Intersection : " << i1 <<  std::endl;
	bool i2= meshPtr->intersectSphere(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.5f, 0.0f, 0.0f), glm::vec3(0.0f, 0.0f, 0.0f), 1.0f);
	std::cout << "p1= [0.0f, 0.0f, 0.0f]  |" << " p2= [0.5f, 0.0f, 0.0f]   |" << " R = 1   |" << " Intersection : " << i2 << std::endl;
	bool i3 = meshPtr->intersectSphere(glm::vec3(1.2f, 0.0f, 0.0f), glm::vec3(1.5f, 0.0f, 0.0f), glm::vec3(0.0f, 0.0f, 0.0f), 1.0f);
	std::cout << "p1= [1.2f, 0.0f, 0.0f]  |" << " p2= [1.5f, 0.0f, 0.0f]   |" << " R = 1   |" << " Intersection : " << i3 << std::endl;
	bool i4 = meshPtr->intersectSphere(glm::vec3(0.5f, 1.5f, 0.0f), glm::vec3(-0.5f, -1.5f, 0.0f), glm::vec3(0.0f, 0.0f, 0.0f), 1.0f);
	std::cout << "p1= [0.5f, 1.5f, 0.0f]  |" << " p2= [-0.5f, -1.5f, 0.0f] |" << " R = 1   |" << " Intersection : " << i4 << std::endl;
	bool i5 = meshPtr->intersectSphere(glm::vec3(1.0f, -1.0f, 0.0f), glm::vec3(1.0f, 1.0f, 0.0f), glm::vec3(0.0f, 0.0f, 0.0f), 1.0f);
	std::cout << "p1= [1.0f, -1.0f, 0.0f] |" << " p2= [1.0f, 1.0f, 0.0f]   |" << " R = 1   |" << " Intersection : " << i5 << std::endl;
	bool i5_2 = meshPtr->intersectSphere(glm::vec3(1.0f, -1.0f, 0.0f), glm::vec3(1.0f, 1.0f, 0.0f), glm::vec3(0.0f, 0.0f, 0.0f), 0.9f);
	std::cout << "p1= [1.0f, -1.0f, 0.0f] |" << " p2= [1.0f, 1.0f, 0.0f]   |" << " R = 0.9 |" << " Intersection : " << i5_2 << std::endl;
	bool i6 = meshPtr->intersectSphere(glm::vec3(1.0f, 1.0f, 1.0f), glm::vec3(0.5f, 0.0f, 0.0f), glm::vec3(-8.0f, 0.0f, 0.0f), 1.0f);
	std::cout << "p1= [1.0f, 1.0f, 1.0f]  |" << " p2= [0.5f, 0.0f, 0.0f]   |" << " R = 1   |" << " Intersection : " << i6 <<  " position du centre r=[-0.8f, 0.0f, 0.0f]"<<std::endl; 
	try {
		MeshLoader::writeEdges("eddes2.off", meshPtr->triangle_Adj_Edges);
	}
	catch (std::exception & e) {
		exitOnCriticalError(std::string("[Error wrtting edges]") + e.what());
	}*/

	meshPtr->m_DecimationPairPositionIndices =MeshDecimation::QEMDecimateForASpecifiedListOfEdges(meshPtr->m_vertexPositions, meshPtr->m_triangleIndices, taille_radius,meshPtr);
	meshPtr->m_vertexDecimationPositions = meshPtr->m_DecimationPairPositionIndices.first;
	meshPtr->m_triangleDecimationIndices = meshPtr->m_DecimationPairPositionIndices.second;
	
	//m_pairEdgesToDecime = meshPtr->WichEdgesDecimeToRespectSizeScale(taille_radius, mDecimationPairPositionIndices.second, mDecimationPairPositionIndices.first);
	//cout << "Le nombre de Edges a Decimer pour un  float_radius : " << taille_radius << " : " << m_pairEdgesToDecime.first << std::endl;
	
#pragma endregion 

#pragma region Voxelisation
	//*******************************Voxelisation******************************************//
	meshPtr->m_VoxelPairPositionIndices = init_Voxelisation(meshPtr->m_vertexPositions, meshPtr->m_triangleIndices, taille_radius);
	meshPtr->m_vertexVoxelPositions = meshPtr->m_VoxelPairPositionIndices.first;
	meshPtr->m_triangleVoxelIndices = meshPtr->m_VoxelPairPositionIndices.second;
	
#pragma endregion 

#pragma region CreationOfFirstStepScaleFixed

	meshPtr->Make_New_List_V2(taille_radius, Triangles, meshPtr);
	/*meshPtr->init_component_Triangle_Potentiel_MacroSurfaces(Triangles, taille_radius);
	meshPtr->Potentiel_To_Real_Macro_Surfaces(taille_radius);

    Real_To_Two_List_Triangle_Position_Indices(meshPtr->m_vertexPositions, meshPtr->m_triangleIndices, meshPtr);
	Two_List_Triangles_Positions_Indices_To_New_List_Triangles_Positions_Indices(taille_radius, meshPtr);
	New_Two_List_To_One( meshPtr); //DataStoreAfter
		
#pragma endregion

#pragma region Interpolation

   CompleteDataInterpolation(meshPtr->real_Macro_List_Vertex, meshPtr->real_Macro_List_Triangles, meshPtr->real_Micro_List_Vertex, meshPtr->real_Micro_List_Triangles, taille_radius);
	
	CompleteTranferGruppe(meshPtr->real_Macro_List_Vertex, meshPtr->real_Macro_List_Triangles, meshPtr->real_Micro_List_Vertex, meshPtr->real_Micro_List_Triangles, taille_radius);
	*/
#pragma endregion

#pragma region WrittingInDocuments
/*	try {
		MeshLoader::writeScaleDataInterpolation(MESH_NAME+"_DataInterpolation.off", CompleteDataInterpolationTranfer,  taille_radius, tranferGruppe,  vPosition2);
	}
	catch (std::exception & e) {
		exitOnCriticalError(std::string("[Error writting DataInterpolation]") + e.what());
	}
	try {
		MeshLoader::loadScaleDataInterpolation(MESH_NAME+"_DataInterpolation.off", meshPtr);
	}
	catch (std::exception & e) {
		exitOnCriticalError(std::string("[Error loading DataInterpolation]") + e.what());
	}
//	MeshLoader::showScaleDataInterpolation(meshPtr);

	

	try {
		MeshLoader::writeDecimationAndVoxel(MESH_NAME+"_DecimationAndVoxel.off", meshPtr);
	}
	catch (std::exception & e) {
		exitOnCriticalError(std::string("[Error writting Decimation and Voxel]") + e.what());
	}

	try {
		MeshLoader::loadDecimationAndVoxel(MESH_NAME+"_DecimationAndVoxel.off", meshPtr);
	}
	catch (std::exception & e) {
		exitOnCriticalError(std::string("[Error writting Decimation and Voxel]") + e.what());
	}*/


	
#pragma endregion
#pragma region Lights 

	
	
	// creation des lights
	lightSource1 = LightSource(center + glm::vec3(0.0, 0.0, 3.0 ), glm::vec3 (0.0, 0.5, 0.5), 3000.f,  glm::vec3(0.1, 0.1, -1.0) ,2.14);
	lightSource1.setAc(0.1f);
	lightSource1.setAq(0.1f);
	lightSource1.setAl(0.2f);
	//pousser les caractéristiques de la source
		shaderProgramPtr->set ("lightSource1.color", lightSource1.getColor());
		shaderProgramPtr->set ("lightSource1.intensity", lightSource1.getIntensity());
		shaderProgramPtr->set ("lightSource1.ac", lightSource1.getAc());
		shaderProgramPtr->set ("lightSource1.al", lightSource1.getAl());
		shaderProgramPtr->set ("lightSource1.aq", lightSource1.getAq());
		shaderProgramPtr->set("lightSource1.coneAngle", lightSource1.getConeAngle());
		
	
#pragma endregion

#pragma region Material
	// Material
	Material material = Material(glm::vec3 (0.4, 0.6, 0.2), 0.01, glm::vec3 (0.90, 0.91, 0.92));

	string dirName = "Resources\\Materials\\Metal\\";
	

	GLuint roughnessTex = material.loadTextureFromFileToGPU(dirName + "Roughness.png");

	GLuint metallicTex = material.loadTextureFromFileToGPU(dirName + "Metallic.png");

	GLuint albedoTex = material.loadTextureFromFileToGPU(dirName + "Base_Color.png");

	GLuint toonTex = material.loadTextureFromFileToGPU(dirName + "XTOOON.png");

	shaderProgramPtr->set ("material.albedoTex", 0u);
	shaderProgramPtr->set ("material.roughnessTex", 1u);
	shaderProgramPtr->set ("material.metallicTex", 2u);
	shaderProgramPtr->set ("material.toonTex", 3u);

	glActiveTexture (GL_TEXTURE0);
	glBindTexture (GL_TEXTURE_2D, albedoTex);

	glActiveTexture (GL_TEXTURE1);
	glBindTexture (GL_TEXTURE_2D, roughnessTex);

	glActiveTexture (GL_TEXTURE2);
	glBindTexture (GL_TEXTURE_2D, metallicTex);


	glActiveTexture (GL_TEXTURE3);
	glBindTexture (GL_TEXTURE_2D, toonTex);
#pragma endregion

#pragma region Xtoon 
	//zMin and zMax for the computation of the detail value
	shaderProgramPtr->set ("z_Min", meshScale);
	shaderProgramPtr->set ("z_Max", meshScale*3);
	shaderProgramPtr->set("Interpolate", interpolate);
#pragma endregion

#pragma region AdjustCamera
	// Adjust the camera to the actual mesh
	
	
	cameraPtr->setTranslation (center + glm::vec3 (0.0, 0.0, 3.0 * meshScale));
	cameraPtr->setNear (meshScale / 100.f);
	cameraPtr->setFar (6.f * meshScale);
	
#pragma endregion 

}

void init (const std::string & meshFilename) {
	initGLFW (); // Windowing system
	initOpenGL (); // OpenGL Context and shader pipeline
	initScene (meshFilename); // Actual scene to render
}

void clear () {
	cameraPtr.reset ();
	meshPtr.reset ();
	shaderProgramPtr.reset ();
	glfwDestroyWindow (windowPtr);
	glfwTerminate ();
}

// The main rendering call
void render () {
	ImGui::Render();
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Erase the color and z buffers.

	// background suivant le mode
	
	glClearColor (color[0], color[1], color[2], color[3] );
	
		

	shaderProgramPtr->use (); // Activate the program to be used for upcoming primitive
	
	
	glm::mat4 projectionMatrix = cameraPtr->computeProjectionMatrix ();
	shaderProgramPtr->set ("projectionMat", projectionMatrix); // Compute the projection matrix of the camera and pass it to the GPU program
	glm::mat4 modelMatrix = meshPtr->computeTransformMatrix ();
	glm::mat4 viewMatrix = cameraPtr->computeViewMatrix ();
	glm::mat4 modelViewMatrix = viewMatrix * modelMatrix;
	glm::mat4 normalMatrix = glm::transpose (glm::inverse (modelViewMatrix));
	shaderProgramPtr->set ("modelViewMat", modelViewMatrix);
	shaderProgramPtr->set ("normalMat", normalMatrix);
	meshPtr->render ();
	shaderProgramPtr->stop ();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

// Update any accessible variable based on the current time
void update (float currentTime) {
	if (oldMEshLoaded != newMEshLoaded) {
		/*oldMEshLoaded = newMEshLoaded;
		LoadAllDataFromFiles(oldMEshLoaded);

		interpolate = 0.f;
		meshPtr->m_vertexPositions = meshPtr->m_vertexOriginPositions;
		meshPtr->m_triangleIndices = meshPtr->m_triangleOriginIndices;
		meshPtr->vertex_CurrentGroupe = meshPtr->m_component_Triangle_Real_MacroSurface_Load;
		meshPtr->recomputePerVertexNormals();
		meshPtr->computePlanarParameterization();
		meshPtr->init(taille_radius, false);*/
	}
	Mode = Mode_;
	if (Mode == -3.f) {
		
		meshPtr->m_vertexPositions = meshPtr->m_CompleteDataInterpolationTranfer.first.second;
		meshPtr->m_triangleIndices = meshPtr->m_CompleteDataInterpolationTranfer.second.second;
		meshPtr->vertex_CurrentGroupe = meshPtr->m_tranferGruppe;
		meshPtr->m_vertexPositions_NEW = meshPtr->m_vPosition2;
		meshPtr->recomputePerVertexNormals();
		meshPtr->computePlanarParameterization();
		meshPtr->init(taille_radius, true);
	}
	if (Mode == 3.0f) {
		meshPtr->vertex_CurrentGroupe= meshPtr->trianglesCurrentGroupToVertexCurrentGroup(meshPtr->component_Triangle_Potentiel_MacroSurface);
		meshPtr->init(taille_radius, false);
	}
	if (Mode == 4.0f) {
		meshPtr->vertex_CurrentGroupe =meshPtr->trianglesCurrentGroupToVertexCurrentGroup(meshPtr->component_Triangle_Real_MacroSurface);
		meshPtr->init(taille_radius, false);
	}
	
	// Animate any entity of the program here
	static const float initialTime = currentTime;
	float dt = currentTime - initialTime;
	// <---- Update here what needs to be animated over time ---->
	shaderProgramPtr->use();

	glm::mat4 matrix = cameraPtr->computeViewMatrix();
	//peut faire varier l'intensité de la lumière
	float intensity = 30.0f;//50.0f * (0.5+cos(dt));
	// variation de la position, de l'angle du cone et de la direction
	
// faire en sorte que la source de lumière soit détaché de la caméra
	lightSource1.TranformPositionandOrientation(matrix);
//pousser les nouvelles caractéristiques au GPU
	
	shaderProgramPtr->set ("lightSource1.intensity", intensity);
	shaderProgramPtr -> set("Mode",Mode);
	shaderProgramPtr->set("lightSource1.position",lightSource1.getPosition() );
	shaderProgramPtr->set("lightSource1.direction", lightSource1.getDirection());
	shaderProgramPtr->set("lightSource1.coneAngle", lightSource1.getConeAngle());
	shaderProgramPtr->set("Interpolate", interpolate);
	shaderProgramPtr->set("Scale", Scale);

	

}

void usage (const char * command) {
	std::cerr << "Usage : " << command << " [<file.off>]" << std::endl;
	
	std::exit (EXIT_FAILURE);
}

int main (int argc, char ** argv) {
	if (argc > 2)
		usage (argv[0]);
	init (argc == 1 ? DEFAULT_MESH_FILENAME : argv[1]); // Your initialization code (user interface, OpenGL states, scene with geometry, material, lights, etc)
	while (!glfwWindowShouldClose (windowPtr)) {
		//SliderFloat3("Model Matrix Translation", &model_matrix_translation.x, 0.0f, 960.0f);
		
		update (static_cast<float> (glfwGetTime ()));
		glfwPollEvents();
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
		ImGui::Begin("GUI", &my_tool_active, ImGuiWindowFlags_MenuBar);
		
		if (ImGui::BeginMenuBar())
		{
			if (ImGui::BeginMenu("Mesh to Show"))
			{
				if (ImGui::MenuItem("robot", "")) { newMEshLoaded = 0; }
				if (ImGui::MenuItem("homer", "")) { newMEshLoaded = 1; }
				if (ImGui::MenuItem("pegaso", "")) { newMEshLoaded = 2; }
				if (ImGui::MenuItem("dancer2", "")) { newMEshLoaded = 3; }
				if (ImGui::MenuItem("bozbezbozzel", "")) { newMEshLoaded = 4; }
				if (ImGui::MenuItem("Close", "")) { my_tool_active = false; }
				ImGui::EndMenu();
			}
			ImGui::EndMenuBar();
		}
		ImGui::Text("This is some useful text.");               // Display some text (you can use a format strings too)
		ImGui::SliderFloat("Scale", &Scale, 0.0f, 1.0f);	
		ImGui::ColorEdit4("color", color);
		if (ImGui::Button("Restart")) {
			interpolate = 0.f;
			meshPtr->m_vertexPositions = meshPtr->m_vertexOriginPositions;
			meshPtr->m_triangleIndices = meshPtr->m_triangleOriginIndices;
			meshPtr->vertex_CurrentGroupe = meshPtr->m_component_Triangle_Real_MacroSurface_Load;
			meshPtr->recomputePerVertexNormals();
			meshPtr->computePlanarParameterization();
			meshPtr->init(taille_radius, false);
		};ImGui::SameLine();
		if (ImGui::Button("Decimate")) {
			interpolate = 0.f;
			meshPtr->m_vertexPositions = meshPtr->m_vertexDecimationPositions;
			meshPtr->m_triangleIndices = meshPtr->m_triangleDecimationIndices;
			meshPtr->recomputePerVertexNormals();
			meshPtr->computePlanarParameterization();
			meshPtr->init(taille_radius, false);
		};ImGui::SameLine();
		if (ImGui::Button("Voxelisation")) {
			interpolate = 0.f;
			meshPtr->m_vertexPositions = meshPtr->m_vertexVoxelPositions;
			meshPtr->m_triangleIndices = meshPtr->m_triangleVoxelIndices;
			meshPtr->recomputePerVertexNormals();
			meshPtr->computePlanarParameterization();
			meshPtr->init(taille_radius, false);
		};ImGui::SameLine();
		if (ImGui::Button("Loading DataInterpolation")) {
			
			interpolate = 1.f;
			meshPtr->m_vertexPositions = meshPtr->m_CompleteDataInterpolationTranfer.first.second;
			meshPtr->m_triangleIndices = meshPtr->m_CompleteDataInterpolationTranfer.second.second;
			meshPtr->vertex_CurrentGroupe = meshPtr->m_tranferGruppe;
			meshPtr->m_vertexPositions_NEW = meshPtr->m_vPosition2;
			meshPtr->recomputePerVertexNormals();
			meshPtr->computePlanarParameterization();
			meshPtr->init(taille_radius, true);
			/*
			meshPtr->m_vertexPositions = CompleteDataInterpolationTranfer.first.second;
			meshPtr->m_triangleIndices = CompleteDataInterpolationTranfer.second.second;
			meshPtr->vertex_CurrentGroupe = tranferGruppe;
			meshPtr->m_vertexPositions_NEW = vPosition2;
			meshPtr->recomputePerVertexNormals();
			meshPtr->computePlanarParameterization();
			meshPtr->init(taille_radius, true);*/
			/*meshPtr->m_vertexPositions = meshPtr->m_DataInterpolationVertexLoad;
			meshPtr->m_triangleIndices = meshPtr->m_DataInterpolationTrianglesLoad;
			meshPtr->vertex_CurrentGroupe = meshPtr->m_vertex_transferGroupe_load;
			meshPtr->m_vertexPositions_NEW = meshPtr->m_vertexPositions_NEWLoad;
			meshPtr->recomputePerVertexNormals();
			meshPtr->computePlanarParameterization();
			meshPtr->init(taille_radius, true);*/
		};ImGui::SameLine();
		
		ImGui::RadioButton("Mode Texture", &Mode_, 0); ImGui::SameLine();
		ImGui::RadioButton("Mode X Toon", &Mode_, 1); 
		ImGui::RadioButton("Mode Profondeur", &Mode_, 2); ImGui::SameLine();
		ImGui::RadioButton("Mode Potentiel Macro", &Mode_, 3); ImGui::SameLine();
		ImGui::RadioButton("Mode Real Macro", &Mode_, 4);ImGui::SameLine();
		
		test = ImGui::GetWindowPos();
		test2 = ImGui::GetWindowSize();
		ImGui::End();
		render ();
		glfwSwapBuffers (windowPtr);
		
		
		
		
	}
	clear ();
	std::cout << " > Quit" << std::endl;
	return EXIT_SUCCESS;
}
