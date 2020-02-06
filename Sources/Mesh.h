#ifndef MESH_H
#define MESH_H

#include <glad/glad.h>
#include <vector>
#include <memory>

#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <queue>
using namespace std;
#include "Transform.h"


class Mesh : public Transform {
public:
	virtual ~Mesh ();
	
	#pragma region  List of attributes needed to display the simple Mesh and fonction to enable loading from files whith MeshLaoder object
	//List of attibutes needed to display the Mesh 
	inline const std::vector<glm::vec3> & vertexPositions () const { return m_vertexPositions; }
	inline std::vector<glm::vec3> & vertexPositions () { return m_vertexPositions; }
	inline const std::vector<glm::vec3> & vertexNormals () const { return m_vertexNormals; }
	inline std::vector<glm::vec3> & vertexNormals () { return m_vertexNormals; }
	inline const std::vector<glm::vec2> & vertexTexCoords () const { return m_vertexTexCoords; }
	inline std::vector<glm::vec2> & vertexTexCoords () { return m_vertexTexCoords; }
	inline const std::vector<glm::uvec3> & triangleIndices () const { return m_triangleIndices; }
	inline std::vector<glm::uvec3> & triangleIndices () { return m_triangleIndices; }

	// attribute to store the vertex origin position and triangles origin indices (need for reset)
	std::vector<glm::vec3> m_vertexOriginPositions;
	std::vector<glm::uvec3> m_triangleOriginIndices;
	//attribute to store the vertex position and triangles indices -> this attribute will be displayed 
	std::vector<glm::vec3> m_vertexPositions;
	std::vector<glm::uvec3> m_triangleIndices;
	//attribute to store the vertex noramls and vertex uv coordinates -> this attribute will be used for displaying  
	std::vector<glm::vec3> m_vertexNormals;
	std::vector<glm::vec2> m_vertexTexCoords;
	#pragma endregion

	#pragma region  List of attributes needed to displaying the mesh analysed just with loading data from files and fonction to enable loading it from files whith MeshLaoder object
	

	//To see each fonction of this attribute check the attribute with the same name without _load at the end at the section [List of attributes needed to store data during the mesh analysis] in Mesh.h files

	std::vector<float> m_vertex_transferGroupe_load;
	std::vector<float> m_component_Triangle_Potentiel_MacroSurface_Load;
	std::vector<float> m_component_Triangle_Real_MacroSurface_Load;
	std::vector<glm::vec3> m_DataInterpolationVertexLoad;
	std::vector<glm::uvec3> m_DataInterpolationTrianglesLoad;
	std::vector<glm::vec3> m_vertexPositions_NEWLoad;
	std::vector<glm::vec3> m_vertexDecimationPositionsLoad;
	std::vector<glm::uvec3> m_triangleDecimationIndicesLoad;
	std::vector<glm::vec3> m_vertexVoxelPositionsLoad;
	std::vector<glm::uvec3> m_triangleVoxelIndicesLoad;
	float m_taille_radius_load;
	inline const std::vector<float> & vertex_transferGroupe_load () const { return m_vertex_transferGroupe_load; }
	inline std::vector<float> & vertex_transferGroupe_load() { return m_vertex_transferGroupe_load; }
	inline const std::vector<glm::vec3> & DataInterpolationVertexLoad () const { return m_DataInterpolationVertexLoad; }
	inline  std::vector<glm::vec3> & DataInterpolationVertexLoad() { return m_DataInterpolationVertexLoad; }
	inline const std::vector<glm::uvec3> & DataInterpolationTrianglesLoad () const { return m_DataInterpolationTrianglesLoad; }
	inline  std::vector<glm::uvec3> & DataInterpolationTrianglesLoad() { return m_DataInterpolationTrianglesLoad; }
	inline const std::vector<glm::vec3> & vertexPositions_NEWLoad () const { return m_vertexPositions_NEWLoad; }
	inline std::vector<glm::vec3> & vertexPositions_NEWLoad()  { return m_vertexPositions_NEWLoad; }
	inline const std::vector<float> & component_Triangle_Potentiel_MacroSurface_Load() const { return m_component_Triangle_Potentiel_MacroSurface_Load; }
	inline  std::vector<float> & component_Triangle_Potentiel_MacroSurface_Load() { return m_component_Triangle_Potentiel_MacroSurface_Load; }
	inline const std::vector<float> & component_Triangle_Real_MacroSurface_Load() const { return m_component_Triangle_Real_MacroSurface_Load; }
	inline std::vector<float> & component_Triangle_Real_MacroSurface_Load() { return m_component_Triangle_Real_MacroSurface_Load; }
	inline const std::vector<glm::vec3> & vertexDecimationPositionsLoad() const { return m_vertexDecimationPositionsLoad; }
	inline  std::vector<glm::vec3> & vertexDecimationPositionsLoad() { return m_vertexDecimationPositionsLoad; }
	inline const std::vector<glm::uvec3> & triangleDecimationIndicesLoad() const { return m_triangleDecimationIndicesLoad; }
	inline  std::vector<glm::uvec3> & triangleDecimationIndicesLoad() { return m_triangleDecimationIndicesLoad; }
	inline const std::vector<glm::vec3> & vertexVoxelPositionsLoad() const { return m_vertexVoxelPositionsLoad; }
	inline  std::vector<glm::vec3> & vertexVoxelPositionsLoad() { return m_vertexVoxelPositionsLoad; }
	inline const std::vector<glm::uvec3> & triangleVoxelIndicesLoad() const { return m_triangleVoxelIndicesLoad; }
	inline  std::vector<glm::uvec3> & triangleVoxelIndicesLoad() { return m_triangleVoxelIndicesLoad; }

#pragma endregion

	#pragma region List of functions needed to simple analysis and render and to init buffer 

	// Compute the parameters of a sphere which bounds the mesh
	void computeBoundingSphere (glm::vec3 & center, float & radius) const;
	// Enable to have the principle caratéristic of the mesh : center, radius... fill in the pair and some explicit display in the console
	pair<glm::vec3, float> analyseBasicGeoStat();
	//Name of function is enough explicit
	void recomputePerVertexNormals(bool angleBased = false);
	//Prepare all the buffer for the GPU, the bool enables the use of the new position of vertex after the decimation : we need to send two vertex positions object
	void init(float taille_radius, bool dblMesh);
	//Name of function is enough explicit
	void render();
	//Name of function is enough explicit
	void clear();
	//Name of function is enough explicit
	void computePlanarParameterization();

# pragma endregion

	//Fill triangle_Adj_Edges attribute with the adjacent triangles indices of each edges (-1 if doesn't exist);
	void computeTriangleAdjEdges();
	// attribute to store the list of adjacents triangle
	std::vector<glm::vec3> triangle_Adj_Edges;

	#pragma region List of fonction needed to check for intersection between figures
	
	//Check if triangle intersect the sphere
	bool triangleIntersectSphere(int absTriangle, glm::vec3 r, float R);
	//Name of function is enough explicit
	bool intersectOrinsideSphere(glm::vec3 x1, glm::vec3 x2, glm::vec3 r, float R);
	//Name of function is enough explicit
	bool intersectSphere(glm::vec3 p1, glm::vec3 p2, glm::vec3 r, float R);

	#pragma endregion

	#pragma region List of functions for the Macro Surface Analysis (included functions of exploration and related triangles count)

	//Check if the triangle is a potentiel macro surface 
	bool IsOnMacroSurface(int absTriangle, float radiusSphere);
	//Init the attribute component_Triangle_Potentiel_MacroSurface;
	void init_component_Triangle_Potentiel_MacroSurfaces(std::vector<int> listTriangle, float taille_radius);
	//Check if the potentiel macro surfaces is a real macro surface (need to have init component Triangle Potentiel MAcro Surfaces attribute
	void Potentiel_To_Real_Macro_Surfaces(float taille_radius);
	//Name of function is enough explicit
	bool IsPotentialMacroSurfacesMacroSurfaces(std::vector<int> Triangles, float radius);	

	//To tranfer triangles current groupe to vertex 
	std::vector<float> trianglesCurrentGroupToVertexCurrentGroup(std::pair<std::vector<int>, std::vector<int>> component_Triangle);

	
	//Check if the list of triangles has just one componnent realated  (need to have init attribute  triangle_Adj_Edges)
	bool hasSingleConnectComponent(std::vector<int> Triangles);
	// Enable the exploration and fill the triangle which have to explore 
	void Explore2(int currentTri,  queue<int> trisToExplore, float radiusSphere, glm::vec3 centerPosisiton);

	// Two attributes need for exploration when we deal With component comnnexe
	std::vector<int> I;
	queue<int> trisToExplore;

# pragma endregion

	#pragma region List of functions needed during mesh analysis whcich is not directed related to a particular section
		// Compute the nombre of edges to decime to respect the size, the list of string is unique -> "vertex indices lower" + " " + "vertex indices higher"
		pair<int,std::vector<string>> WichEdgesDecimeToRespectSizeScale(float taille_radius, std::vector<glm::uvec3> TrianglesIndices, std::vector<glm::vec3> VertexPositions);
		//attribute to store the list od edgest to decim
		std::vector<string> listEdges;

		//Check if element is in Queue
		bool containsQueue(queue<int> Queue, int element);
		// To check if the three vectors is in the list of Vertex and return the index find
		pair<bool[3],int[3]> containsVector(std::vector<glm::vec3> listVertex, glm::vec3 elements_0, glm::vec3 elements_1, glm::vec3 elements_2);
	# pragma endregion

	#pragma region	List of attributes to store the data index for all the buffer needed
	GLuint m_vao = 0;
	GLuint m_posVbo = 0;
	GLuint m_normalVbo = 0;
	GLuint m_texCoordVbo = 0;
	GLuint m_ibo = 0;
	GLuint m_connect = 0;
	GLuint m_posVbo_2 = 0;
	#pragma endregion

	#pragma region List of attributes need to store data from decimation and voxelisation to display (Section of code that is independent from the mesh analysis)
	std::vector<glm::vec3> m_vertexDecimationPositions;
	std::vector<glm::uvec3> m_triangleDecimationIndices;
	std::vector<glm::vec3> m_vertexVoxelPositions;
	std::vector<glm::uvec3> m_triangleVoxelIndices;

#pragma endregion
	
	#pragma region List of attributes needed to store data during the mesh analysis 

		// attribute to store the new position of vertex after the decimation  
		std::vector<glm::vec3> m_vertexPositions_NEW;

		// attribute to store the state of the vertex which enables the display [red/green for Potentiel Macro/MicroSurfaces] OR [red/green for Real Macro/MicroSurfaces] OR [state required for Interpolation (new index OR -1  OR -2)] 
		//-> linked with the gui interface (bool interpolate/ float Mode)
		std::vector<float> vertex_CurrentGroupe;

		// attribute to store the information for the potentiel macro surfaces triangles [List index Potentiel Macro Surfaces / List index Potentiel Micro Surfaces]
		std::pair<std::vector<int>, std::vector<int>> component_Triangle_Potentiel_MacroSurface;

		// attribute to store the information for the real macro surfaces triangles [List index Real Macro Surfaces / List index Real Micro Surfaces]
		std::pair<std::vector<int>, std::vector<int>> component_Triangle_Real_MacroSurface;

		// attribute to store vertex positions and triangles index of real macro surfaces and micro surfaces
		std::vector<glm::vec3> real_Macro_List_Vertex;
		std::vector<glm::uvec3> real_Macro_List_Triangles;
		std::vector<glm::vec3> real_Micro_List_Vertex;
		std::vector<glm::uvec3> real_Micro_List_Triangles;

		// attribute to store data after decimation on macro surface and voxelisation on micro surface
		pair< std::vector<glm::vec3>, std::vector<glm::uvec3>> m_DecimationPairPositionIndices;
		pair< std::vector<glm::vec3>, std::vector<glm::uvec3>> m_VoxelPairPositionIndices;

		// attribute to store the data  m_DecimationPairPositionIndices / m_VoxelPairPositionIndices in one pair with coherent index
		pair< pair<int, std::vector<glm::vec3>>, pair<int, std::vector<glm::uvec3>>> m_DataStoreAfterVoxDec;

		//attribute to store all the data in one attribute to enable interpolation with coherent index : Old Macro Vertex positions / Old Micro Vertex Positions / New Micro Vertex
		pair< pair<int, std::vector<glm::vec3>>, pair<int, std::vector<glm::uvec3>>> m_CompleteDataInterpolationTranfer;

		// attribute to store the temporary index of the decime vertex -> need to be modify deeply to work
		std::vector<float> vertex_transferGroupe;

		//attribute to store the state of each vertex [int>=0 -> index of the new vertex positions ;int==-1 vertex doesn't change in this step ; int==-1,5 old vertex belongs to micro surface ; int==-2 -> new vertex will belong to micro surface]
		std::vector<float> m_tranferGruppe;

		//Needed to send this to the buffer
		std::vector<glm::vec3> m_vPosition2;

	#pragma endregion 

	#pragma region List of functions that is used during the mesh analysis

		// The function that provide the complete analysis of the mesh : use the function in the right order to complete data analysis
		void Make_New_List_V2(float taille_radius, std::vector<int> TrianglesIndices,std::shared_ptr<Mesh> meshPtr);
	
		// From component_Triangles_Real_Macro_Surfaces to four attributes List real_[Macro/Micro]_List_Vertex + real_[Macro/Micro]_List_Triangles
		void Real_To_Two_List_Triangle_Position_Indices_V2();

		// Fill of two Pair <glm::vec3 List Position Vertex ; glm:uvec3 List Triangles Indices> : 
		//            - m_DecimationPairPositionIndices after decimation of the real_Macro_List_[Vertex/Tirangles]    + fill m_tranferGroupe for preservation of index during decimation
		//            - m_VoxelPairPositionIndices after voxelisation to the real_Micro_List_[Vertex/Tirangles]
		// 
		void Two_List_Triangles_Positions_Indices_To_New_List_Triangles_Positions_Indices_V2(float taille_radius, std::shared_ptr<Mesh> meshPtr);
	
		//From m_DecimationPairPositionIndices/m_VoxelPairPositionIndices to one attribute m_ CompleteDataInterpolatioTransfer whiwh index coherent
		void New_Two_List_To_One_V2();

		// From vertex_tranferGroupe (from QEM) which have the index tranfer during decimation to m_tranferGruppe 
		//  m_tranferGruppe -> [int>=0 -> index of the new vertex positions ;int==-1 vertex doesn't change in this step ; int==-1,5 old vertex belongs to micro surface ; int==-2 -> new vertex will belong to micro surface]
		// Fill m_vertexPositions_NEW to know the new positions of vertex
		void CompleteDataInterpolation_V2(float taille_radius);

		//modify m_tranferGruppe to arrange the index (if A became B and B became C, we need to put index of A became C directly + solve index changement due to disparition of index during edges collapse)
		// Finish to fill m_vertexPositions_NEW
		void CompleteTranferGruppe_V2(float taille_radius);

	#pragma endregion

};

#endif // MESH_H
