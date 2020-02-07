// ----------------------------------------------
// Base code for practical computer graphics
// assignments.
//
// Copyright (C) 2018 Tamy Boubekeur
// All rights reserved.
// Code deeply change by Paul-Ernest Martin
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


static const std::string DEFAULT_MESH_FILENAME("Resources/mesh_collection/ogre.off");
float taille_radius = 1.5f; 
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

static float Scale=0.f;   // need for geometry shader to knwo the percent of the scale
static float interpolate = 0; //bool to the fragment shader to knwo if you need to interpolate

static float color[4] = { 0.5f,0.5f,0.5f,0.5f }; // background color
static int Mode_ = 0; //mode to render (can be change in the menu imgui)

//Rendering mode (0 : PBR, 1 : toon shading, 2 : x-toon shading)
static float Mode = 0.f;

void clear ();




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

	// canceled -> Use it when you load the attributes from the files 
	else if (action == GLFW_PRESS && key == GLFW_KEY_N) {
		
		/*//les 3 modes possibles
		Mode = Mode + 1.0f;
		if (Mode ==  5.0f) {
			Mode = 0.0f;
		}
		if (Mode == 3.0f) {
			meshPtr->vertex_CurrentGroupe = meshPtr->m_component_Triangle_Potentiel_MacroSurface_Load;
			meshPtr->init(taille_radius,false);
		}
		if (Mode == 4.0f) {
			meshPtr->vertex_CurrentGroupe = meshPtr->m_component_Triangle_Real_MacroSurface_Load;
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
			meshPtr->vertex_CurrentGroupe = meshPtr->m_component_Triangle_Potentiel_MacroSurface_Load;
			meshPtr->init(taille_radius,false);
		}
		if (Mode == 4.0f) {
			meshPtr->vertex_CurrentGroupe = meshPtr->m_component_Triangle_Real_MacroSurface_Load;
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
		meshPtr->vertex_CurrentGroupe = meshPtr->m_component_Triangle_Potentiel_MacroSurface_Load;
		meshPtr->recomputePerVertexNormals();
		meshPtr->init(taille_radius,false);
		
	}
	else if (action == GLFW_PRESS && key == GLFW_KEY_D) {
		std::cout << "D" << std::endl;
		meshPtr->m_vertexPositions =meshPtr->m_vertexDecimationPositionsLoad;
		meshPtr->m_triangleIndices = meshPtr->m_triangleDecimationIndicesLoad;
		meshPtr->recomputePerVertexNormals();
		meshPtr->init(taille_radius,false);
		
	}
	else if (action == GLFW_PRESS && key == GLFW_KEY_V) {
		std::cout << "V"<<std::endl;
		meshPtr->m_vertexPositions = meshPtr->m_vertexVoxelPositionsLoad;
		meshPtr->m_triangleIndices = meshPtr->m_triangleVoxelIndicesLoad;
		meshPtr->recomputePerVertexNormals();
		meshPtr->computePlanarParameterization();
		meshPtr->init(taille_radius,false);
	}

	else if (action == GLFW_PRESS && key == GLFW_KEY_P) {
		std::cout << "P" << std::endl;
		meshPtr->m_vertexPositions = meshPtr->m_DataInterpolationVertexLoad;
		meshPtr->m_triangleIndices = meshPtr->m_DataInterpolationTrianglesLoad;
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
		meshPtr->m_vertexPositions = meshPtr->m_DataInterpolationVertexLoad;
		meshPtr->m_triangleIndices = meshPtr->m_DataInterpolationTrianglesLoad;
		meshPtr->vertex_CurrentGroupe = meshPtr->m_vertex_transferGroupe_load;
		meshPtr->m_vertexPositions_NEW = meshPtr->m_vertexPositions_NEWLoad;
		meshPtr->recomputePerVertexNormals();
		meshPtr->computePlanarParameterization();
		meshPtr->init(taille_radius,true);*/
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

#pragma endregion
	pair<glm::vec3, float> tmp_cen;
	tmp_cen=meshPtr->analyseBasicGeoStat( taille_radius);
	center = tmp_cen.first;
	meshScale = tmp_cen.second;

#pragma region HalfEdges

	
	meshPtr->computeTriangleAdjEdges();

	
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
	
#pragma region BufferInitialisation
	//************Initiation du Mesh pour préparer les buffers par défault**************//
	meshPtr->init(taille_radius,false);
	//************Initiation du Mesh pour préparer les buffers par défault**************/
#pragma endregion 
	//Pour permettre de reset le mesh : stockage des données d'origine dans un attribut de Mesh
	meshPtr->m_vertexOriginPositions = meshPtr->m_vertexPositions;
	meshPtr->m_triangleOriginIndices = meshPtr->m_triangleIndices;
	
#pragma region Intersection

	//if you want to check if the intersection fonction work

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
	std::cout << "p1= [1.0f, 1.0f, 1.0f]  |" << " p2= [0.5f, 0.0f, 0.0f]   |" << " R = 1   |" << " Intersection : " << i6 <<  " position du centre r=[-0.8f, 0.0f, 0.0f]"<<std::endl;*/
#pragma region
	
#pragma region Decimation
	
	
	meshPtr->m_DecimationPairPositionIndices =MeshDecimation::QEMDecimateForASpecifiedListOfEdges(meshPtr->m_vertexPositions, meshPtr->m_triangleIndices, taille_radius,meshPtr);
	meshPtr->m_vertexDecimationPositions = meshPtr->m_DecimationPairPositionIndices.first;
	meshPtr->m_triangleDecimationIndices = meshPtr->m_DecimationPairPositionIndices.second;
	
	
#pragma endregion 

#pragma region Voxelisation

	//*******************************Voxelisation******************************************//
	meshPtr->m_VoxelPairPositionIndices = init_Voxelisation(meshPtr->m_vertexPositions, meshPtr->m_triangleIndices, taille_radius);
	meshPtr->m_vertexVoxelPositions = meshPtr->m_VoxelPairPositionIndices.first;
	meshPtr->m_triangleVoxelIndices = meshPtr->m_VoxelPairPositionIndices.second;
	
#pragma endregion 

#pragma region CreationOfFirstStepScaleFixed

	meshPtr->Make_New_List_V2(taille_radius, Triangles, meshPtr);

#pragma endregion

#pragma region WrittingInDocuments	

	//MeshLoader::WriteAllDataToFiles("robot", meshPtr, taille_radius);

	//MeshLoader::LoadAllDataFromFiles("robot", meshPtr);

	//MeshLoader::ShowAllDataLoaded(meshPtr);*/
	
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
	//zMin and zMax is needed for the toon shading
	shaderProgramPtr->set ("z_Min", meshScale);
	shaderProgramPtr->set ("z_Max", meshScale*3);
	// to know if the interpolation need to be done
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

	// background in fonction of the mode
	
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
	
	float intensity = 30.0f;//50.0f * (0.5+cos(dt));

	
// lights not related to the camera
	lightSource1.TranformPositionandOrientation(matrix);

	
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
		
		
		update (static_cast<float> (glfwGetTime ()));
		glfwPollEvents();
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		ImGui::Begin("GUI", &my_tool_active, ImGuiWindowFlags_MenuBar);
		
		if (ImGui::BeginMenuBar())
		{
			if (ImGui::BeginMenu("Mesh to Show"))
			{  /*
				if (ImGui::MenuItem("robot", "")) { newMEshLoaded = 0; }
				if (ImGui::MenuItem("homer", "")) { newMEshLoaded = 1; }
				if (ImGui::MenuItem("pegaso", "")) { newMEshLoaded = 2; }
				if (ImGui::MenuItem("dancer2", "")) { newMEshLoaded = 3; }
				if (ImGui::MenuItem("bozbezbozzel", "")) { newMEshLoaded = 4; }
				if (ImGui::MenuItem("Close", "")) { my_tool_active = false; }*/
				ImGui::EndMenu();
			}
			ImGui::EndMenuBar();
		}
		ImGui::Text("This is Computer Science");               
		ImGui::SliderFloat("Scale", &Scale, 0.0f, 1.0f);	
		ImGui::ColorEdit4("color", color);
		if (ImGui::Button("Restart")) {
			interpolate = 0.f;
			meshPtr->m_vertexPositions = meshPtr->m_vertexOriginPositions;
			meshPtr->m_triangleIndices = meshPtr->m_triangleOriginIndices;
			meshPtr->vertex_CurrentGroupe = meshPtr->trianglesCurrentGroupToVertexCurrentGroup(meshPtr->component_Triangle_Potentiel_MacroSurface);
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
