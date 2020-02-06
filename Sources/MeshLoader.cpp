#include "MeshLoader.h"

#include <iostream>
#include <fstream>
#include <exception>
#include <ios>

using namespace std;

void MeshLoader::loadOFF (const std::string & filename, std::shared_ptr<Mesh> meshPtr) {
	std::cout << " > Start loading mesh <" << filename << ">" << std::endl;
    meshPtr->clear ();
	ifstream in (filename.c_str ());
    if (!in)
        throw std::ios_base::failure ("[Mesh Loader][loadOFF] Cannot open " + filename);
	string offString;
    unsigned int sizeV, sizeT, tmp;
    in >> offString >> sizeV >> sizeT >> tmp;
	
    auto & P = meshPtr->vertexPositions ();
    auto & T = meshPtr->triangleIndices ();
    P.resize (sizeV);
    T.resize (sizeT);
    size_t tracker = (sizeV + sizeT)/20;
    std::cout << " > [" << std::flush;
    for (unsigned int i = 0; i < sizeV; i++) {
    	if (i % tracker == 0)
    		std::cout << "-" << std::flush;
        in >> P[i][0] >> P[i][1] >> P[i][2];
    }
    int s;
    for (unsigned int i = 0; i < sizeT; i++) {
    	if ((sizeV + i) % tracker == 0)
    		std::cout << "-" << std::flush;
        in >> s;
        for (unsigned int j = 0; j < 3; j++)
            in >> T[i][j];
    }
    std::cout << "]" << std::endl;
    in.close ();
    meshPtr->vertexNormals ().resize (P.size (), glm::vec3 (0.f, 0.f, 1.f));
    meshPtr->vertexTexCoords ().resize (P.size (), glm::vec2 (0.f, 0.f));
    meshPtr->recomputePerVertexNormals ();
		meshPtr->computePlanarParameterization();
    std::cout << " > Mesh <" << filename << "> loaded" <<  std::endl;
}


//Write List of indices of adjacent triangles in the file
void MeshLoader::writeEdges(const std::string & filename,std::shared_ptr<Mesh> meshPtr) {
	std::cout << " > Start writting edges <" << filename << ">" << std::endl;
	
	ofstream out(filename.c_str());
	if (!out)
		throw std::ios_base::failure("[Edges writer] Cannot open " + filename);
	size_t tracker = meshPtr->triangle_Adj_Edges.size()/ 10;
	std::cout << " > [" << std::flush;
	for (unsigned int i = 0; i < meshPtr->triangle_Adj_Edges.size(); i++) {
		if (i % tracker == 0)
			std::cout << "-" << std::flush;
		out << meshPtr->triangle_Adj_Edges[i][0] << " " << meshPtr->triangle_Adj_Edges[i][1] << " " << meshPtr->triangle_Adj_Edges[i][2] << endl;
	}
	
	std::cout << "]" << std::endl;
	out.close();
	
	std::cout << " > Edges <" << filename << "> writed" << std::endl;
}
// Load attribut triangle_Adj_Edges from mesh object
void MeshLoader::loadEdges(const std::string & filename, std::shared_ptr<Mesh> meshPtr) {
	std::cout << " > Start loading Edges <" << filename << ">" << std::endl;
	ifstream in(filename.c_str());
	if (!in)
		throw std::ios_base::failure("[Edges Loader][loadOFF] Cannot open " + filename);
	auto & P = meshPtr->triangle_Adj_Edges;
	int size = meshPtr->m_triangleIndices.size();
	std::cout << " size : " << size << std::endl;
	P.resize(size);
	
	size_t tracker = (P.size()) / 20;
	std::cout << " > [" << std::flush;
	for (unsigned int i = 0; i < size; i++) {
		if (i % tracker == 0)
			std::cout << "-" << std::flush;
		in >> P[i][0] >> P[i][1] >> P[i][2];
	}

	
	std::cout << "]" << std::endl;
	in.close();
	std::cout << " > Edges <" << filename << "> loaded" << std::endl;
}


//Write in files the ListofTriangles and Vertex from DataInterpolation, size_radius used to analyse the mesh, the transfer data between thow sizes, the new position of vertex after the decimation
void MeshLoader::writeScaleDataInterpolation(const std::string & filename,std::shared_ptr<Mesh> meshPtr,float &taille_radius) {
	std::cout << " > Start writting DataInterpolation <" << filename << ">" << std::endl;

	ofstream out(filename.c_str());
	if (!out)
		throw std::ios_base::failure("[Edges writer] Cannot open " + filename);
	size_t tracker = (meshPtr->m_CompleteDataInterpolationTranfer.first.second.size()+ meshPtr->m_CompleteDataInterpolationTranfer.second.second.size() )/ 10;
	out << taille_radius << " " << meshPtr->m_CompleteDataInterpolationTranfer.first.second.size() << " " << meshPtr->m_CompleteDataInterpolationTranfer.second.second.size() << std::endl;
	std::cout << " > [" << std::flush;
	for (unsigned int i = 0; i < meshPtr->m_CompleteDataInterpolationTranfer.first.second.size(); i++) {
		if (i % tracker == 0) {
			std::cout << "-" << std::flush;
		}
		out << meshPtr->m_CompleteDataInterpolationTranfer.first.second[i][0] << " " << meshPtr->m_CompleteDataInterpolationTranfer.first.second[i][1] << " " << meshPtr->m_CompleteDataInterpolationTranfer.first.second[i][2] << endl;
	}
	for (unsigned int i = 0; i < meshPtr->m_CompleteDataInterpolationTranfer.second.second.size(); i++) {
		if (i % tracker == 0) {
			std::cout << "-" << std::flush;
		}
		out << meshPtr->m_CompleteDataInterpolationTranfer.second.second[i][0] << " " << meshPtr->m_CompleteDataInterpolationTranfer.second.second[i][1] << " " << meshPtr->m_CompleteDataInterpolationTranfer.second.second[i][2] << endl;
	}
	for (unsigned int i = 0; i < meshPtr->m_vertexPositions_NEW.size(); i++) {
		if (i % tracker == 0) {
			std::cout << "-" << std::flush;
		}
		out << meshPtr->m_vertexPositions_NEW[i][0] << " " << meshPtr->m_vertexPositions_NEW[i][1] << " " << meshPtr->m_vertexPositions_NEW[i][2] << endl;
	}
	for (unsigned int i = 0; i < meshPtr->m_CompleteDataInterpolationTranfer.first.second.size(); i++) {
		if (i % tracker == 0) {
			std::cout << "-" << std::flush;
		}
		out << meshPtr->m_tranferGruppe[i] << std::flush;
	}
	
	std::cout << "]" << std::endl;
	out.close();

	std::cout << " > DataInterpolation <" << filename << "> writed" << std::endl;
}
//Load data from file to fill the attributes from mesh object : DataInterpolationVertexLoad/DataInterpolationTrianglesLoad/vertex_transferGroupe_load/vertexPositions_NEWLoad/size_radius
void MeshLoader::loadScaleDataInterpolation(const std::string & filename, std::shared_ptr<Mesh> meshPtr) {
	std::cout << " > Start loading DataInterpolation <" << filename << ">" << std::endl;
	
	ifstream in(filename.c_str());
	if (!in)
		throw std::ios_base::failure("[Mesh Loader][loadOFF] Cannot open " + filename);
	float taille_radius;
	int sizeV, sizeT;
	in >> taille_radius >> sizeV >> sizeT;
	auto & P = meshPtr->DataInterpolationVertexLoad ();
	auto & T = meshPtr->DataInterpolationTrianglesLoad ();
	auto & TG = meshPtr->vertex_transferGroupe_load ();
	auto & VN = meshPtr->vertexPositions_NEWLoad ();
	
	meshPtr->m_taille_radius_load = taille_radius;
	P.resize(sizeV);
	T.resize(sizeT);
	TG.resize(sizeV);
	VN.resize(sizeV);

	size_t tracker = (sizeV + sizeT) / 20;
	std::cout << " > [" << std::flush;
	for (unsigned int i = 0; i < sizeV; i++) {
		if (i % tracker == 0)
			std::cout << "-" << std::flush;
		in >> P[i][0] >> P[i][1] >> P[i][2];
	}
	int s;
	for (unsigned int i = 0; i < sizeT; i++) {
		if ((sizeV + i) % tracker == 0)
			std::cout << "-" << std::flush;
	
		in >> T[i][0] >> T[i][1] >> T[i][2];
	}
	for (unsigned int i = 0; i < sizeV; i++) {
		if (i % tracker == 0)
			std::cout << "-" << std::flush;
		in >> VN[i][0] >> VN[i][1] >> VN[i][2];
	}
	std::cout << " nbr Col " << sizeT + sizeV * 2<<std::endl;
	for (unsigned int i = 0; i < sizeV; i++) {
		if (i % tracker == 0)
			std::cout << "-" << std::flush;
		in >> TG[i] ;
		
	}
	std::cout << "]" << std::endl;
	in.close();
	meshPtr->vertexNormals().resize(P.size(), glm::vec3(0.f, 0.f, 1.f));
	meshPtr->vertexTexCoords().resize(P.size(), glm::vec2(0.f, 0.f));
	meshPtr->recomputePerVertexNormals();
	meshPtr->computePlanarParameterization();
	std::cout << " > Mesh <" << filename << "> loaded" << std::endl;
}
//Show Data from m_vertexPositions_NEWLoad/m_DataInterpolationTrianglesLoad/m_DataInterpolationVertexLoad/m_vertex_transferGroupe_load
void MeshLoader::showScaleDataInterpolation(std::shared_ptr<Mesh> meshPtr){
	std::cout << " Taille radius : " << meshPtr->m_taille_radius_load << std::endl;
	
	std::cout << " Taille Nbr [Vertex / Triangles]: [" << meshPtr->m_DataInterpolationVertexLoad.size()<< "/"<< meshPtr->m_DataInterpolationTrianglesLoad.size()<<"]"<<std::endl;
	/*std::cout << "List Vertex Position : " << std::endl;
	for (size_t i = 0; i < meshPtr->m_DataInterpolationVertexLoad.size(); i++)
	{
		std::cout << "[" << meshPtr->m_DataInterpolationVertexLoad[i][0] << "/" << meshPtr->m_DataInterpolationVertexLoad[i][1] << "/" << meshPtr->m_DataInterpolationVertexLoad[i][2] << "]" << " " << std::flush;
	}
	std::cout<< ""<<std::endl;
	std::cout << "List Triangles " << std::endl;
	for (size_t i = 0; i < meshPtr->m_DataInterpolationTrianglesLoad.size(); i++)
	{
		std::cout << "[" << meshPtr->m_DataInterpolationTrianglesLoad[i][0] << "/" << meshPtr->m_DataInterpolationTrianglesLoad[i][1] << "/" << meshPtr->m_DataInterpolationTrianglesLoad[i][2] << "]" << " " << std::flush;
	}
	std::cout << "" << std::endl;

	std::cout << "List Vertex New : " << meshPtr->m_vertexPositions_NEWLoad.size()<< std::endl;
	for (size_t i = 0; i < meshPtr->m_vertexPositions_NEWLoad.size(); i++)
	{
		std::cout << "[" << meshPtr->m_vertexPositions_NEWLoad[i][0] << "/" << meshPtr->m_vertexPositions_NEWLoad[i][1] << "/" << meshPtr->m_vertexPositions_NEWLoad[i][2] << "]" << " " << std::flush;
	}
	std::cout << "" << std::endl;*/
	std::cout << "List Tranfer Gruppe : " << meshPtr->m_vertex_transferGroupe_load.size()<<std::endl;
	for (size_t i = 0; i < meshPtr->m_vertex_transferGroupe_load.size()/100; i++)
	{
		std::cout  << meshPtr->m_vertex_transferGroupe_load[i] << "/"  << std::flush;
	}
	std::cout << "" << std::endl;
}


//Write  Potentiel and real macro surfaces current groupe object : Need to have component_Triangle_Potentiel_MacroSurface and component_Triangle_Real_MacroSurface initialise 
void MeshLoader::writePotentielAndRealMacrosurface(const std::string & filename, std::shared_ptr<Mesh> meshPtr) {
	std::cout << " > Start writting Potentiel and Real Macrosurfaces <" << filename << ">" << std::endl;
	std::vector<float> tmp= meshPtr->trianglesCurrentGroupToVertexCurrentGroup(meshPtr->component_Triangle_Potentiel_MacroSurface);
	ofstream out(filename.c_str());
	if (!out)
		throw std::ios_base::failure("[Edges writer] Cannot open " + filename);
	out << tmp.size() << std::endl;
	size_t tracker = 2*tmp.size() / 10;
	std::cout << " > [" << std::flush;
	for (unsigned int i = 0; i < tmp.size(); i++) {
		if (i % tracker == 0)
			std::cout << "-" << std::flush;
		out << tmp[i] << " " << flush;
	}
	out << std::endl;
	tmp = meshPtr->trianglesCurrentGroupToVertexCurrentGroup(meshPtr->component_Triangle_Real_MacroSurface);
	for (unsigned int i = 0; i < tmp.size(); i++) {
		if (i % tracker == 0)
			std::cout << "-" << std::flush;
		out << tmp[i] << " " << flush;
	}
	std::cout << "]" << std::endl;
	out.close();

	std::cout << " > Potentiel and Real Macrosurfaces <" << filename << "> writed" << std::endl;

}
//Load  component_Triangle_Potentiel_MacroSurface_Load and component_Triangle_Real_MacroSurface_Load which is currentGroupe of vertex
void MeshLoader::loadPotentielAndRealMacrosurface(const std::string & filename, std::shared_ptr<Mesh> meshPtr) {
	std::cout << " > Start loading Potentiel and Real Macro Surfaces <" << filename << ">" << std::endl;
	ifstream in(filename.c_str());
	if (!in)
		throw std::ios_base::failure("[Potentiel And Real Macro Surfaces Loader][loadOFF] Cannot open " + filename);
	auto & PMS = meshPtr->component_Triangle_Potentiel_MacroSurface_Load ();
	auto & RMS = meshPtr->component_Triangle_Real_MacroSurface_Load();
	int r;
	in >> r;
	PMS.resize(r);
	RMS.resize(r);
	size_t tracker = (PMS.size()) / 20;
	std::cout << " > [" << std::flush;
	for (unsigned int i = 0; i < PMS.size(); i++) {
		if (i % tracker == 0)
			std::cout << "-" << std::flush;
		in >> PMS[i];
	}
	
	for (unsigned int i = 0; i < PMS.size(); i++) {
		if (i % tracker == 0)
			std::cout << "-" << std::flush;
		in >> RMS[i];
	}
	std::cout << "]" << std::endl;
	in.close();
	std::cout << " > Potentiel And Real Macro Surfaces <" << filename << "> loaded" << std::endl;
}


// Write the attributes vertexDecimationPositions/ triangleDecimationIndices/ vertexVoxelPositions / triangleVoxelIndices from files to Mesh input
void MeshLoader::writeDecimationAndVoxel(const std::string & filename, std::shared_ptr<Mesh> meshPtr) {
	std::cout << " > Start writting Decimation And Voxel <" << filename << ">" << std::endl;

	ofstream out(filename.c_str());
	if (!out)
		throw std::ios_base::failure("[Decimation And Voxel writer] Cannot open " + filename);
	size_t tracker = (meshPtr->m_vertexDecimationPositions.size() + meshPtr->m_triangleDecimationIndices.size() + meshPtr->m_triangleVoxelIndices.size()+ meshPtr->m_vertexVoxelPositions.size()) / 10;
	out << meshPtr->m_vertexDecimationPositions.size() <<" "<< meshPtr->m_triangleDecimationIndices.size() <<" "<< meshPtr->m_vertexVoxelPositions.size() << " " << meshPtr->m_triangleVoxelIndices.size();
	std::cout << " > [" << std::flush;
	for (unsigned int i = 0; i < meshPtr->m_vertexDecimationPositions.size(); i++) {
		if (i % tracker == 0)
			std::cout << "-" << std::flush;
		out << meshPtr->m_vertexDecimationPositions[i][0] << " " << meshPtr->m_vertexDecimationPositions[i][1] << " " << meshPtr->m_vertexDecimationPositions[i][2] << endl;
	}
	for (unsigned int i = 0; i < meshPtr->m_triangleDecimationIndices.size(); i++) {
		if (i % tracker == 0)
			std::cout << "-" << std::flush;
		out << meshPtr->m_triangleDecimationIndices[i][0] << " " << meshPtr->m_triangleDecimationIndices[i][1] << " " << meshPtr->m_triangleDecimationIndices[i][2] << endl;
	}
	for (unsigned int i = 0; i < meshPtr->m_vertexVoxelPositionsLoad.size(); i++) {
		if (i % tracker == 0)
			std::cout << "-" << std::flush;
		out << meshPtr->m_vertexVoxelPositions[i][0] << " " << meshPtr->m_vertexVoxelPositions[i][1] << " " << meshPtr->m_vertexVoxelPositions[i][2] << endl;
	}
	for (unsigned int i = 0; i < meshPtr->m_triangleVoxelIndices.size(); i++) {
		if (i % tracker == 0)
			std::cout << "-" << std::flush;
		out << meshPtr->m_triangleVoxelIndices[i][0] << " " << meshPtr->m_triangleVoxelIndices[i][1] << " " << meshPtr->m_triangleVoxelIndices[i][2] << endl;
	}
	std::cout << "]" << std::endl;
	out.close();

	std::cout << " > Decimation and Voxel <" << filename << "> writed" << std::endl;
}
//Load attributes vertexDecimationPositionsLoad/ triangleDecimationIndicesLoad/vertexVoxelPositionsLoad/triangleVoxelIndicesLoad from files to Mesh input
void MeshLoader::loadDecimationAndVoxel(const std::string & filename, std::shared_ptr<Mesh> meshPtr) {
	std::cout << " > Start loading Decimation and Voxel <" << filename << ">" << std::endl;
	ifstream in(filename.c_str());
	if (!in)
		throw std::ios_base::failure("[Decimation and Voxel Loader][loadOFF] Cannot open " + filename);
	
	auto & VPD = meshPtr->vertexDecimationPositionsLoad();
	auto & ITD = meshPtr->triangleDecimationIndicesLoad();
	auto & VPV = meshPtr->vertexVoxelPositionsLoad();
	auto & ITV = meshPtr->triangleVoxelIndicesLoad();
	int vpd, itd, vpv, itv;
	in >> vpd >> itd >> vpv >> itv;
	VPD.resize(vpd);
	ITD.resize(itd);
	VPV.resize(vpv);
	ITV.resize(itv);
	size_t tracker = (VPD.size()+ITD.size()+VPV.size()+ITV.size()) / 20;
	std::cout << " > [" << std::flush;
	for (unsigned int i = 0; i < vpd; i++) {
		if (i % tracker == 0)
			std::cout << "-" << std::flush;
		in >> VPD[i][0] >> VPD[i][1] >> VPD[i][2];
	}

	for (unsigned int i = 0; i < itd; i++) {
		if (i % tracker == 0)
			std::cout << "-" << std::flush;
		in >> ITD[i][0] >> ITD[i][1] >> ITD[i][2];
	}
	for (unsigned int i = 0; i < vpv; i++) {
		if (i % tracker == 0)
			std::cout << "-" << std::flush;
		in >> VPV[i][0] >> VPV[i][1] >> VPV[i][2];
	}
	for (unsigned int i = 0; i < itv; i++) {
		if (i % tracker == 0)
			std::cout << "-" << std::flush;
		in >> ITV[i][0] >> ITV[i][1] >> ITV[i][2];
	}

	std::cout << "]" << std::endl;
	in.close();
	std::cout << " > Decimation and Voxel <" << filename << "> loaded" << std::endl;
}

void MeshLoader::WriteAllDataToFiles(std::string Name, std::shared_ptr<Mesh> meshPtr) {
	try {
		MeshLoader::writeEdges(Name + "_edges.off", meshPtr);
	}
	catch (std::exception & e) {
		//exitOnCriticalError(std::string("[Error loading edges]") + e.what());
	}
	try {
		MeshLoader::writePotentielAndRealMacrosurface(Name + "_PotentielAndRealSurfaces.off", meshPtr);
	}
	catch (std::exception & e) {
		//exitOnCriticalError(std::string("[Error loading edges]") + e.what());
	}
	try {
		MeshLoader::writeDecimationAndVoxel(Name + "_DecimationAndVoxel.off", meshPtr);
	}
	catch (std::exception & e) {
		//	exitOnCriticalError(std::string("[Error loading edges]") + e.what());
	}
	try {
		MeshLoader::writeScaleDataInterpolation(Name + "_DataInterpolation.off", meshPtr);
	}
	catch (std::exception & e) {
		//exitOnCriticalError(std::string("[Error loading edges]") + e.what());
	}

}