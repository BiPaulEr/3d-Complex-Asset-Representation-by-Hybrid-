#define _USE_MATH_DEFINES

#include "Mesh.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <cstdlib>
#include <string>
#include <queue>
#include "mdmeshDecimator.h"
#include "voxelizer.h"
using namespace std;

// In mesh.h explanations of each function is provided



Mesh::~Mesh () {
	clear ();
}
void Mesh::computeBoundingSphere (glm::vec3 & center, float & radius) const {
	center = glm::vec3 (0.0);
	radius = 0.f;
	for (const auto & p : m_vertexPositions)
		center += p;
	center /= m_vertexPositions.size ();
	for (const auto & p : m_vertexPositions)
		radius = std::max (radius, distance (center, p));
}
void Mesh::recomputePerVertexNormals (bool angleBased) {
	m_vertexNormals.clear ();
	// Change the following code to compute a proper per-vertex normal
	float agl0;
	float agl1;
	float agl2;
	m_vertexNormals.resize (m_vertexPositions.size (), glm::vec3 (0.0, 0.0, 0.0));

	glm::vec3 p0;
	glm::vec3 p1;
	glm::vec3 p2;
	
	glm::vec3 normal;
	std::vector<float> vertexAnglesSum;
	vertexAnglesSum.resize(m_vertexPositions.size(), 0.0);
	std::vector<glm::vec3> vertexNormSum;
	vertexNormSum.resize(m_vertexPositions.size (), glm::vec3 (0.0, 0.0, 0.0));

	//pour l'intégralité des triangles faire :

	for (int i = 0; i < m_triangleIndices.size (); i++) {
	//	std::cout <<"["<< m_triangleIndices[i][0] << "/" << m_triangleIndices[i][1] << "/" << m_triangleIndices[i][2] << "] " << std::flush;
		p0 = m_vertexPositions[m_triangleIndices[i][0]];
		p1 = m_vertexPositions[m_triangleIndices[i][1]];
		p2 = m_vertexPositions[m_triangleIndices[i][2]];

		if (angleBased) {
			agl0 = acos(dot(normalize(p1-p0), normalize(p2-p0)));
			agl1 = acos(dot(normalize(p2-p1), normalize(p0-p1)));
			agl2 = acos(dot(normalize(p0-p2), normalize(p1-p2)));
		}// si angleBased est faux on fait juste une moyenne des normes
		else {
			agl0 = 1;
			agl1 = 1;
			agl2 = 1;
		}

		normal = normalize(cross(p1-p0, p2-p0));
		// somme pondérées par les angles des normales pour les trois sommets
		vertexNormSum[m_triangleIndices[i][0]] += agl0*normal;
		vertexNormSum[m_triangleIndices[i][1]] += agl1*normal;
		vertexNormSum[m_triangleIndices[i][2]] += agl2*normal;
		// somme les angles pour les trois sommets
		vertexAnglesSum[m_triangleIndices[i][0]] += agl0;
		vertexAnglesSum[m_triangleIndices[i][1]] += agl1;
		vertexAnglesSum[m_triangleIndices[i][2]] += agl2;
	}

	for (int i = 0; i < m_vertexNormals.size (); i++) {
		//appplique la formule de pondération
		m_vertexNormals[i] = normalize(vertexNormSum[i]/vertexAnglesSum[i]);

	}
}
void Mesh::init (float taille_reduce,bool dblMesh) {
	
	float translation = taille_reduce* 10.0f;
	m_vertexPositions.push_back(glm::vec3(0, 0, translation));
	m_vertexPositions.push_back(glm::vec3(0, taille_reduce, translation));
	m_vertexPositions.push_back(glm::vec3(0, 0, taille_reduce+ translation));
										  
	m_vertexPositions.push_back(glm::vec3(0, 0, translation));
	m_vertexPositions.push_back(glm::vec3(taille_reduce,0, translation));
	m_vertexPositions.push_back(glm::vec3(0, 0, translation+taille_reduce));

	m_vertexNormals.push_back(glm::vec3(0, 0, 0));
	m_vertexNormals.push_back(glm::vec3(0, 0, 0));
	m_vertexNormals.push_back(glm::vec3(0, 0, 0));
	m_vertexNormals.push_back(glm::vec3(0, 0, 0));
	m_vertexNormals.push_back(glm::vec3(0, 0, 0));
	m_vertexNormals.push_back(glm::vec3(0, 0, 0));

	m_vertexTexCoords.push_back(glm::vec3(0, 0, 0));
	m_vertexTexCoords.push_back(glm::vec3(0, 0, 0));
	m_vertexTexCoords.push_back(glm::vec3(0, 0, 0));
	m_vertexTexCoords.push_back(glm::vec3(0, 0, 0));
	m_vertexTexCoords.push_back(glm::vec3(0, 0, 0));
	m_vertexTexCoords.push_back(glm::vec3(0, 0, 0));

	m_triangleIndices.push_back(glm::vec3(m_vertexPositions.size() - 1, m_vertexPositions.size() - 2, m_vertexPositions.size() - 3));
	m_triangleIndices.push_back(glm::vec3(m_vertexPositions.size() - 3, m_vertexPositions.size() - 4, m_vertexPositions.size() - 5));

	vertex_CurrentGroupe.push_back(-1);
	vertex_CurrentGroupe.push_back(-1);
	vertex_CurrentGroupe.push_back(-1);
	vertex_CurrentGroupe.push_back(-1);
	vertex_CurrentGroupe.push_back(-1);
	vertex_CurrentGroupe.push_back(-1);

	if (!dblMesh) {
		m_vertexPositions_NEW.resize(m_vertexPositions.size(), glm::vec3(0, 0, 0));
	}
	else {
		m_vertexPositions_NEW.push_back(glm::vec3(0, 0, 0));
		m_vertexPositions_NEW.push_back(glm::vec3(0, 0, 0));
		m_vertexPositions_NEW.push_back(glm::vec3(0, 0, 0));
		m_vertexPositions_NEW.push_back(glm::vec3(0, 0, 0));
		m_vertexPositions_NEW.push_back(glm::vec3(0, 0, 0));
		m_vertexPositions_NEW.push_back(glm::vec3(0, 0, 0));
	}
	glCreateBuffers (1, &m_posVbo); // Generate a GPU buffer to store the positions of the vertices
	size_t vertexBufferSize = sizeof (glm::vec3) * m_vertexPositions.size (); // Gather the size of the buffer from the CPU-side vector
	glNamedBufferStorage (m_posVbo, vertexBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT); // Create a data store on the GPU
	glNamedBufferSubData (m_posVbo, 0, vertexBufferSize, m_vertexPositions.data ()); // Fill the data store from a CPU array

	glCreateBuffers (1, &m_normalVbo); // Same for normal
	glNamedBufferStorage (m_normalVbo, vertexBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT);
	glNamedBufferSubData (m_normalVbo, 0, vertexBufferSize, m_vertexNormals.data ());

	glCreateBuffers (1, &m_texCoordVbo); // Same for texture coordinates
	size_t texCoordBufferSize = sizeof (glm::vec2) * m_vertexTexCoords.size ();
	glNamedBufferStorage (m_texCoordVbo, texCoordBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT);
	glNamedBufferSubData (m_texCoordVbo, 0, texCoordBufferSize, m_vertexTexCoords.data ());

	glCreateBuffers (1, &m_ibo); // Same for the index buffer, that stores the list of indices of the triangles forming the mesh
	size_t indexBufferSize = sizeof (glm::uvec3) * m_triangleIndices.size ();
	glNamedBufferStorage (m_ibo, indexBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT);
	glNamedBufferSubData (m_ibo, 0, indexBufferSize, m_triangleIndices.data ());

	
	//Component use to show the current Group of the triangles
	glCreateBuffers(1, &m_connect); // to show the component of each traingles
	size_t componentBufferSize = sizeof(float) *vertex_CurrentGroupe.size();
	glNamedBufferStorage(m_connect, componentBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT);
	glNamedBufferSubData(m_connect, 0, componentBufferSize, vertex_CurrentGroupe.data ());

	//Deuxieme pos
	glCreateBuffers(1, &m_posVbo_2);
	size_t vertex_2_BufferSize = sizeof(glm::vec3) * m_vertexPositions_NEW.size();
	glNamedBufferStorage(m_posVbo_2, vertex_2_BufferSize, NULL, GL_DYNAMIC_STORAGE_BIT);
	glNamedBufferSubData(m_posVbo_2, 0, vertex_2_BufferSize, m_vertexPositions_NEW.data());

	glCreateVertexArrays (1, &m_vao); // Create a single handle that joins together attributes (vertex positions, normals) and connectivity (triangles indices)
	glBindVertexArray (m_vao);

	glEnableVertexAttribArray (0);
	glBindBuffer (GL_ARRAY_BUFFER, m_posVbo);
	glVertexAttribPointer (0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof (GLfloat), 0);
	glEnableVertexAttribArray (1);
	glBindBuffer (GL_ARRAY_BUFFER, m_normalVbo);
	glVertexAttribPointer (1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof (GLfloat), 0);
	glEnableVertexAttribArray (2);
	glBindBuffer (GL_ARRAY_BUFFER, m_texCoordVbo);
	glVertexAttribPointer (2, 2, GL_FLOAT, GL_FALSE, 2 * sizeof (GLfloat), 0);
	

	//Component
	
	glEnableVertexAttribArray(3);
	glBindBuffer(GL_ARRAY_BUFFER, m_connect);
	glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE,  sizeof(GLfloat), 0);
	//2e pos
	glEnableVertexAttribArray(4);
	glBindBuffer(GL_ARRAY_BUFFER, m_posVbo_2);
	glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo);
	glBindVertexArray (0); // Desactive the VAO just created. Will be activated at rendering time.
}
void Mesh::render () {
	glBindVertexArray (m_vao); // Activate the VAO storing geometry data
	glDrawElements (GL_TRIANGLES, static_cast<GLsizei> (m_triangleIndices.size () * 3), GL_UNSIGNED_INT, 0); // Call for rendering: stream the current GPU geometry through the current GPU program
}
void Mesh::clear () {
	m_vertexPositions.clear ();
	m_vertexNormals.clear ();
	m_vertexTexCoords.clear ();
	m_triangleIndices.clear ();
	if (m_vao) {
		glDeleteVertexArrays (1, &m_vao);
		m_vao = 0;
	}
	if(m_posVbo) {
		glDeleteBuffers (1, &m_posVbo);
		m_posVbo = 0;
	}
	if (m_normalVbo) {
		glDeleteBuffers (1, &m_normalVbo);
		m_normalVbo = 0;
	}
	if (m_texCoordVbo) {
		glDeleteBuffers (1, &m_texCoordVbo);
		m_texCoordVbo = 0;
	}
	if (m_ibo) {
		glDeleteBuffers (1, &m_ibo);
		m_ibo = 0;
	}
}
void Mesh::computePlanarParameterization() {
	m_vertexTexCoords.clear();
	m_vertexTexCoords.resize(m_vertexPositions.size(), glm::vec2(0.f, 0.f));
	//coordonnéees des indices trangles des meshs 
	float Min_x = m_vertexPositions[0][0];
	float Max_x = m_vertexPositions[0][0];
	float Min_y = m_vertexPositions[0][1];
	float Max_y = m_vertexPositions[0][1];
	//recherche maximum minimum
	for (int i = 0; i < m_vertexPositions.size();i++) {
		if (Min_x > m_vertexPositions[i][0]) {
			Min_x = m_vertexPositions[i][0];
		}
		if (Max_x < m_vertexPositions[i][0]) {
			Max_x = m_vertexPositions[i][0];
		}
		if (Min_y > m_vertexPositions[i][1]) {
			Min_y = m_vertexPositions[i][1];
		}
		if (Max_y < m_vertexPositions[i][1]) {
			Max_y = m_vertexPositions[i][1];
		}
	}
	for (int i = 0; i < m_vertexPositions.size();i++) {
		float x = m_vertexPositions[i][0];
		float y = m_vertexPositions[i][1];
		x = (x - Min_x) / (Max_x - Min_x);
		y = (y - Min_y) / (Max_y - Min_y);
		glm::vec2 vec = glm::vec2(x, y);
		//associe les coordonnées dans la texture à chaque sommet 
		m_vertexTexCoords[i] = vec;
	}

}
pair<glm::vec3, float> Mesh::analyseBasicGeoStat() {
	float Min_dist = 1000.f;
	float Max_dist = 0.f;
	float sum_Dist = 0.f;
	float meshScale;
	glm::vec3 center = glm::vec3(0, 0, 0);
	int cmp = 0;
	glm::vec3 p0, p1, p2;
	float dist0, dist1, dist2;
	for (int i = 0; i < m_vertexPositions.size(); i++) {
		center = center + m_vertexPositions[i];
	}
	center /= m_vertexPositions.size();
	for (int i = 0; i < m_vertexPositions.size(); i++) {
		meshScale = std::max(meshScale, distance(m_vertexPositions[i], center));
	}
	for (int i = 0; i < m_triangleIndices.size(); i++) {
		p0 = m_vertexPositions[m_triangleIndices[i][0]];
		p1 = m_vertexPositions[m_triangleIndices[i][1]];
		p2 = m_vertexPositions[m_triangleIndices[i][2]];
		dist0 = glm::distance(p0, p1);
		dist1 = glm::distance(p2, p1);
		dist2 = glm::distance(p0, p2);
		cmp = cmp + 3;
		sum_Dist += dist0 + dist1 + dist2;
		if (dist0 > Max_dist) {
			Max_dist = dist0;
		}
		if (dist1 > Max_dist) {
			Max_dist = dist1;
		}
		if (dist2 > Max_dist) {
			Max_dist = dist2;
		}
		if (dist0 < Min_dist) {
			Min_dist = dist0;
		}
		if (dist1 < Min_dist) {
			Max_dist = dist1;
		}
		if (dist2 < Min_dist) {
			Min_dist = dist2;
		}
	}
	sum_Dist = sum_Dist / cmp;
	std::cout << "Max dist = " << Max_dist << " Min dist = " << Min_dist << " Avarage = " << sum_Dist << " center = [" << center[0] << "/" << center[1] << "/" << center[2] << "]" << " MesScale = " << meshScale << std::endl;
	pair<glm::vec3, float> result;
	result.first = center;
	result.second = meshScale;
	return result;
}



void Mesh::computeTriangleAdjEdges() {
	int i1, i2, i3;
	int x1, x2, x3;
	int t1, t2, t3;
	std::cout << "Begin Compute Triangle Edges : " << std::endl;
	std::cout << "[" << std::flush;
	int tracker = m_triangleIndices.size() / 20;
	triangle_Adj_Edges.resize(m_triangleIndices.size());
	for (int i = 0; i < m_triangleIndices.size();i++) {
		if (i%tracker == 0) {
			std::cout << "-" << std::flush;
		}
		x1 = m_triangleIndices[i][0];
		x2 = m_triangleIndices[i][1];
		x3 = m_triangleIndices[i][2];
		triangle_Adj_Edges[i][0] = -1;
		triangle_Adj_Edges[i][1] = -1;
		triangle_Adj_Edges[i][2] = -1;
		for (int t = 0; t < m_triangleIndices.size(); t++) {
			if (t != i) {
				t1 = m_triangleIndices[t][0];
				t2 = m_triangleIndices[t][1];
				t3 = m_triangleIndices[t][2];
				if ((x1 == t1 || x1 == t2 || x1 == t3) && (x2 == t1 || x2 == t2 || x2 == t3)) {
					triangle_Adj_Edges[i][0] = t;
				}
				if ((x2 == t1 || x2 == t2 || x2 == t3) && (x3 == t1 || x3 == t2 || x3 == t3)) {
					triangle_Adj_Edges[i][1] = t;
				}
				if ((x1 == t1 || x1 == t2 || x1 == t3) && (x3 == t1 || x3 == t2 || x3 == t3)) {
					triangle_Adj_Edges[i][2] = t;
				}

			}
		}

	}
	std::cout << "]" << std::endl;
	std::cout << "End Compute Triangle Edges" << std::endl;

}



bool Mesh::intersectOrinsideSphere(glm::vec3 p1, glm::vec3 p2, glm::vec3 r, float R) {
	float a = pow((p2[0] - p1[0]), 2.0f) + pow((p2[1] - p1[1]), 2.0f) + pow((p2[2] - p1[2]), 2.0f);
	float b = 2.0f * ((p2[0] - p1[0])*(p1[0] - r[0]) + (p2[1] - p1[1])*(p1[1] - r[1]) + (p2[2] - p1[2])*(p1[2] - r[2]));
	float c = pow(r[0], 2.0f) + pow(r[1], 2.0f) + pow(r[2], 2.0f) + pow(p1[0], 2.0f) + pow(p1[1], 2.0f) + pow(p1[2], 2.0f) - 2.0f* (r[0] * p1[0] + r[1] * p1[1] + r[2] * p1[2]) - pow(R, 2.0f);
	float delta = pow(b, 2.0f) - (4.0f*a*c);
	if (delta < 0.0f) {
		return false;
	}

	float u1 = (-b + pow(delta, 0.5f)) / (2.0f * a);
	float u2 = (-b - pow(delta, 0.5f)) / (2.0f * a);
	if (((u1 < 0.0f) && (u2 < 0.0f)) || ((u1 > 1.0f) && (u2 > 1.0f))) {
		return false;
	}
	else {
		return true;
	}
}
bool Mesh::intersectSphere(glm::vec3 p1, glm::vec3 p2, glm::vec3 r, float R) {
	float a = pow((p2[0] - p1[0]), 2.0f) + pow((p2[1] - p1[1]), 2.0f) + pow((p2[2] - p1[2]), 2.0f);
	float b = 2.0f * ((p2[0] - p1[0])*(p1[0] - r[0]) + (p2[1] - p1[1])*(p1[1] - r[1]) + (p2[2] - p1[2])*(p1[2] - r[2]));
	float c = pow(r[0], 2.0f) + pow(r[1], 2.0f) + pow(r[2], 2.0f) + pow(p1[0], 2.0f) + pow(p1[1], 2.0f) + pow(p1[2], 2.0f) - 2.0f* (r[0] * p1[0] + r[1] * p1[1] + r[2] * p1[2]) - pow(R, 2.0f);
	float delta = pow(b, 2.0f) - (4.0f*a*c);
	if (delta < 0.0f) {
		return false;
	}

	float u1 = (-b + pow(delta, 0.5f)) / (2.0f * a);
	float u2 = (-b - pow(delta, 0.5f)) / (2.0f * a);
	if (((u1 < 0.0f) && (u2 < 0.0f)) || ((u1 > 1.0f) && (u2 > 1.0f))) {
		return false;
	}
	if (((u1 < 0.0f) || (u2 < 0.0f)) && ((u1 > 1.0f) || (u2 > 1.0f))) {
		return false;
	}

	else

	{
		return true;
	}
}
bool Mesh::triangleIntersectSphere(int absTriangle, glm::vec3 r, float R) {
	glm::vec3 pos1 = m_vertexPositions[m_triangleIndices[absTriangle][0]];
	glm::vec3 pos2 = m_vertexPositions[m_triangleIndices[absTriangle][1]];
	glm::vec3 pos3 = m_vertexPositions[m_triangleIndices[absTriangle][2]];
	if (intersectSphere(pos1, pos2, r, R) || intersectSphere(pos2, pos3, r, R) || intersectSphere(pos3, pos1, r, R)) {
		return true;
	}
	else {
		return false;
	}
}




bool Mesh::IsPotentialMacroSurfacesMacroSurfaces(std::vector<int> Triangles, float radius) {
	glm::vec3 center = glm::vec3(0.f, 0.f, 0.f);

	//std::cout << "Position Vertex : " << std::endl;
	for (size_t i = 0; i < Triangles.size(); i++)
	{
		for (size_t t = 0; t < 3; t++)
		{
			//	std::cout << " [" << m_vertexPositions[m_triangleIndices[Triangles[i]][t]][0] << "," << m_vertexPositions[m_triangleIndices[Triangles[i]][t]][1] << "," << m_vertexPositions[m_triangleIndices[Triangles[i]][t]][1] << "] "<<std::flush;
			center += m_vertexPositions[m_triangleIndices[Triangles[i]][t]];
		}
	}
	//std::cout <<"/n" << "Position center : " << std::endl;
	center = (center /= (Triangles.size() * 3));
	//std::cout << "[" << center[0]<<","<< center[1]<<","<<center[2]<<"]"<< std::endl;
	for (size_t i = 0; i < Triangles.size(); i++)
	{
		for (size_t t = 0; t < 3; t++)
		{
			if (radius < distance(center, m_vertexPositions[m_triangleIndices[Triangles[i]][t]])) {
				//std::cout << "Resultat : " << true <<"\n" << std::endl;
				return true;

			}

		}
	}
	return false;
}
bool Mesh::IsOnMacroSurface(int absTriangle, float radiusSphere) {
	glm::vec3 pos1 = m_vertexPositions[m_triangleIndices[absTriangle][0]];
	glm::vec3 pos2 = m_vertexPositions[m_triangleIndices[absTriangle][1]];
	glm::vec3 pos3 = m_vertexPositions[m_triangleIndices[absTriangle][2]];
	glm::vec3 barycenter = (pos1 + pos2 + pos3);
	barycenter = barycenter / 3.0f;

	I.clear();
	while (!trisToExplore.empty()) {
		trisToExplore.pop();
	}

	int currentTri;
	trisToExplore.push(absTriangle);
	queue<int> exploredTris;
	while (!trisToExplore.empty()) {
		currentTri = trisToExplore.front();
		trisToExplore.pop();

		Mesh::Explore2(currentTri, exploredTris, radiusSphere, barycenter);
		exploredTris.push(currentTri);

	}
	//std::cout << "taille de I "<< I.size() << std::endl;
	return (hasSingleConnectComponent(I));
}
void Mesh::init_component_Triangle_Potentiel_MacroSurfaces(std::vector<int> listTriangle, float taille_radius) {
	std::cout << "[Potentiel Macro Surface Loading]" << std::endl;
	std::cout << "<                               >" << std::endl;
	size_t tracker = listTriangle.size() / 20;
	std::cout << "[" << std::flush;


	component_Triangle_Potentiel_MacroSurface.first = listTriangle;
	component_Triangle_Potentiel_MacroSurface.second.resize(listTriangle.size(), -1);
	bool temp;

	for (int i = 0; i < listTriangle.size();i++) {
		if (i % tracker == 0) {
			std::cout << "-" << std::flush;
		}
		temp = IsOnMacroSurface(listTriangle[i], taille_radius);
		if (temp == 0) {
			component_Triangle_Potentiel_MacroSurface.second[listTriangle[i]] = 0.f;
		}
		else {
			component_Triangle_Potentiel_MacroSurface.second[listTriangle[i]] = 1.f;
			//std::cout << "Macro surface detected" << std::endl;
		}

	}

	std::cout << "]" << std::endl;
	std::cout << "<                               >" << std::endl;
	std::cout << "[Potentiel Macro Surface Loaded]" << std::endl;
	std::cout << "" << std::endl;
	std::cout << "" << std::endl;
}
void Mesh::Potentiel_To_Real_Macro_Surfaces(float taille_radius) {
	std::cout << "[Real Macro Surface Loading] " << std::endl;
	std::cout << "<                          >" << std::endl;

	std::cout << "[" << std::flush;

	component_Triangle_Real_MacroSurface = component_Triangle_Potentiel_MacroSurface;
	std::vector<int> listTriangleToCheck;
	for (size_t i = 0; i < component_Triangle_Potentiel_MacroSurface.first.size(); i++)
	{

		if (component_Triangle_Potentiel_MacroSurface.second[i] == 1.f) {
			listTriangleToCheck.push_back(component_Triangle_Potentiel_MacroSurface.first[i]);
		}
	}
	size_t tracker = listTriangleToCheck.size() / 20;
	int recherche = 0;
	bool trouve = false;
	int currentGroup = 0;
	int t2;
	std::queue<int> Q;
	std::vector<int> listTriangleConnect;
	std::vector<int> comp(listTriangleToCheck.size(), -1);
	for (int t = 0; t < listTriangleToCheck.size();t++) {
		if (t % tracker == 0) {
			std::cout << "-" << std::flush;
		}
		if (comp[t] == -1) { //t pas encore visité
			listTriangleConnect.clear();
			Q.push(t);
			listTriangleConnect.push_back(listTriangleToCheck[t]);
			while (!Q.empty()) {
				t2 = Q.front();
				Q.pop();
				comp[t2] = currentGroup;
				for (int t3 = 0;t3 < 3;t3++) {
					if (triangle_Adj_Edges[listTriangleToCheck[t2]][t3] != -1) {
						// l'abj de Triangles ne correspond pas à l'indice de triangles_AdjEdges
						recherche = 0;
						trouve = false;
						while ((recherche < listTriangleToCheck.size()) && (!trouve)) {
							if (listTriangleToCheck[recherche] == triangle_Adj_Edges[listTriangleToCheck[t2]][t3]) {
								trouve = true;
							}
							else {
								recherche++;
							}
						}
						if (trouve) {
							if (comp[recherche] == -1) {
								comp[recherche] = currentGroup;
								Q.push(recherche);
								listTriangleConnect.push_back(listTriangleToCheck[recherche]);
							}
						}
					}
				}

			}

			currentGroup++;

			if (!IsPotentialMacroSurfacesMacroSurfaces(listTriangleConnect, taille_radius)) {

				for (size_t i = 0; i < listTriangleConnect.size(); i++)
				{
					component_Triangle_Real_MacroSurface.second[listTriangleConnect[i]] = 0.f;
				}
			}
		}


	}

	std::cout << "]" << std::endl;
	std::cout << "<                          >" << std::endl;
	std::cout << "[Real Macro Surface Loaded]" << std::endl;
	std::cout << "" << std::endl;
	std::cout << "" << std::endl;
}
bool Mesh::hasSingleConnectComponent(std::vector<int> Triangles) {
	int recherche = 0;
	bool trouve = false;
	int currentGroup = 0;
	int t2;
	std::queue<int> Q;
	std::vector<int> comp(Triangles.size(), -1);
	for (int t = 0; t < Triangles.size();t++){
		if (comp[t] == -1) { //t pas encore visité
			Q.push(t);
			while (!Q.empty()) {
				t2 = Q.front();
				Q.pop();
				comp[t2] = currentGroup;
				for (int t3 = 0;t3 < 3;t3++) {
					if (triangle_Adj_Edges[Triangles[t2]][t3] != -1) {			
						// l'abj de Triangles ne correspond pas à l'indice de triangles_AdjEdges
						recherche = 0;
						trouve = false;
						while ((recherche < Triangles.size()) && (!trouve)) {
							if (Triangles[recherche] == triangle_Adj_Edges[Triangles[t2]][t3]) {
								trouve = true;
							}
							else {
								recherche++;
							}
						}
						if (trouve) {
							if (comp[recherche] == -1) {	
								comp[recherche] = currentGroup;
								Q.push(recherche);
							}
						}
					}
				}
	
			}

			currentGroup++;
		}
		
	}
	// si on veut afficher les components connexes de tous les triangles
	/*component_Triangles = comp;
	std::cout << " Le nombre de component connexes trouve est de : " << currentGroup << std::endl;
	*/
	if (currentGroup == 1) {
		return true;
	}
	else {
		return false;
	}

}
std::vector<float> Mesh::trianglesCurrentGroupToVertexCurrentGroup(std::pair<std::vector<int>, std::vector<int>> component_Triangle) {
	std::vector<float> tmp;
	tmp.resize(m_vertexPositions.size(), -1);
	for (int i = 0; i < component_Triangle.first.size();i++) {
		float valueofComponent= (float)component_Triangle.second[i];
		tmp[m_triangleIndices[component_Triangle.first[i]][0]] = valueofComponent;
		tmp[m_triangleIndices[component_Triangle.first[i]][1]] = valueofComponent;
		tmp[m_triangleIndices[component_Triangle.first[i]][2]] = valueofComponent;
	}
	return tmp;
}
void Mesh::Explore2(int currentTri, queue<int> exploredTris, float radiusSphere, glm::vec3 centerPosition) {
	glm::vec3 pos1 = m_vertexPositions[m_triangleIndices[currentTri][0]];
	glm::vec3 pos2 = m_vertexPositions[m_triangleIndices[currentTri][1]];
	glm::vec3 pos3 = m_vertexPositions[m_triangleIndices[currentTri][2]];
	int adjTri;
	if (intersectOrinsideSphere(pos1, pos2, centerPosition, radiusSphere)) {
		adjTri = triangle_Adj_Edges[currentTri][0];
		if (adjTri != -1) {
			if (!containsQueue(exploredTris, adjTri) && !containsQueue(trisToExplore, adjTri)) {
				trisToExplore.push(adjTri);
			}
		}
	}
	if (intersectOrinsideSphere(pos2, pos3, centerPosition, radiusSphere)) {
		adjTri = triangle_Adj_Edges[currentTri][1];
		if (adjTri != -1) {
			if (!containsQueue(exploredTris, adjTri) && !containsQueue(trisToExplore, adjTri)) {
				trisToExplore.push(adjTri);
			}
		}
	}
	if (intersectOrinsideSphere(pos3, pos1, centerPosition, radiusSphere)) {
		adjTri = triangle_Adj_Edges[currentTri][2];
		if (adjTri != -1) {
			if (!containsQueue(exploredTris, adjTri) && !containsQueue(trisToExplore, adjTri)) {
				trisToExplore.push(adjTri);
			}
		}
	}
	if (triangleIntersectSphere(currentTri, centerPosition, radiusSphere)) {
		I.push_back(currentTri);
	}

}




bool Mesh::containsQueue(queue<int> Queue, int elements) {
	while (!Queue.empty()) {
		if (Queue.front() == elements) {
			return true;
					}
		else {
			Queue.pop();
		}
	}
	return false;
}
pair<bool[3], int[3]> Mesh::containsVector(std::vector<glm::vec3> listVertex, glm::vec3 elements_0, glm::vec3 elements_1, glm::vec3 elements_2) {
	int length = listVertex.size();
	int out = 0;
	glm::vec3 tmp;
	pair<bool[3], int[3]> resultats;
	resultats.first[0] = false;
	resultats.first[1] = false;
	resultats.first[2] = false;
	resultats.second[0] = -1;
	resultats.second[1] = -1;
	resultats.second[2] = -1;
	for (size_t i = 0; i < length; i++)
	{
		tmp = listVertex[i];
		if ((tmp[0] == elements_0[0]) && (tmp[1] == elements_0[1]) && (tmp[2] = elements_0[2])) {
			resultats.first[0] = true;
			out++;
			resultats.second[0] = i;
		}
		if ((tmp[0] == elements_1[0]) && (tmp[1] == elements_1[1]) && (tmp[2] = elements_1[2])) {
			resultats.first[1] = true;
			out++;
			resultats.second[1] = i;
		}
		if ((tmp[0] == elements_2[0]) && (tmp[1] == elements_2[1]) && (tmp[2] = elements_2[2])) {
			resultats.first[2] = true;
			out++;
			resultats.second[2] = i;
		}
		if (out > 2) {
			return resultats;
		}

	}
	return resultats;
}
pair<int,std::vector<string>> Mesh::WichEdgesDecimeToRespectSizeScale(float taille_radius,std::vector<glm::uvec3> TrianglesIndices,std::vector<glm::vec3> VertexPositions) {
	pair<int, std::vector<string>> result;
	string tmp;
	std::vector<string> listVerTexAllReadyCheck;
	int comptNumbreOfVertexToDecime=0;
	glm::vec3 X,Y,Z;
	int x, y, z,s;
	for (size_t i = 0; i < TrianglesIndices.size(); i++)
	{
		x = TrianglesIndices[i][0];
		y = TrianglesIndices[i][1];
		z = TrianglesIndices[i][2];
		if (x > y) {
			s = x;
			x = y;
			y = s;
		}
		if (y > z) {
			s = y;
			y = z;
			z = s;
			if (x > y) {
				s = x;
				x = y;
				y = s;
			}
		}
		
		X = VertexPositions[x];
		Y = VertexPositions[y];
		Z = VertexPositions[z];
	
		if (distance(X, Y) < taille_radius) {
			tmp = to_string(x) + " "+to_string(y);
			
			if (!(std::find(listVerTexAllReadyCheck.begin(), listVerTexAllReadyCheck.end(), tmp) != listVerTexAllReadyCheck.end())) {
				listVerTexAllReadyCheck.push_back(tmp);
				comptNumbreOfVertexToDecime++;
			}
		}
		if (distance(Y, Z) < taille_radius) {
			tmp = to_string(y) + " " + to_string(z);
			
			if (!(std::find(listVerTexAllReadyCheck.begin(), listVerTexAllReadyCheck.end(), tmp) != listVerTexAllReadyCheck.end())) {
				listVerTexAllReadyCheck.push_back(tmp);
				comptNumbreOfVertexToDecime++;
			}
		}
		if (distance(X, Z) < taille_radius) {
			tmp = to_string(x) + " " + to_string(z);
			if (!(std::find(listVerTexAllReadyCheck.begin(), listVerTexAllReadyCheck.end(), tmp) != listVerTexAllReadyCheck.end())) {
				listVerTexAllReadyCheck.push_back(tmp);
				comptNumbreOfVertexToDecime++;
			}
		}
	}
	result.second=listVerTexAllReadyCheck;
	result.first = comptNumbreOfVertexToDecime;
	return result;
}

#pragma region newlist
void Mesh::Real_To_Two_List_Triangle_Position_Indices_V2() {
	
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
	for (size_t i = 0; i < component_Triangle_Real_MacroSurface.first.size(); i++)
	{
		indiceTriangle = component_Triangle_Real_MacroSurface.first[i];
		vertexPosition_1 = m_vertexOriginPositions[m_triangleOriginIndices[indiceTriangle][0]];
		vertexPosition_2 = m_vertexOriginPositions[m_triangleOriginIndices[indiceTriangle][1]];
		vertexPosition_3 = m_vertexOriginPositions[m_triangleOriginIndices[indiceTriangle][2]];

		if (component_Triangle_Real_MacroSurface.second[i] == 1) {

			pair<bool[3], int[3]> rsult(containsVector(listMacroVertexPositions, vertexPosition_1, vertexPosition_2, vertexPosition_3));

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
			if ((indice_Triangle[0] < 0) || (indice_Triangle[1] < 0) || (indice_Triangle[2] < 0)) {
				std::cout << "Macro_inf_[" << indice_Triangle[0] << "," << indice_Triangle[1] << "," << indice_Triangle[2] << "] " << std::flush;
			}
		}
		else {
			pair<bool[3], int[3]> rsult(containsVector(listMicroVertexPositions, vertexPosition_1, vertexPosition_2, vertexPosition_3));
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

	std::cout << " [MacroPositions.size()/MacroTriangles.size()]=[" << listMacroVertexPositions.size() << " / " << listMacroTrianglesIndices.size() << "]" << std::endl;
	std::cout << "[MicroPositions.size() / MicroTriangles.size()] = [" << listMicroVertexPositions.size() << " / " << listMicroTrianglesIndices.size() << "]" << std::endl;

	
	real_Macro_List_Vertex = listMacroVertexPositions;
	real_Macro_List_Triangles = listMacroTrianglesIndices;
	real_Micro_List_Vertex = listMicroVertexPositions;
	real_Micro_List_Triangles = listMicroTrianglesIndices;
	std::cout << "<                               >" << std::endl;
	std::cout << "[Real to Two List Split Loaded]" << std::endl;
	std::cout << "" << std::endl;
	std::cout << "" << std::endl;

}

void Mesh::Two_List_Triangles_Positions_Indices_To_New_List_Triangles_Positions_Indices_V2(float taille_radius, std::shared_ptr<Mesh> meshPtr) {
	std::cout << "[Decimation And Voxelisation on the Two List Split Loading]" << std::endl;
	std::cout << "<                                                        >" << std::endl;


	//decimation
	std::cout << "" << std::endl;
	std::cout << "[---------------------------------------]" << std::endl;
	std::cout << "Begin  Decimation : " << "[" << real_Macro_List_Vertex.size() << "/" << real_Macro_List_Triangles.size() << "]" << std::endl;

	m_DecimationPairPositionIndices.first.clear();
	m_DecimationPairPositionIndices.second.clear();
	m_DecimationPairPositionIndices = MeshDecimation::QEMDecimateForASpecifiedListOfEdges(real_Macro_List_Vertex, real_Macro_List_Triangles, taille_radius, meshPtr);


	std::cout << "End Decimation  : " << "[" << m_DecimationPairPositionIndices.first.size() << "/" << m_DecimationPairPositionIndices.second.size() << "]" << std::endl;
	std::cout << "[---------------------------------------]" << std::endl;
	std::cout << "" << std::endl;
	

		//voxelisation
	std::cout << "[---------------------------------------]" << std::endl;
	m_VoxelPairPositionIndices.first.clear();
	m_VoxelPairPositionIndices.second.clear();
	std::cout << "Begin Voxelization : " << "[" << real_Micro_List_Vertex.size() << "/" << real_Micro_List_Triangles.size() << "]" << std::endl;
	m_VoxelPairPositionIndices = init_Voxelisation(real_Micro_List_Vertex, real_Micro_List_Triangles, taille_radius / 2);
	std::cout << "End Voxelization : " << "[" << m_VoxelPairPositionIndices.first.size() << "/" << m_VoxelPairPositionIndices.second.size() << "]" << std::endl;
	std::cout << "[---------------------------------------]" << std::endl;
	std::cout << "<                                                        >" << std::endl;
	std::cout << "[Decimation And Voxelisation on the Two List Split Loaded]" << std::endl;
	std::cout << "" << std::endl;
	std::cout << "" << std::endl;
}

void Mesh::New_Two_List_To_One_V2() {
	std::cout << "[New List Concatenated Loading]" << std::endl;
	std::cout << "<                             >" << std::endl;
	std::cout << " OldList VerTex nbr et Triangles [Macro/Micro] : " << "{[" << real_Macro_List_Vertex.size() << "/" << real_Micro_List_Vertex.size() << "],[" << real_Macro_List_Triangles.size() << "/" << real_Micro_List_Triangles.size() << "]}" << std::endl;
	
	
	int sizePos = m_DecimationPairPositionIndices.first.size();
	int sizeInd = m_DecimationPairPositionIndices.second.size();
	m_DataStoreAfterVoxDec.first.first = sizePos;
	m_DataStoreAfterVoxDec.second.first = sizeInd;
	for (size_t i = 0; i < sizePos; i++)
	{

		m_DataStoreAfterVoxDec.first.second.push_back(m_DecimationPairPositionIndices.first[i]);
	}
	for (size_t i = 0; i < m_VoxelPairPositionIndices.first.size(); i++)
	{
		m_DataStoreAfterVoxDec.first.second.push_back(m_VoxelPairPositionIndices.first[i]);

	}
	for (size_t i = 0; i < sizeInd; i++)
	{
		m_DataStoreAfterVoxDec.second.second.push_back(m_DecimationPairPositionIndices.second[i]);
		//std::cout << "[" << DataStoreAfterVoxDec.second.second[i][0] << "," << DataStoreAfterVoxDec.second.second[i][1] << "," << DataStoreAfterVoxDec.second.second[i][2] << "] " << std::flush;
	}
	for (size_t i = 0; i < m_VoxelPairPositionIndices.second.size(); i++)
	{
		m_DataStoreAfterVoxDec.second.second.push_back(m_VoxelPairPositionIndices.second[i] += sizePos);
	}
	m_DecimationPairPositionIndices.first.clear();
	m_DecimationPairPositionIndices.second.clear();
	m_VoxelPairPositionIndices.first.clear();
	m_VoxelPairPositionIndices.second.clear();
	std::cout << " NewList : " << "{[" << m_DataStoreAfterVoxDec.first.first << "->" << m_DataStoreAfterVoxDec.first.second.size() << "],[" << m_DataStoreAfterVoxDec.second.first << "->" << m_DataStoreAfterVoxDec.second.second.size() << "]}" << std::endl;
	std::cout << "<                             >" << std::endl;
	std::cout << "[New List Concatenated Loaded]" << std::endl;

}

void Mesh::Make_New_List_V2(float taille_radius, std::vector<int> TrianglesIndices, std::shared_ptr<Mesh> meshPtr) {

	
	std::cout << "[Potentiel Macro Surface Loading]" << std::endl;
	std::cout << "<                               >" << std::endl;
	init_component_Triangle_Potentiel_MacroSurfaces(TrianglesIndices, taille_radius);
	std::cout << "<                               >" << std::endl;
	std::cout << "[Potentiel Macro Surface Loaded]" << std::endl;
	std::cout << "" << std::endl;
	std::cout << "" << std::endl;


	std::cout << "[Real Macro Surface Loading] " << std::endl;
	std::cout << "<                          >" << std::endl;
	Potentiel_To_Real_Macro_Surfaces(taille_radius);
	std::cout << "<                          >" << std::endl;
	std::cout << "[Real Macro Surface Loaded]" << std::endl;
	std::cout << "" << std::endl;
	std::cout << "" << std::endl;

	std::cout << "[Real to Two List Split Loading]" << std::endl;
	std::cout << "<                              >" << std::endl;
	Real_To_Two_List_Triangle_Position_Indices_V2();
	std::cout << "<                               >" << std::endl;
	std::cout << "[Real to Two List Split Loaded]" << std::endl;
	std::cout << "" << std::endl;
	std::cout << "" << std::endl;


	std::cout << "[Decimation And Voxelisation on the Two List Split Loading]" << std::endl;
	std::cout << "<                                                        >" << std::endl;
	Two_List_Triangles_Positions_Indices_To_New_List_Triangles_Positions_Indices_V2(taille_radius,meshPtr);
	std::cout << "<                                                        >" << std::endl;
	std::cout << "[Decimation And Voxelisation on the Two List Split Loaded]" << std::endl;
	std::cout << "" << std::endl;
	std::cout << "" << std::endl;

	std::cout << "[New List Concatenated Loading]" << std::endl;
	std::cout << "<                             >" << std::endl;
	//std::cout << " OldList : " << "{[" << TemporaireStockage.first.first.size() << "/" << TemporaireStockage.second.first.size() << "],[" << TemporaireStockage.first.second.size() << "/" << TemporaireStockage.second.second.size() << "]}" << std::endl;
	New_Two_List_To_One_V2();
	std::cout << " NewList : " << "{[" << m_DataStoreAfterVoxDec.first.first << "->" << m_DataStoreAfterVoxDec.first.second.size() << "],[" << m_DataStoreAfterVoxDec.second.first << "->" << m_DataStoreAfterVoxDec.second.second.size() << "]}" << std::endl;
	std::cout << "<                             >" << std::endl;
	std::cout << "[New List Concatenated Loaded]" << std::endl;

	std::cout << "[Complete Data Interpolation Loading]" << std::endl;
	std::cout << "<                             >" << std::endl;
	CompleteDataInterpolation_V2(taille_radius);
	std::cout << "<                             >" << std::endl;
	std::cout << "[Complete Data Interpolation Loaded]" << std::endl;

	std::cout << "[Complete Data Interpolation Loading]" << std::endl;
	std::cout << "<                             >" << std::endl;
	CompleteTranferGruppe_V2(taille_radius);
	std::cout << "<                             >" << std::endl;
	std::cout << "[Complete Data Interpolation Loaded]" << std::endl;

}
void Mesh::CompleteDataInterpolation_V2( float taille_radius) {
	pair< pair<std::vector<glm::vec3>, std::vector<glm::uvec3>>, pair< std::vector<glm::vec3>, std::vector<glm::uvec3>>> TemporaireStockage;
	TemporaireStockage.first.first = real_Macro_List_Vertex;
	TemporaireStockage.first.second = real_Macro_List_Triangles;
	TemporaireStockage.second.first = real_Micro_List_Vertex;
	TemporaireStockage.second.second = real_Micro_List_Triangles;

	m_CompleteDataInterpolationTranfer.first.first = TemporaireStockage.first.first.size();
	m_CompleteDataInterpolationTranfer.first.second = TemporaireStockage.first.first;
	m_CompleteDataInterpolationTranfer.second.first = TemporaireStockage.first.second.size();
	m_CompleteDataInterpolationTranfer.second.second = TemporaireStockage.first.second;
	int size_old_d = m_CompleteDataInterpolationTranfer.first.second.size();

	for (size_t i = 0; i < TemporaireStockage.second.first.size(); i++)
	{
		m_CompleteDataInterpolationTranfer.first.second.push_back(TemporaireStockage.second.first[i]);
	}
	std::cout << "OldDEcim + OldVox : " << m_CompleteDataInterpolationTranfer.first.second.size() << std::endl;
	int size_old_d_v = m_CompleteDataInterpolationTranfer.first.second.size();
	glm::uvec3 indice_triangle_tmp;
	for (size_t i = 0; i < TemporaireStockage.second.second.size(); i++)
	{
		indice_triangle_tmp = TemporaireStockage.second.second[i];
		indice_triangle_tmp += size_old_d;
		m_CompleteDataInterpolationTranfer.second.second.push_back(indice_triangle_tmp);
		if ((indice_triangle_tmp[0] > size_old_d_v) || (indice_triangle_tmp[1] > size_old_d_v) || (indice_triangle_tmp[2] > size_old_d_v)) {
			std::cout << "OldDEcim + OldVox traingles : Max  " << std::endl;
		}
	}
	std::cout << "OldDEcim + OldVox traingles : " << m_CompleteDataInterpolationTranfer.second.second.size() << std::endl;


	for (size_t i = m_DataStoreAfterVoxDec.first.first; i < m_DataStoreAfterVoxDec.first.second.size(); i++)
	{
		m_CompleteDataInterpolationTranfer.first.second.push_back(m_DataStoreAfterVoxDec.first.second[i]);
	}
	int size_old_d_v_n_v = m_CompleteDataInterpolationTranfer.first.second.size();
	std::cout << "OldDecim + OldVox + New Vox : " << m_CompleteDataInterpolationTranfer.first.second.size() << std::endl;
	int size_new = m_DataStoreAfterVoxDec.first.first;

	for (size_t i = m_DataStoreAfterVoxDec.second.first; i < m_DataStoreAfterVoxDec.second.second.size(); i++)
	{
		indice_triangle_tmp = (m_DataStoreAfterVoxDec.second.second[i] += (size_old_d_v - size_new));
		m_CompleteDataInterpolationTranfer.second.second.push_back(indice_triangle_tmp);
		if ((indice_triangle_tmp[0] > size_old_d_v_n_v) || (indice_triangle_tmp[2] > size_old_d_v_n_v) || (indice_triangle_tmp[2] > size_old_d_v_n_v)) {
			std::cout << "OldDEcim + OldVox traingles + Newtriangles : Max  " << std::endl;
		}
		//std::cout << "[" << CompleteDataInterpolationTranfer.second.second[i][0] << "," << CompleteDataInterpolationTranfer.second.second[i][1] << "," << CompleteDataInterpolationTranfer.second.second[i][2] << "] " << std::flush;
	}
	std::cout << "OldDecim + OldVox + New Vox Triangle : " << m_CompleteDataInterpolationTranfer.second.second.size() << std::endl;
	std::cout << "\n" << m_CompleteDataInterpolationTranfer.first.second.size();
}
void Mesh::CompleteTranferGruppe_V2(float taille_radius) {
	pair< pair<std::vector<glm::vec3>, std::vector<glm::uvec3>>, pair< std::vector<glm::vec3>, std::vector<glm::uvec3>>> TemporaireStockage;
	TemporaireStockage.first.first = real_Macro_List_Vertex;
	TemporaireStockage.first.second = real_Macro_List_Triangles;
	TemporaireStockage.second.first = real_Micro_List_Vertex;
	TemporaireStockage.second.second = real_Micro_List_Triangles;
	std::cout << "vexter Gruppe : " << vertex_transferGroupe.size();
	int cmp = 0, cmp_ = 0;
	int final_, tmp_ = 0;
	m_tranferGruppe.resize(vertex_transferGroupe.size(), -1);
	std::cout << "[" << std::endl;
	std::vector<int> tmp_tmp;

	for (size_t i = 0; i < vertex_transferGroupe.size(); i++)
	{
		//std::cout << "[" << i << "->" << meshPtr->vertex_transferGroupe[i] << "] " << std::flush;
		tmp_ = vertex_transferGroupe[i];
		if (tmp_ != -1) {
			cmp_++;
			while (tmp_ != -1) {
				final_ = tmp_;
				tmp_ = vertex_transferGroupe[tmp_];
			}
			tmp_tmp.push_back(cmp_);
			m_tranferGruppe[i] = final_;

		}
		else {
			cmp++;
			tmp_tmp.push_back(cmp_);
		}


	}
	for (size_t i = 0; i < vertex_transferGroupe.size(); i++)
	{

		if (vertex_transferGroupe[i] != -1) {
			m_tranferGruppe[i] = m_tranferGruppe[i] - tmp_tmp[m_tranferGruppe[i]];
		}

	}
	std::cout << "]" << std::endl;
	std::cout << "nbr de -1 : " << cmp << "  " << std::endl;



	for (size_t i = vertex_transferGroupe.size(); i < m_CompleteDataInterpolationTranfer.first.second.size(); i++)
	{
		if (i < TemporaireStockage.first.first.size() + TemporaireStockage.second.first.size()) {
			m_tranferGruppe.push_back(-1.5f);
		}
		else {
			m_tranferGruppe.push_back(-2.f);
		}


	}


	m_vertexPositions_NEW.resize(m_CompleteDataInterpolationTranfer.first.second.size(), glm::vec3(0.f, 0.f, 0.f));
	for (size_t i = 0; i < TemporaireStockage.first.first.size(); i++)
	{
		if (m_tranferGruppe[i] > -0.5f) {
			m_vertexPositions_NEW[i] = m_DataStoreAfterVoxDec.first.second[m_tranferGruppe[i]];
		}
		else {
			//meshPtr->m_vertexPositions_NEW[i] = glm::vec3(0, 0, 0);
		}
	}
	m_vPosition2 = m_vertexPositions_NEW;
	TemporaireStockage.first.first.clear();
	TemporaireStockage.first.second.clear();
	TemporaireStockage.second.first.clear();
	TemporaireStockage.second.second.clear();
}
#pragma endregion
	