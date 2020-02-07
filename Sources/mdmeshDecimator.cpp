


#define _CRT_SECURE_NO_WARNINGS
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include "mdMeshDecimator.h"
#include <glm/glm.hpp>
#include "Mesh.h"
namespace MeshDecimation
{

	int MeshDecimator::GetTriangle(int v1, int v2, int v3) const
	{
		int i, j, k;
		int idTriangle;
		for (size_t it = 0; it < m_vertices[v1].m_triangles.Size(); ++it)
		{
			idTriangle = m_vertices[v1].m_triangles[it];
			i = m_triangles[idTriangle].X();
			j = m_triangles[idTriangle].Y();
			k = m_triangles[idTriangle].Z();
			if ((i == v1 && j == v2 && k == v3) || (i == v1 && j == v3 && k == v2) ||
				(i == v2 && j == v1 && k == v3) || (i == v2 && j == v3 && k == v1) ||
				(i == v3 && j == v2 && k == v1) || (i == v3 && j == v1 && k == v2))
			{
				return idTriangle;
			}
		}
		return -1;
	}
	int MeshDecimator::GetEdge(int v1, int v2) const
	{
		int idEdge;
		for (size_t it = 0; it < m_vertices[v1].m_edges.Size(); ++it)
		{
			idEdge = m_vertices[v1].m_edges[it];
			if ((m_edges[idEdge].m_v1 == v1 && m_edges[idEdge].m_v2 == v2) ||
				(m_edges[idEdge].m_v1 == v2 && m_edges[idEdge].m_v2 == v1))
			{
				return idEdge;
			}
		}
		return -1;
	}
	void MeshDecimator::EdgeCollapse(int v1, int v2)
	{
		transferInd[v2] = (float)v1;
		int u, w;
		int shift;
		int idTriangle;
		for (size_t itT = 0; itT < m_vertices[v2].m_triangles.Size(); ++itT)
		{
			idTriangle = m_vertices[v2].m_triangles[itT];
			if (m_triangles[idTriangle].X() == v2)
			{
				shift = 0;
				u = m_triangles[idTriangle].Y();
				w = m_triangles[idTriangle].Z();
			}
			else if (m_triangles[idTriangle].Y() == v2)
			{
				shift = 1;
				u = m_triangles[idTriangle].X();
				w = m_triangles[idTriangle].Z();
			}
			else
			{
				shift = 2;
				u = m_triangles[idTriangle].X();
				w = m_triangles[idTriangle].Y();
			}

			if ((u == v1) || (w == v1))
			{
				m_trianglesTags[idTriangle] = false;
				m_vertices[u].m_triangles.Erase(idTriangle);
				m_vertices[w].m_triangles.Erase(idTriangle);
				m_nTriangles--;
			}
			else if (GetTriangle(v1, u, w) == -1)
			{
				m_vertices[v1].m_triangles.Insert(idTriangle);
				m_triangles[idTriangle][shift] = v1;
			}
			else
			{
				m_trianglesTags[idTriangle] = false;
				m_vertices[u].m_triangles.Erase(idTriangle);
				m_vertices[w].m_triangles.Erase(idTriangle);
				m_nTriangles--;
			}
		}
		int idEdge;
		for (size_t itE = 0; itE < m_vertices[v2].m_edges.Size(); ++itE)
		{
			idEdge = m_vertices[v2].m_edges[itE];
			w = (m_edges[idEdge].m_v1 == v2) ? m_edges[idEdge].m_v2 : m_edges[idEdge].m_v1;
			if (w == v1)
			{
				m_edges[idEdge].m_tag = false;
				m_vertices[w].m_edges.Erase(idEdge);
				m_nEdges--;
			}
			else if (GetEdge(v1, w) == -1)
			{
				if (m_edges[idEdge].m_v1 == v2)    m_edges[idEdge].m_v1 = v1;
				else                            m_edges[idEdge].m_v2 = v1;
				m_vertices[v1].m_edges.Insert(idEdge);
			}
			else
			{
				m_edges[idEdge].m_tag = false;
				m_vertices[w].m_edges.Erase(idEdge);
				m_nEdges--;
			}
		}
		m_vertices[v2].m_tag = false;
		m_nVertices--;
		// update boundary edges
		SArray<int, 64> incidentVertices;
		incidentVertices.PushBack(v1);
		for (size_t itE = 0; itE < m_vertices[v1].m_edges.Size(); ++itE)
		{
			incidentVertices.PushBack((m_edges[idEdge].m_v1 != v1) ? m_edges[idEdge].m_v1 : m_edges[idEdge].m_v2);
			idEdge = m_vertices[v1].m_edges[itE];
			m_edges[idEdge].m_onBoundary = (IsBoundaryEdge(m_edges[idEdge].m_v1, m_edges[idEdge].m_v2) != -1);
		}
		// update boundary vertices
		int idVertex;
		for (size_t itV = 0; itV < incidentVertices.Size(); ++itV)
		{
			idVertex = incidentVertices[itV];
			m_vertices[idVertex].m_onBoundary = false;
			for (size_t itE = 0; itE < m_vertices[idVertex].m_edges.Size(); ++itE)
			{
				idEdge = m_vertices[idVertex].m_edges[itE];
				if (m_edges[idEdge].m_onBoundary)
				{
					m_vertices[idVertex].m_onBoundary = true;
					break;
				}
			}
		}
	}
	int MeshDecimator::IsBoundaryEdge(int v1, int v2) const
	{
		int commonTri = -1;
		int itTriangle1, itTriangle2;
		for (size_t itT1 = 0; itT1 < m_vertices[v1].m_triangles.Size(); ++itT1)
		{
			itTriangle1 = m_vertices[v1].m_triangles[itT1];
			for (size_t itT2 = 0; itT2 < m_vertices[v2].m_triangles.Size(); ++itT2)
			{
				itTriangle2 = m_vertices[v2].m_triangles[itT2];
				if (itTriangle1 == itTriangle2)
				{
					if (commonTri == -1)
					{
						commonTri = itTriangle1;
					}
					else
					{
						return -1;
					}
				}
			}
		}
		return commonTri;
	}
	bool MeshDecimator::IsBoundaryVertex(int v) const
	{
		int idEdge;
		for (size_t itE = 0; itE < m_vertices[v].m_edges.Size(); ++itE)
		{
			idEdge = m_vertices[v].m_edges[itE];
			if (IsBoundaryEdge(m_edges[idEdge].m_v1, m_edges[idEdge].m_v2) != -1) return true;
		}
		return false;
	}
	void MeshDecimator::GetMeshData(Vec3<Float> * points, Vec3<int> * triangles) const
	{
		int * map = new int[m_nPoints];
		int counter = 0;
		for (size_t v = 0; v < m_nPoints; ++v)
		{
			if (m_vertices[v].m_tag)
			{
				points[counter] = m_points[v];
				map[v] = counter++;
			}
		}
		counter = 0;
		for (size_t t = 0; t < m_nInitialTriangles; ++t)
		{
			if (m_trianglesTags[t])
			{
				triangles[counter].X() = map[m_triangles[t].X()];
				triangles[counter].Y() = map[m_triangles[t].Y()];
				triangles[counter].Z() = map[m_triangles[t].Z()];
				counter++;
			}
		}
		delete[] map;
	}	
	void MeshDecimator::InitializePriorityQueue()
	{
		double progressOld = -1.0;
		double progress = 0.0;
		char msg[1024];
		double ptgStep = 1.0;
		int v1, v2;
		MDEdgePriorityQueue pqEdge;
		size_t nE = m_edges.size();
		size_t tracker = nE / 100;
		std::cout << "\n" << "[" << std::flush;
		
		
		for (size_t e = 0; e < nE; ++e)
		{if ((e%tracker) == 0) {
			
			std::cout << "-" <<std::flush;
			
		}
			

			if (m_edges[e].m_tag)
			{
				v1 = m_edges[e].m_v1;
				v2 = m_edges[e].m_v2;
				if ((!m_ecolManifoldConstraint) || (ManifoldConstraint(v1, v2)))
				{
					pqEdge.m_qem = m_edges[e].m_qem = ComputeEdgeCost(v1, v2, m_edges[e].m_pos);
					pqEdge.m_name = static_cast<int>(e);
					m_pqueue.push(pqEdge);
				}
			}
		}
		std::cout << "]" << std::endl;
	}
	double MeshDecimator::ComputeEdgeCost(int v1, int v2, Vec3<Float> & newPos) const
	{
		double Q[10];
		double M[12];
		Vec3<double> pos;
		for (int i = 0; i < 10; ++i) Q[i] = m_vertices[v1].m_Q[i] + m_vertices[v2].m_Q[i];
		M[0] = Q[0]; // (0, 0)
		M[1] = Q[1]; // (0, 1) 
		M[2] = Q[2]; // (0, 2)
		M[3] = Q[3]; // (0, 3)
		M[4] = Q[1]; // (1, 0)
		M[5] = Q[4]; // (1, 1)
		M[6] = Q[5]; // (1, 2)
		M[7] = Q[6]; // (1, 3)
		M[8] = Q[2]; // (2, 0)
		M[9] = Q[5]; // (2, 1)
		M[10] = Q[7]; // (2, 2);
		M[11] = Q[8]; // (2, 3);
		double det = M[0] * M[5] * M[10] + M[1] * M[6] * M[8] + M[2] * M[4] * M[9]
			- M[0] * M[6] * M[9] - M[1] * M[4] * M[10] - M[2] * M[5] * M[8];
		if (det != 0.0)
		{
			double d = 1.0 / det;
			pos.X() = d * (M[1] * M[7] * M[10] + M[2] * M[5] * M[11] + M[3] * M[6] * M[9]
				- M[1] * M[6] * M[11] - M[2] * M[7] * M[9] - M[3] * M[5] * M[10]);
			pos.Y() = d * (M[0] * M[6] * M[11] + M[2] * M[7] * M[8] + M[3] * M[4] * M[10]
				- M[0] * M[7] * M[10] - M[2] * M[4] * M[11] - M[3] * M[6] * M[8]);
			pos.Z() = d * (M[0] * M[7] * M[9] + M[1] * M[4] * M[11] + M[3] * M[5] * M[8]
				- M[0] * M[5] * M[11] - M[1] * M[7] * M[8] - M[3] * M[4] * M[9]);
			newPos.X() = static_cast<Float>(pos.X());
			newPos.Y() = static_cast<Float>(pos.Y());
			newPos.Z() = static_cast<Float>(pos.Z());
		}
		else
		{
			const Float w = static_cast<Float>(0.5f);
			newPos = w * m_points[v1] + w * m_points[v2];
			pos.X() = static_cast<double>(newPos.X());
			pos.Y() = static_cast<double>(newPos.Y());
			pos.Z() = static_cast<double>(newPos.Z());
		}

		double qem = pos.X()  * (Q[0] * pos.X() + Q[1] * pos.Y() + Q[2] * pos.Z() + Q[3]) +
			pos.Y()  * (Q[1] * pos.X() + Q[4] * pos.Y() + Q[5] * pos.Z() + Q[6]) +
			pos.Z()  * (Q[2] * pos.X() + Q[5] * pos.Y() + Q[7] * pos.Z() + Q[8]) +
			(Q[3] * pos.X() + Q[6] * pos.Y() + Q[8] * pos.Z() + Q[9]);

		Vec3<Float> d1;
		Vec3<Float> d2;
		Vec3<Float> n1;
		Vec3<Float> n2;
		Vec3<Float> oldPosV1 = m_points[v1];
		Vec3<Float> oldPosV2 = m_points[v2];

		SArray<int, SARRAY_DEFAULT_MIN_SIZE> triangles = m_vertices[v1].m_triangles;
		int idTriangle;
		for (size_t itT = 0; itT < m_vertices[v2].m_triangles.Size(); ++itT)
		{
			idTriangle = m_vertices[v2].m_triangles[itT];
			triangles.Insert(idTriangle);
		}
		int a[3];
		for (size_t itT = 0; itT != triangles.Size(); ++itT)
		{
			idTriangle = triangles[itT];
			a[0] = m_triangles[idTriangle].X();
			a[1] = m_triangles[idTriangle].Y();
			a[2] = m_triangles[idTriangle].Z();

			d1 = m_points[a[1]] - m_points[a[0]];
			d2 = m_points[a[2]] - m_points[a[0]];
			n1 = d1 ^ d2;

			m_points[v1] = newPos;
			m_points[v2] = newPos;

			d1 = m_points[a[1]] - m_points[a[0]];
			d2 = m_points[a[2]] - m_points[a[0]];
			n2 = d1 ^ d2;

			m_points[v1] = oldPosV1;
			m_points[v2] = oldPosV2;

			n1.Normalize();
			n2.Normalize();
			if (n1*n2 < 0.0)
			{
				return std::numeric_limits<double>::max();
			}
		}
		if (m_ecolManifoldConstraint && !ManifoldConstraint(v1, v2))
		{
			return std::numeric_limits<double>::max();
		}
		return qem;
	}
	bool MeshDecimator::ManifoldConstraint(int v1, int v2) const
	{
		std::set<int> vertices;
		int a, b;
		int idEdge1;
		int idEdge2;
		int idEdgeV1V2;
		for (size_t itE1 = 0; itE1 < m_vertices[v1].m_edges.Size(); ++itE1)
		{
			idEdge1 = m_vertices[v1].m_edges[itE1];
			a = (m_edges[idEdge1].m_v1 == v1) ? m_edges[idEdge1].m_v2 : m_edges[idEdge1].m_v1;
			vertices.insert(a);
			if (a != v2)
			{
				for (size_t itE2 = 0; itE2 < m_vertices[v2].m_edges.Size(); ++itE2)
				{
					idEdge2 = m_vertices[v2].m_edges[itE2];
					b = (m_edges[idEdge2].m_v1 == v2) ? m_edges[idEdge2].m_v2 : m_edges[idEdge2].m_v1;
					vertices.insert(b);
					if (a == b)
					{
						if (GetTriangle(v1, v2, a) == -1)
						{
							return false;
						}
					}
				}
			}
			else
			{
				idEdgeV1V2 = idEdge1;
			}
		}
		if (vertices.size() <= 4 || (m_vertices[v1].m_onBoundary && m_vertices[v2].m_onBoundary && !m_edges[idEdgeV1V2].m_onBoundary))
		{
			return false;
		}
		return true;
	}


	// function that have been modified 
	MeshDecimator::MeshDecimator(void)
	{
		m_triangles = 0;
		m_points = 0;
		m_nPoints = 0;
		m_nInitialTriangles = 0;
		m_nVertices = 0;
		m_nTriangles = 0;
		m_nEdges = 0;
		m_trianglesTags = 0;
		m_ecolManifoldConstraint = true;
		m_callBack = 0;
		m_taille_radius = 0;
		targetError = 0;
	}
	MeshDecimator::~MeshDecimator(void)
	{
		ReleaseMemory();
	}
	void MeshDecimator::ReleaseMemory()
	{
		delete[] m_trianglesTags;
		std::vector< MDVertex > emptyVertices(0);
		m_vertices.swap(emptyVertices);
		std::vector<MDEdge> emptyEdges(0);
		m_edges.swap(emptyEdges);
		m_pqueue = std::priority_queue<
			MDEdgePriorityQueue,
			std::vector<MDEdgePriorityQueue>,
			std::less<MDEdgePriorityQueue> >();
		m_triangles = 0;
		m_points = 0;
		m_nPoints = 0;
		m_nInitialTriangles = 0;
		m_nVertices = 0;
		m_nTriangles = 0;
		m_nEdges = 0;
		m_trianglesTags = 0;
		transferInd.clear();
		m_taille_radius = 0;
	}
	void MeshDecimator::Initialize(size_t nVertices, Float taille_radius, size_t nTriangles, Vec3<Float> * points, Vec3<int> * triangles)
	{
		m_nVertices = nVertices;
		m_nTriangles = nTriangles;
		m_nInitialTriangles = nTriangles;
		m_points = points;
		m_nPoints = nVertices;
		m_triangles = triangles;
		m_taille_radius = taille_radius;
		m_trianglesTags = new bool[m_nTriangles];
		m_edges.reserve(3 * m_nTriangles);
		m_vertices.resize(m_nVertices);
		transferInd.resize(m_nVertices, -1);
		targetError = 0.1f;
		for (size_t v = 0; v < m_nVertices; ++v)
		{
			m_vertices[v].m_tag = true;
		}
		int tri[3];
		MDEdge edge;
		edge.m_tag = true;
		edge.m_onBoundary = true;
		int nEdges = 0;
		int idEdge;
		int nTris = static_cast<int>(m_nTriangles);
		for (int t = 0; t < nTris; ++t)
		{
			tri[0] = m_triangles[t].X();
			tri[1] = m_triangles[t].Y();
			tri[2] = m_triangles[t].Z();
			m_trianglesTags[t] = true;
			for (int k = 0; k < 3; ++k)
			{
				edge.m_v1 = tri[k];
				edge.m_v2 = tri[(k + 1) % 3];
				m_vertices[edge.m_v1].m_triangles.Insert(t);
				idEdge = GetEdge(edge.m_v1, edge.m_v2);
				if (idEdge == -1)
				{
					m_edges.push_back(edge);
					m_vertices[edge.m_v1].m_edges.Insert(nEdges);
					m_vertices[edge.m_v2].m_edges.Insert(nEdges);
					++nEdges;
				}
				else
				{
					m_edges[idEdge].m_onBoundary = false;
				}
			}
		}
		m_nEdges = static_cast<size_t>(nEdges);
		for (size_t v = 0; v < m_nVertices; ++v)
		{
			m_vertices[v].m_onBoundary = false;
			for (size_t itE = 0; itE < m_vertices[v].m_edges.Size(); ++itE)
			{
				idEdge = m_vertices[v].m_edges[itE];
				if (m_edges[idEdge].m_onBoundary)
				{
					m_vertices[v].m_onBoundary = true;
					break;
				}
			}
		}
	}
	void MeshDecimator::InitializeQEM()
	{
		Vec3<Float> coordMin = m_points[0];
		Vec3<Float> coordMax = m_points[0];
		Vec3<Float> coord;

		for (size_t p = 1; p < m_nPoints; ++p)
		{
			coord = m_points[p];
			if (coordMin.X() > coord.X()) coordMin.X() = coord.X();
			if (coordMin.Y() > coord.Y()) coordMin.Y() = coord.Y();
			if (coordMin.Z() > coord.Z()) coordMin.Z() = coord.Z();
			if (coordMax.X() < coord.X()) coordMax.X() = coord.X();
			if (coordMax.Y() < coord.Y()) coordMax.Y() = coord.Y();
			if (coordMax.Z() < coord.Z()) coordMax.Z() = coord.Z();
		}
		coordMax -= coordMin;
		m_diagBB = coordMax.GetNorm();

		int i, j, k;
		Vec3<Float> n;
		Float d = 0;
		Float area = 0;
		size_t tracker = m_nPoints / 20;
		std::cout << "\n" << "[" << std::flush;


		for (size_t v = 0; v < m_nPoints; ++v)
		{
			if ((v%tracker) == 0) {
				std::cout << "-" << std::flush;
			}
			memset(m_vertices[v].m_Q, 0, 10 * sizeof(Float));
			int idTriangle;
			for (size_t itT = 0; itT < m_vertices[v].m_triangles.Size(); ++itT)
			{
				idTriangle = m_vertices[v].m_triangles[itT];
				i = m_triangles[idTriangle].X();
				j = m_triangles[idTriangle].Y();
				k = m_triangles[idTriangle].Z();
				n = (m_points[j] - m_points[i]) ^ (m_points[k] - m_points[i]);
				area = n.GetNorm();
				n.Normalize();
				d = -(m_points[v] * n);
				m_vertices[v].m_Q[0] += area * (n.X() * n.X());
				m_vertices[v].m_Q[1] += area * (n.X() * n.Y());
				m_vertices[v].m_Q[2] += area * (n.X() * n.Z());
				m_vertices[v].m_Q[3] += area * (n.X() * d);
				m_vertices[v].m_Q[4] += area * (n.Y() * n.Y());
				m_vertices[v].m_Q[5] += area * (n.Y() * n.Z());
				m_vertices[v].m_Q[6] += area * (n.Y() * d);
				m_vertices[v].m_Q[7] += area * (n.Z() * n.Z());
				m_vertices[v].m_Q[8] += area * (n.Z() * d);
				m_vertices[v].m_Q[9] += area * (d     * d);
			}
		}
		std::cout << "]" << std::endl;
		Vec3<Float> u1, u2;
		const Float w = static_cast<Float>(1000);
		int t, v1, v2, v3;
		tracker = m_edges.size() / 20;
		std::cout << "\n" << "[" << std::flush;


		for (size_t e = 0; e < m_edges.size(); ++e)
		{
			if ((e%tracker) == 0) {
				std::cout << "-" << std::flush;
			}
			v1 = m_edges[e].m_v1;
			v2 = m_edges[e].m_v2;
			t = IsBoundaryEdge(v1, v2);
			if (t != -1)
			{
				if (m_triangles[t].X() != v1 && m_triangles[t].X() != v2) v3 = m_triangles[t].X();
				else if (m_triangles[t].Y() != v1 && m_triangles[t].Y() != v2) v3 = m_triangles[t].Y();
				else                                                           v3 = m_triangles[t].Z();
				u1 = m_points[v2] - m_points[v1];
				u2 = m_points[v3] - m_points[v1];
				area = w * (u1^u2).GetNorm();
				u1.Normalize();
				n = u2 - (u2 * u1) * u1;
				n.Normalize();

				d = -(m_points[v1] * n);
				m_vertices[v1].m_Q[0] += area * (n.X() * n.X());
				m_vertices[v1].m_Q[1] += area * (n.X() * n.Y());
				m_vertices[v1].m_Q[2] += area * (n.X() * n.Z());
				m_vertices[v1].m_Q[3] += area * (n.X() * d);
				m_vertices[v1].m_Q[4] += area * (n.Y() * n.Y());
				m_vertices[v1].m_Q[5] += area * (n.Y() * n.Z());
				m_vertices[v1].m_Q[6] += area * (n.Y() * d);
				m_vertices[v1].m_Q[7] += area * (n.Z() * n.Z());
				m_vertices[v1].m_Q[8] += area * (n.Z() * d);
				m_vertices[v1].m_Q[9] += area * (d * d);

				d = -(m_points[v2] * n);
				m_vertices[v2].m_Q[0] += area * (n.X() * n.X());
				m_vertices[v2].m_Q[1] += area * (n.X() * n.Y());
				m_vertices[v2].m_Q[2] += area * (n.X() * n.Z());
				m_vertices[v2].m_Q[3] += area * (n.X() * d);
				m_vertices[v2].m_Q[4] += area * (n.Y() * n.Y());
				m_vertices[v2].m_Q[5] += area * (n.Y() * n.Z());
				m_vertices[v2].m_Q[6] += area * (n.Y() * d);
				m_vertices[v2].m_Q[7] += area * (n.Z() * n.Z());
				m_vertices[v2].m_Q[8] += area * (n.Z() * d);
				m_vertices[v2].m_Q[9] += area * (d * d);
			}
		}
		std::cout << "]" << std::endl;
	}
	bool MeshDecimator::EdgeCollapse(double & qem)
	{
		MDEdgePriorityQueue currentEdge;
		int v1, v2;
		bool done = false;
		do
		{
			done = false;
			if (m_pqueue.size() == 0)
			{
				done = true;
				break;
			}
			else
			{
				currentEdge = m_pqueue.top();
				m_pqueue.pop();
			}
		} while ((!m_edges[currentEdge.m_name].m_tag) || (m_edges[currentEdge.m_name].m_qem != currentEdge.m_qem));

		if (done) return false;
		v1 = m_edges[currentEdge.m_name].m_v1;
		v2 = m_edges[currentEdge.m_name].m_v2;

		qem = currentEdge.m_qem;
		//ADD condition
		if (IsDistanceBelowFLoatRadius(m_points[v1], m_points[v2])) {
			EdgeCollapse(v1, v2);
			m_points[v1] = m_edges[currentEdge.m_name].m_pos;
			for (int k = 0; k < 10; k++) m_vertices[v1].m_Q[k] += m_vertices[v2].m_Q[k];
			
			// Update priority queue
			int idEdge;
			int a, b;
			SArray<int, SARRAY_DEFAULT_MIN_SIZE> incidentVertices;
			for (size_t itE = 0; itE < m_vertices[v1].m_edges.Size(); ++itE)
			{
				idEdge = m_vertices[v1].m_edges[itE];
				a = m_edges[idEdge].m_v1;
				b = m_edges[idEdge].m_v2;
				incidentVertices.PushBack((a != v1) ? a : b);
				MDEdgePriorityQueue pqEdge;
				pqEdge.m_qem = m_edges[idEdge].m_qem = ComputeEdgeCost(a, b, m_edges[idEdge].m_pos);
				pqEdge.m_name = idEdge;
				m_pqueue.push(pqEdge);
			}
			int idVertex;
			for (size_t itV = 0; itV < incidentVertices.Size(); ++itV)
			{
				idVertex = incidentVertices[itV];
				for (size_t itE = 0; itE < m_vertices[idVertex].m_edges.Size(); ++itE)
				{
					idEdge = m_vertices[idVertex].m_edges[itE];
					a = m_edges[idEdge].m_v1;
					b = m_edges[idEdge].m_v2;
					if (a != v1 && b != v1)
					{
						MDEdgePriorityQueue pqEdge;
						pqEdge.m_qem = m_edges[idEdge].m_qem = ComputeEdgeCost(a, b, m_edges[idEdge].m_pos);
						pqEdge.m_name = idEdge;
						m_pqueue.push(pqEdge);
					}
				}
			}
		}
		else {
			if (!m_pqueue.empty()) {
				m_pqueue.pop();
			}
		}
		return true;
	}
	bool MeshDecimator::Decimate()
	{
		double qem = 0.0;
		double progressOld = -1.0;
		double progress = 0.0;
		char msg[1024];
		double ptgStep = 1.0;

		std::cout <<"\n"<< "Initialize QEM" << std::flush;
		InitializeQEM();
		std::cout << "Initialize priority queue" << std::flush;
		InitializePriorityQueue();
		std::cout << "Simplification" << std::flush;
		double invDiag = 1.0 / m_diagBB;
		std::cout << "\n" << "["<<std::flush;
		size_t tracker = m_pqueue.size() / 20;
		int nvr = 1;
		while ((m_pqueue.size() > 0) &&
			(m_nEdges > 0)&&(qem < targetError))
		{
			if ((tracker*20-nvr*tracker) <m_pqueue.size()) {
				nvr++;
				std::cout << '-' << std::flush;
			}
			
			
			if (!EdgeCollapse(qem)) break;
			if (qem < 0.0) qem = 0.0;
			else           qem = sqrt(qem) * invDiag;
		}
		std::cout << "]" << std::endl;
		
		return true;
	}
	std::pair< std::vector<glm::vec3>, std::vector<glm::uvec3>>  QEMDecimateForASpecifiedListOfEdges(std::vector<glm::vec3> VertexPositions, std::vector<glm::uvec3> trianglesIndices, float taille_radius, std::shared_ptr<Mesh> meshPtr) {
		std::vector<float> tranfer_Ind;
		std::cout << "Data Transfer Loading" << std::flush;
		MeshDecimation::Vec3<MeshDecimation::Float> * points = new MeshDecimation::Vec3<MeshDecimation::Float>[VertexPositions.size()];
		for (size_t i = 0; i < VertexPositions.size(); i++)
		{
			points[i] = MeshDecimation::Vec3<MeshDecimation::Float>(VertexPositions[i][0], VertexPositions[i][1], VertexPositions[i][2]);
		}

		MeshDecimation::Vec3<int> * triangles = new MeshDecimation::Vec3<int>[trianglesIndices.size()];
		for (size_t i = 0; i < trianglesIndices.size(); i++)
		{

			triangles[i] = MeshDecimation::Vec3<int>(trianglesIndices[i][0], trianglesIndices[i][1], trianglesIndices[i][2]);
		}
		int nPoints = VertexPositions.size();
		int nTriangles = trianglesIndices.size();

		std::cout << "---------->Data Tranfert Loaded" << std::endl;
		// decimate mesh

		MeshDecimation::MeshDecimator myMDecimator;
		std::cout << "Initialisation Loading" << std::flush;
		myMDecimator.Initialize(nPoints, taille_radius,nTriangles, points, triangles);
		std::cout << "---------->Initialisation Loaded" << std::endl;
		std::cout << "Decimation Loading" << std::flush;
		myMDecimator.Decimate();
		std::cout << "---------->Decimation Loaded" << std::endl;
		points[VertexPositions.size()].~Vec3();
		triangles[trianglesIndices.size()].~Vec3();

		// allocate memory for decimated mesh
		size_t nVerticesDecimatedMesh = myMDecimator.GetNVertices();
		size_t nTrianglesDecimatedMesh = myMDecimator.GetNTriangles();
		tranfer_Ind = myMDecimator.GetTransferVertex();
		meshPtr->vertex_transferGroupe = tranfer_Ind;
		MeshDecimation::Vec3 < MeshDecimation::Float > * decimatedMeshPoints = new MeshDecimation::Vec3<MeshDecimation::Float>[nVerticesDecimatedMesh];
		MeshDecimation::Vec3<int> * decimatedMeshTriangles = new MeshDecimation::Vec3<int>[nTrianglesDecimatedMesh];
		std::cout << "Get Mesh Loading" << std::flush;
		// retrieve decimated mesh
		myMDecimator.GetMeshData(decimatedMeshPoints, decimatedMeshTriangles);
		std::cout << "---------->Get Mesh Loaded" << std::endl;
		VertexPositions.clear();
		std::cout << "Tranfert Data Loading" << std::flush;
		for (size_t i = 0; i < nVerticesDecimatedMesh; i++)
		{
			VertexPositions.push_back(glm::vec3(decimatedMeshPoints[i][0], decimatedMeshPoints[i][1], decimatedMeshPoints[i][2]));
		}

		trianglesIndices.clear();
	
		std::cout << "---------->Transfer Data Loaded" << std::endl;
		std::cout << "Free data Loading" << std::flush;
		myMDecimator.~MeshDecimator();
		std::cout << "1" << std::flush;
		decimatedMeshPoints[nVerticesDecimatedMesh].~Vec3();
		std::cout << "2" << std::flush;
		decimatedMeshTriangles[nTrianglesDecimatedMesh].~Vec3();

		std::cout << "---------->Free data Loaded" << std::endl;
		std::cout << "Pair Creation Loading" << std::flush;
		std::pair< std::vector<glm::vec3>, std::vector<glm::uvec3>> result(VertexPositions, trianglesIndices);
		VertexPositions.~vector();
		trianglesIndices.~vector();
		std::cout << "---------->Pair Creation Loaded" << std::endl;
		return result;
	}
	bool MeshDecimator::IsDistanceBelowFLoatRadius(Vec3<Float> v1, Vec3<Float> v2) {
		glm::vec3 v1_, v2_;
		v1_ = glm::vec3(v1[0], v1[1], v1[2]);
		v2_ = glm::vec3(v2[0], v2[1], v2[2]);
		if (glm::distance(v1_, v2_) < m_taille_radius) {
			return true;
		}
		else {
			return false;
		}

	}
}