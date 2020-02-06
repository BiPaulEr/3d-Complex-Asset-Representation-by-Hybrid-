
/* Copyright (c) 2011 Khaled Mamou 
   Code has been modify by Paul-Ernest Martin to be adapted to the specifity of the use by this project 
   Attributes, functions and conditions have been modified or added
*/
#pragma once
#ifndef MD_MESH_DECEMATOR_H
#define MD_MESH_DECEMATOR_H
#include <queue>
#include <set>
#include <vector>
#include <limits>
#include "mdVector.h"
#include "mdSArray.h"
#include "Mesh.h"

namespace MeshDecimation
{
	//    typedef double Float;
	typedef float Float;
	struct MDVertex
	{
		SArray<int, SARRAY_DEFAULT_MIN_SIZE>    m_edges;
		SArray<int, SARRAY_DEFAULT_MIN_SIZE>    m_triangles;
		Float                                   m_Q[10];
		// 0 1 2 3
		//   4 5 6
		//     7 8
		//       9
		bool                                    m_tag;
		bool                                    m_onBoundary;
	};

	struct MDEdge
	{
		int                                     m_v1;
		int                                     m_v2;
		double                                  m_qem;
		Vec3<Float>                             m_pos;
		bool                                    m_onBoundary;
		bool                                    m_tag;
	};
	struct MDEdgePriorityQueue
	{
		int                                     m_name;
		double                                  m_qem;
		inline    friend bool                   operator<(const MDEdgePriorityQueue & lhs, const MDEdgePriorityQueue & rhs) { return (lhs.m_qem > rhs.m_qem); }
		inline    friend bool                   operator>(const MDEdgePriorityQueue & lhs, const MDEdgePriorityQueue & rhs) { return (lhs.m_qem < rhs.m_qem); }
	};
	typedef void(*CallBackFunction)(const char * msg);


	// Function added to enable the transfer of data from the mesh analysis to the decimation QEM librairy and freeing the data linked
	std::pair< std::vector<glm::vec3>, std::vector<glm::uvec3>>  QEMDecimateForASpecifiedListOfEdges(std::vector<glm::vec3> VertexPositions, std::vector<glm::uvec3> trianglesIndices, float taille_radius,  std::shared_ptr<Mesh> meshPtr);

	class MeshDecimator
	{
	public:
		//! Sets the call-back function
		//! @param callBack pointer to the call-back function
		void                                    SetCallBack(CallBackFunction  callBack) { m_callBack = callBack; }
		//! Gives the call-back function
		//! @return pointer to the call-back function
		const CallBackFunction                  GetCallBack() const { return m_callBack; }

		inline void                             SetEColManifoldConstraint(bool ecolManifoldConstraint) { m_ecolManifoldConstraint = ecolManifoldConstraint; }
		inline size_t                           GetNVertices()const { return m_nVertices; };
		inline size_t                           GetNTriangles() const { return m_nTriangles; };
		inline size_t                           GetNEdges() const { return m_nEdges; };
		inline std::vector<float>               GetTransferVertex() const { return transferInd; };
		void                                    GetMeshData(Vec3<Float> * points, Vec3<int> * triangles) const;
		void                                    ReleaseMemory();


		//edges collapse has been modify to enable the attributes added
		void                                    Initialize(size_t nVertices,  Float taille_radius,size_t nTriangles,Vec3<Float> *  points,Vec3<int> * triangles);
		//decimatee has been modify to add a condition
		bool                                    Decimate();
		//add to enable the condition of decimation
		bool IsDistanceBelowFLoatRadius(Vec3<Float> v1, Vec3<Float> v2);

		// Modify to adapt to added attributes and tranfer of data
		MeshDecimator(void);
		~MeshDecimator(void);
	private:
		void                                    EdgeCollapse(int v1, int v2);
		int                                     GetTriangle(int v1, int v2, int v3) const;
		int                                     GetEdge(int v1, int v2) const;
		int                                     IsBoundaryEdge(int v1, int v2) const;
		bool                                    IsBoundaryVertex(int v) const;
		void                                    InitializePriorityQueue();

		void                                    InitializeQEM();
		bool                                    ManifoldConstraint(int v1, int v2) const;
		double                                  ComputeEdgeCost(int v1, int v2, Vec3<Float> & pos) const;
		//edges collapse has been modify
		bool                                    EdgeCollapse(double & error);
		
	private:

		// add to enable transfer of index of each vertex decimed
		std::vector<float>						transferInd;
		// add to enable a condition in the size of the edges
		Float								    m_taille_radius;
		Vec3<int> *                             m_triangles;
		Vec3<Float> *                           m_points;
		size_t                                  m_nPoints;
		size_t                                  m_nInitialTriangles;
		size_t                                  m_nVertices;
		size_t                                  m_nTriangles;
		size_t                                  m_nEdges;
		double                                  m_diagBB;
		std::vector<MDVertex>                   m_vertices;
		std::vector<MDEdge>                     m_edges;
		Float									targetError;
		std::priority_queue<
			MDEdgePriorityQueue,
			std::vector<MDEdgePriorityQueue>,
			std::less<MDEdgePriorityQueue> >   m_pqueue;
		CallBackFunction                        m_callBack;                    //>! call-back function
		bool *                                  m_trianglesTags;
		bool                                    m_ecolManifoldConstraint;
	};
}
#endif