/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2025 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#include <vector>
#include <unordered_map>
#include <array>
#include <FECore/vec3d.h>
#include "feamr_api.h"

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ (h2 << 1);
    }
};

enum meshSizes { SMALL = 4, FULL = 6};

class FEAMR_API SpherePointsGenerator
{
public:
    static size_t GetNumNodes(int n);
    static size_t GetNumFaces(int n);
    static std::vector<vec3d>& GetNodes(int n);
    static std::vector<std::array<int, 3>>& GetFaces(int n);

private:
	SpherePointsGenerator() {}
    static void Instantiate();

    void GeneratePoints(int n);

    vec3d midpoint(const vec3d& a, const vec3d& b);

    // Create modified icosahedron which has a belt around the equator
    void create_icosahedron(std::vector<vec3d>& vertices, std::vector<std::array<int, 3>>& faces);

    void rotate(std::vector<vec3d>& vertices, const vec3d& axis, double angle_rad);

    // Aligns the equator of the sphere with the x-y plane
    void align_equator(std::vector<vec3d>& vertices);

    // Get or create midpoint between two vertices
    int get_midpoint_index(int i1, int i2, std::vector<vec3d>& vertices, 
        std::unordered_map<std::pair<int, int>, int, pair_hash>& midpoint_cache);

    // Subdivide icosahedral mesh
    void subdivide(std::vector<vec3d>& vertices, std::vector<std::array<int, 3>>& faces, int depth);

    
private:
    static SpherePointsGenerator* m_instance;

    std::unordered_map<int, std::vector<vec3d>> m_nodes;
    std::unordered_map<int, std::vector<std::array<int, 3>>> m_faces;
};
