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

#include <algorithm>
#include <cassert>
#include "SpherePointsGenerator.h"

SpherePointsGenerator* SpherePointsGenerator::m_instance = nullptr;

void SpherePointsGenerator::Instantiate()
{
    if (m_instance == nullptr) {
        m_instance = new SpherePointsGenerator();
    }
}

size_t SpherePointsGenerator::GetNumNodes(int n)
{
    Instantiate();

    m_instance->GeneratePoints(n);

    return m_instance->m_nodes[n].size();
}

size_t SpherePointsGenerator::GetNumFaces(int n)
{
    Instantiate();

    m_instance->GeneratePoints(n);

    return m_instance->m_faces[n].size();
}

std::vector<vec3d>& SpherePointsGenerator::GetNodes(int n)
{
    Instantiate();

    m_instance->GeneratePoints(n);

    return m_instance->m_nodes[n];
}

std::vector<std::array<int, 3>>& SpherePointsGenerator::GetFaces(int n)
{
    Instantiate();

    m_instance->GeneratePoints(n);

    return m_instance->m_faces[n];
}

void SpherePointsGenerator::GeneratePoints(int n)
{
    if(m_nodes.find(n) == m_nodes.end() || m_faces.find(n) == m_faces.end())
    {
        m_nodes[n] = std::vector<vec3d>();
        m_faces[n] = std::vector<std::array<int, 3>>();

        std::vector<vec3d>& vertices = m_nodes[n];
        std::vector<std::array<int, 3>>& faces = m_faces[n];

        // Create icosahedron
        create_icosahedron(vertices, faces);
    
        // Align equator with x-y plane
        align_equator(vertices);
        
        // Subdivide
        subdivide(vertices, faces, n);
    }
}

vec3d SpherePointsGenerator::midpoint(const vec3d& a, const vec3d& b)
{
    vec3d out = (a + b).normalized();

    // Ensure that points on the equator are exactly on the x-y plane
    if(abs(out.z) < 1e-5)
    {
        out.z = 0;
    }

    return out;
}

void SpherePointsGenerator::create_icosahedron(std::vector<vec3d>& vertices, std::vector<std::array<int, 3>>& faces)
{
    double phi = (1.0 + sqrt(5.0)) / 2.0;
    
    // Match Python vertices exactly
    std::vector<vec3d> verts =
    {
        {0, 1, phi}, {0, -1, phi}, {0, 1, -phi}, {0, -1, -phi},
        {1, phi, 0}, {-1, phi, 0}, {1, -phi, 0}, {-1, -phi, 0},
        {phi, 0, 1}, {-phi, 0, 1}, {phi, 0, -1}, {-phi, 0, -1}
    };
    
    // Normalize vertices
    for (auto& v : verts)
    {
        v.Normalize();
    }
    
    vertices = verts;
    
    // Use same face connectivity as in Python
    faces =
    {
        {0, 1, 8}, {0, 8, 4}, {0, 4, 5}, {0, 5, 9}, {0, 9, 1},
        {1, 9, 7}, {1, 7, 6}, {1, 6, 8}, {8, 6, 10}, {8, 10, 4},
        {4, 10, 2}, {4, 2, 5}, {5, 2, 11}, {5, 11, 9}, {9, 11, 7},
        {3, 2, 10}, {3, 10, 6}, {3, 6, 7}, {3, 7, 11}, {3, 11, 2}
    };
}

// Rotation function to match Python's implementation
void SpherePointsGenerator::rotate(std::vector<vec3d>& vertices, const vec3d& axis, double angle_rad)
{
    // Normalize axis
    vec3d norm_axis = axis.normalized();
    
    // Create rotation matrix using Rodrigues' formula (same as Python)
    double cos_a = cosf(angle_rad);
    double sin_a = sinf(angle_rad);
    double one_minus_cos = 1.0f - cos_a;
    
    // K matrix (skew-symmetric cross-product matrix)
    double K[3][3] =
    {
        {0, -norm_axis.z, norm_axis.y},
        {norm_axis.z, 0, -norm_axis.x},
        {-norm_axis.y, norm_axis.x, 0}
    };
    
    // K^2 matrix
    double K2[3][3] =
    {
        {-norm_axis.y*norm_axis.y - norm_axis.z*norm_axis.z, norm_axis.x*norm_axis.y, norm_axis.x*norm_axis.z},
        {norm_axis.x*norm_axis.y, -norm_axis.x*norm_axis.x - norm_axis.z*norm_axis.z, norm_axis.y*norm_axis.z},
        {norm_axis.x*norm_axis.z, norm_axis.y*norm_axis.z, -norm_axis.x*norm_axis.x - norm_axis.y*norm_axis.y}
    };
    
    // R = I + sin(θ)K + (1-cos(θ))K²
    double R[3][3] =
    {
        {1.0f + sin_a*K[0][0] + one_minus_cos*K2[0][0], sin_a*K[0][1] + one_minus_cos*K2[0][1], sin_a*K[0][2] + one_minus_cos*K2[0][2]},
        {sin_a*K[1][0] + one_minus_cos*K2[1][0], 1.0f + sin_a*K[1][1] + one_minus_cos*K2[1][1], sin_a*K[1][2] + one_minus_cos*K2[1][2]},
        {sin_a*K[2][0] + one_minus_cos*K2[2][0], sin_a*K[2][1] + one_minus_cos*K2[2][1], 1.0f + sin_a*K[2][2] + one_minus_cos*K2[2][2]}
    };
    
    // Apply rotation to each vertex
    for (auto& v : vertices)
    {
        double x = v.x, y = v.y, z = v.z;
        v.x = R[0][0]*x + R[0][1]*y + R[0][2]*z;
        v.y = R[1][0]*x + R[1][1]*y + R[1][2]*z;
        v.z = R[2][0]*x + R[2][1]*y + R[2][2]*z;
    }
}

// Match Python's align_equator implementation
void SpherePointsGenerator::align_equator(std::vector<vec3d>& vertices)
{
    // Step 1: Detect "belt" axis
    double abs_mean_x = 0, abs_mean_y = 0, abs_mean_z = 0;
    for (const auto& v : vertices)
    {
        abs_mean_x += abs(v.x);
        abs_mean_y += abs(v.y);
        abs_mean_z += abs(v.z);
    }
    
    int n = (int)vertices.size();
    abs_mean_x /= n;
    abs_mean_y /= n;
    abs_mean_z /= n;
    
    int belt_axis = (abs_mean_x < abs_mean_y && abs_mean_x < abs_mean_z) ? 0 : 
                   (abs_mean_y < abs_mean_z) ? 1 : 2;
    
    // Step 2: Find vertices in the belt
    std::vector<vec3d> band_vertices;
    std::vector<double> abs_coords;
    
    for (const auto& v : vertices)
    {
        double coord = (belt_axis == 0) ? v.x : (belt_axis == 1) ? v.y : v.z;
        abs_coords.push_back(abs(coord));
    }
    
    // Calculate median of absolute coordinates (approximate)
    std::vector<double> sorted_abs_coords = abs_coords;
    std::sort(sorted_abs_coords.begin(), sorted_abs_coords.end());
    double belt_band = sorted_abs_coords[sorted_abs_coords.size() / 2];
    
    for (size_t i = 0; i < vertices.size(); i++)
    {
        double coord = (belt_axis == 0) ? vertices[i].x : (belt_axis == 1) ? vertices[i].y : vertices[i].z;
        if (abs(abs(coord) - belt_band) < 1e-3f)
        {
            band_vertices.push_back(vertices[i]);
        }
    }
    
    if (band_vertices.empty())
    {
        
        assert(false); // No vertices in belt
    }
    
    // Pick first band vertex and normalize
    vec3d target_vertex = band_vertices[0].normalized();
    
    // Step 3: Rotate this vertex to Z axis
    vec3d z_axis(0, 0, 1);
    vec3d axis = vec3d(
        target_vertex.y * z_axis.z - target_vertex.z * z_axis.y,
        target_vertex.z * z_axis.x - target_vertex.x * z_axis.z,
        target_vertex.x * z_axis.y - target_vertex.y * z_axis.x
    );
    
    axis.Normalize();
    double dot = target_vertex.x * z_axis.x + target_vertex.y * z_axis.y + target_vertex.z * z_axis.z;
    double angle = acos(std::clamp(dot, -1.0, 1.0));
    
    // Rotate all vertices
    rotate(vertices, axis, angle);
    
    // Ensure that points on the equator are exactly on the x-y plane
    for(int index = 0; index < vertices.size(); index++)
    {
        if(abs(vertices[index].z) < 1e-5)
        {
            vertices[index].z = 0;
        }
    }
}

int SpherePointsGenerator::get_midpoint_index(int i1, int i2, std::vector<vec3d>& vertices, 
                       std::unordered_map<std::pair<int, int>, int, pair_hash>& midpoint_cache)
{
    std::pair<int, int> key = (i1 < i2) ? std::make_pair(i1, i2) : std::make_pair(i2, i1);
    
    if (midpoint_cache.count(key))
    {
        return midpoint_cache[key];
    }
    
    // Calculate midpoint 
    vec3d mid = midpoint(vertices[i1], vertices[i2]);
    vertices.push_back(mid);
    int idx = (int)vertices.size() - 1;
    midpoint_cache[key] = idx;
    return idx;
}

void SpherePointsGenerator::subdivide(std::vector<vec3d>& vertices, std::vector<std::array<int, 3>>& faces, int depth)
{
    for (int i = 0; i < depth; ++i)
    {
        std::unordered_map<std::pair<int, int>, int, pair_hash> midpoint_cache;
        std::vector<std::array<int, 3>> new_faces;
        
        for (const auto& tri : faces)
        {
            int v1 = tri[0], v2 = tri[1], v3 = tri[2];
            
            // Get midpoints
            int a = get_midpoint_index(v1, v2, vertices, midpoint_cache);
            int b = get_midpoint_index(v2, v3, vertices, midpoint_cache);
            int c = get_midpoint_index(v3, v1, vertices, midpoint_cache);
            
            // Create new faces 
            new_faces.push_back({v1, a, c});
            new_faces.push_back({v2, b, a});
            new_faces.push_back({v3, c, b});
            new_faces.push_back({a, b, c});
        }
        
        faces = new_faces;
    }
}