#pragma once

#include "glm.hpp"
#include "vector"


enum IntersectionType
{
    NONE = 0,
    UNDEFINED,
    FLOOR,
    CEILING,
    LEFT_WALL,
    RIGHT_WALL,
    BACK_WALL,
    BOX_CEILING,
    BOX_LEFT_WALL,
    BOX_RIGHT_WALL,
    BOX_FRONT_WALL,
    SPHERE,
    LIGHT
};

struct SRay
{
    glm::vec3 m_start;
    glm::vec3 m_dir;
};

struct SLight {
    glm::vec3 m_pos;
    glm::vec3 m_dir;
};

struct SCamera
{
    glm::vec3 m_pos;          // Camera position and orientation
    glm::vec3 m_forward;      // Orthonormal basis
    glm::vec3 m_right;
    glm::vec3 m_up;

    glm::vec2 m_viewAngle;    // View angles, rad
    glm::uvec2 m_resolution;  // Image resolution: w, h

    std::vector<glm::vec3> m_pixels;  // Pixel array
};

struct SMesh
{
    std::vector<glm::vec3> m_vertices;  // vertex positions
    std::vector<glm::uvec3> m_triangles;  // vetrex indices
};

struct SSphere {
    glm::vec3 centre; //sphere centre
    glm::vec3 r; //sphere radius
};

struct SPhotonMap {
    std::vector<float> m_pixels;
    int ResX;
    int ResY;
    int Rad;
};

struct SLens {
    float R1;
    float R2;
    float n;
};

struct SIntrsectionObj
    {
        float distance;
        IntersectionType mesh_type;
        unsigned int triangle_num;
        
        SIntrsectionObj() : distance(0), mesh_type(NONE), triangle_num(0) {}
    };
