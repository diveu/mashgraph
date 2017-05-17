#include "Scene.h"
#include "glm.hpp"
#include "vector"

using namespace glm;

CScene::CScene(glm::vec3 LightPos, glm::vec3 LightDir)
{
    const int sceneSize = 40;
    const int box_size = 10;
    const float sphere_rad = 17;
    
    sphere_.centre = vec3(sphere_rad, sceneSize - sphere_rad, sphere_rad / 2.0);
    sphere_.r = vec3(sphere_rad, sphere_rad, sphere_rad) / sqrt(3.0f);
    float a = M_PI / 8.0f;
    glm::mat3x3 rot = mat3x3 (cos(a), -sin(a), 0,
                              sin(a),  cos(a), 0,
                                  0 ,      0 , 1 );
    
    light_.m_vertices.push_back(LightPos);
    light_.m_vertices.push_back(LightPos + 3 * rot * LightDir);
    light_.m_vertices.push_back(LightPos + 6 * rot * LightDir);
    light_.m_vertices.push_back(LightPos - 3 * rot * LightDir - LightDir);

    a = -a;
    rot = mat3x3 (cos(a), -sin(a), 0,
                  sin(a),  cos(a), 0,
                      0 ,      0 , 1 );
    light_.m_vertices.push_back(LightPos + 3 * rot * LightDir);
    light_.m_vertices.push_back(LightPos + 6 * rot * LightDir);
    light_.m_vertices.push_back(LightPos - 3 * rot * LightDir - LightDir);
    
    for (uint j = 0; j < light_.m_vertices.size(); j++)
        light_.m_vertices[j] += vec3 (0, 0, -1);
    
    light_.m_triangles.push_back(uvec3(0, 2, 5));
    light_.m_triangles.push_back(uvec3(1, 6, 3));
    light_.m_triangles.push_back(uvec3(1, 4, 3));
    
    box_right_wall_.m_vertices.push_back(vec3(box_size, -4 * box_size, -box_size));
    box_right_wall_.m_vertices.push_back(vec3(box_size, box_size, -box_size));
    box_right_wall_.m_vertices.push_back(vec3(box_size, box_size,  box_size));
    box_right_wall_.m_vertices.push_back(vec3(box_size, -4 * box_size, box_size));
    
    box_left_wall_.m_vertices.push_back(vec3(-box_size, -4 * box_size, -box_size));
    box_left_wall_.m_vertices.push_back(vec3(-box_size, box_size, -box_size));
    box_left_wall_.m_vertices.push_back(vec3(-box_size, box_size,  box_size));
    box_left_wall_.m_vertices.push_back(vec3(-box_size, -4 * box_size, box_size));
    
    box_front_wall_.m_vertices.push_back(vec3(-box_size, -4 * box_size, box_size));
    box_front_wall_.m_vertices.push_back(vec3(-box_size, box_size, box_size));
    box_front_wall_.m_vertices.push_back(vec3(box_size, box_size, box_size));
    box_front_wall_.m_vertices.push_back(vec3(box_size, -4 * box_size, box_size));
    
    box_ceiling_.m_vertices.push_back(vec3(-box_size, -4 * box_size, -box_size));
    box_ceiling_.m_vertices.push_back(vec3(-box_size, -4 * box_size, box_size));
    box_ceiling_.m_vertices.push_back(vec3(box_size, -4 * box_size, box_size));
    box_ceiling_.m_vertices.push_back(vec3(box_size, -4 * box_size, -box_size));
    
    const vec3 box_shift = vec3(-2 * box_size, sceneSize - box_size, - 3 * box_size / 1.5);
    
    for (uint i = 0; i < box_right_wall_.m_vertices.size(); i++)
        box_right_wall_.m_vertices[i] += box_shift;
    for (uint i = 0; i < box_left_wall_.m_vertices.size(); i++)
        box_left_wall_.m_vertices[i] += box_shift;
    for (uint i = 0; i < box_front_wall_.m_vertices.size(); i++)
        box_front_wall_.m_vertices[i] += box_shift;
    for (uint i = 0; i < box_ceiling_.m_vertices.size(); i++)
        box_ceiling_.m_vertices[i] += box_shift;
    
    box_right_wall_.m_triangles.push_back(uvec3(0, 1, 2));
    box_right_wall_.m_triangles.push_back(uvec3(0, 2, 3));
    box_left_wall_.m_triangles.push_back(uvec3(0, 1, 2));
    box_left_wall_.m_triangles.push_back(uvec3(0, 2, 3));
    box_front_wall_.m_triangles.push_back(uvec3(0, 1, 2));
    box_front_wall_.m_triangles.push_back(uvec3(0, 2, 3));
    box_ceiling_.m_triangles.push_back(uvec3(0, 1, 2));
    box_ceiling_.m_triangles.push_back(uvec3(0, 2, 3));
    
    
    ceiling_.m_vertices.push_back(vec3(sceneSize, -sceneSize, sceneSize));
    ceiling_.m_vertices.push_back(vec3(sceneSize, -sceneSize, -sceneSize));
    ceiling_.m_vertices.push_back(vec3(-sceneSize, -sceneSize, -sceneSize));
    ceiling_.m_vertices.push_back(vec3(-sceneSize, -sceneSize, sceneSize));
    ceiling_.m_triangles.push_back(uvec3(0, 2, 1));
    ceiling_.m_triangles.push_back(uvec3(3, 2, 0));
    
    floor_.m_vertices.push_back(vec3(sceneSize, sceneSize, sceneSize));
    floor_.m_vertices.push_back(vec3(sceneSize, sceneSize, -sceneSize));
    floor_.m_vertices.push_back(vec3(-sceneSize, sceneSize, -sceneSize));
    floor_.m_vertices.push_back(vec3(-sceneSize, sceneSize, sceneSize));
    floor_.m_triangles.push_back(uvec3(0, 2, 1));
    floor_.m_triangles.push_back(uvec3(3, 2, 0));
    
    left_wall_.m_vertices.push_back(vec3(-sceneSize, -sceneSize, sceneSize));
    left_wall_.m_vertices.push_back(vec3(-sceneSize, -sceneSize, -sceneSize));
    left_wall_.m_vertices.push_back(vec3(-sceneSize, sceneSize, -sceneSize));
    left_wall_.m_vertices.push_back(vec3(-sceneSize, sceneSize, sceneSize));
    left_wall_.m_triangles.push_back(uvec3(0, 2, 1));
    left_wall_.m_triangles.push_back(uvec3(3, 2, 0));
    
    right_wall_.m_vertices.push_back(vec3(sceneSize, -sceneSize, sceneSize));
    right_wall_.m_vertices.push_back(vec3(sceneSize, -sceneSize, -sceneSize));
    right_wall_.m_vertices.push_back(vec3(sceneSize, sceneSize, -sceneSize));
    right_wall_.m_vertices.push_back(vec3(sceneSize, sceneSize, sceneSize));
    right_wall_.m_triangles.push_back(uvec3(0, 2, 1));
    right_wall_.m_triangles.push_back(uvec3(3, 2, 0));
    
    back_wall_.m_vertices.push_back(vec3(-sceneSize, -sceneSize, -sceneSize));
    back_wall_.m_vertices.push_back(vec3(sceneSize, -sceneSize, -sceneSize));
    back_wall_.m_vertices.push_back(vec3(sceneSize, sceneSize, -sceneSize));
    back_wall_.m_vertices.push_back(vec3(-sceneSize, sceneSize, -sceneSize));
    back_wall_.m_triangles.push_back(uvec3(0, 2, 1));
    back_wall_.m_triangles.push_back(uvec3(3, 2, 0));
}
