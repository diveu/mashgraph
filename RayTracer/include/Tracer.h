#pragma once

#include "glm.hpp"
#include "EasyBMP.h"

#include "Types.h"
#include "Scene.h"

#include "string"


class CTracer
{
public:
    SRay MakeRay(glm::uvec2 pixelPos, glm::vec2 pl);   // Create ray for specified pixel
    glm::vec3 TraceRay(SRay ray); // Trace ray, compute its color
    void RenderImage(int xRes, int yRes, int AA);
    void SaveImageToFile(std::string fileName);
    void lighting(SPhotonMap &m_points, int n, float s);
    float FresnelReflectance(float ci, float n);
public:
    SCamera m_camera;
    SPhotonMap m_photonMap;
    CScene* m_pScene;
    SLight m_light;
    SLens m_lens;
    glm::vec3 right_width;
    glm::vec3 up_height;
    SIntrsectionObj CheckRayCrossMesh(SMesh& mesh, SRay& ray);
    SIntrsectionObj CheckRayCrossSphere(SSphere& sphere, SRay& ray);
    glm::vec3 ProcessRay(SRay& ray, SIntrsectionObj& IntersectedObj);
};
