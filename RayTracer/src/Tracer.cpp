#include "EasyBMP.h"

#include "Tracer.h"

using namespace glm;

SRay CTracer::MakeRay(uvec2 pixelPos, glm::vec2 pl)
{
    SRay Ray;
    Ray.m_start = m_camera.m_pos;
    float x = (pixelPos.x + 0.5f +  pl.x) / m_camera.m_resolution.x - 0.5f;
    float y = (pixelPos.y + 0.5f + pl.y) / m_camera.m_resolution.y - 0.5f;
    vec3 ray_dir = m_camera.m_forward + right_width * (2 * x) + up_height * (2 * y);
    Ray.m_dir = glm::normalize(ray_dir);
    
    return Ray;
}

glm::vec3 CTracer::TraceRay(SRay ray)
{
    SIntrsectionObj IntersectedObj, tempIntersectedObj;

    //checking every object
    tempIntersectedObj = CheckRayCrossMesh(m_pScene->floor_, ray);
    if (tempIntersectedObj.mesh_type == UNDEFINED)
        if (tempIntersectedObj.distance < IntersectedObj.distance || IntersectedObj.mesh_type == NONE)
        {
            IntersectedObj = tempIntersectedObj;
            IntersectedObj.mesh_type = FLOOR;
        }
    tempIntersectedObj = CheckRayCrossMesh(m_pScene->back_wall_, ray); 
    if (tempIntersectedObj.mesh_type == UNDEFINED)
        if (tempIntersectedObj.distance < IntersectedObj.distance || IntersectedObj.mesh_type == NONE)
        {
            IntersectedObj = tempIntersectedObj;
            IntersectedObj.mesh_type = BACK_WALL;
        }
    tempIntersectedObj = CheckRayCrossMesh(m_pScene->ceiling_, ray);
    if (tempIntersectedObj.mesh_type == UNDEFINED)
        if (tempIntersectedObj.distance < IntersectedObj.distance || IntersectedObj.mesh_type == NONE)
        {
            IntersectedObj = tempIntersectedObj;
            IntersectedObj.mesh_type = CEILING;
        }
    tempIntersectedObj = CheckRayCrossMesh(m_pScene->left_wall_, ray);
    if (tempIntersectedObj.mesh_type == UNDEFINED)
        if (tempIntersectedObj.distance < IntersectedObj.distance || IntersectedObj.mesh_type == NONE)
        {
            IntersectedObj = tempIntersectedObj;
            IntersectedObj.mesh_type = LEFT_WALL;
        }
    tempIntersectedObj = CheckRayCrossMesh(m_pScene->right_wall_, ray);
    if (tempIntersectedObj.mesh_type == UNDEFINED)
        if (tempIntersectedObj.distance < IntersectedObj.distance || IntersectedObj.mesh_type == NONE)
        {
            IntersectedObj = tempIntersectedObj;
            IntersectedObj.mesh_type = RIGHT_WALL;
        }
    
    tempIntersectedObj = CheckRayCrossMesh(m_pScene->box_ceiling_, ray);
    if (tempIntersectedObj.mesh_type == UNDEFINED)
        if (tempIntersectedObj.distance < IntersectedObj.distance || IntersectedObj.mesh_type == NONE)
        {
            IntersectedObj = tempIntersectedObj;
            IntersectedObj.mesh_type = BOX_CEILING;
        }
    tempIntersectedObj = CheckRayCrossMesh(m_pScene->box_right_wall_, ray);
    if (tempIntersectedObj.mesh_type == UNDEFINED)
        if (tempIntersectedObj.distance < IntersectedObj.distance || IntersectedObj.mesh_type == NONE)
        {
            IntersectedObj = tempIntersectedObj;
            IntersectedObj.mesh_type = BOX_RIGHT_WALL;
        }
    tempIntersectedObj = CheckRayCrossMesh(m_pScene->box_left_wall_, ray);
    if (tempIntersectedObj.mesh_type == UNDEFINED)
        if (tempIntersectedObj.distance < IntersectedObj.distance || IntersectedObj.mesh_type == NONE)
        {
            IntersectedObj = tempIntersectedObj;
            IntersectedObj.mesh_type = BOX_LEFT_WALL;
        }
    tempIntersectedObj = CheckRayCrossMesh(m_pScene->box_front_wall_, ray);
    if (tempIntersectedObj.mesh_type == UNDEFINED)
        if (tempIntersectedObj.distance < IntersectedObj.distance || IntersectedObj.mesh_type == NONE)
        {
            IntersectedObj = tempIntersectedObj;
            IntersectedObj.mesh_type = BOX_FRONT_WALL;
        }
    tempIntersectedObj = CheckRayCrossSphere(m_pScene->sphere_, ray);
    if (tempIntersectedObj.mesh_type == UNDEFINED)
        if (tempIntersectedObj.distance < IntersectedObj.distance || IntersectedObj.mesh_type == NONE)
        {
            IntersectedObj = tempIntersectedObj;
            IntersectedObj.mesh_type = SPHERE;
        }
    tempIntersectedObj = CheckRayCrossMesh(m_pScene->light_, ray);
    if (tempIntersectedObj.mesh_type == UNDEFINED)
        if (tempIntersectedObj.distance < IntersectedObj.distance || IntersectedObj.mesh_type == NONE)
        {
            IntersectedObj = tempIntersectedObj;
            IntersectedObj.mesh_type = LIGHT;
        }

    vec3 color = ProcessRay(ray, IntersectedObj);
    return color;
}


SIntrsectionObj CTracer::CheckRayCrossMesh(SMesh& mesh, SRay& ray) //planes
{
    SIntrsectionObj IntersectedObj;
    
    for (uint i = 0; i < mesh.m_triangles.size(); i++)
    {
        uvec3 triangle = mesh.m_triangles[i];
        vec3 a = mesh.m_vertices[triangle.x];
        vec3 T = ray.m_start - a;
        vec3 b = mesh.m_vertices[triangle.y];
        vec3 c = mesh.m_vertices[triangle.z];
        vec3 E1 = b - a;
        vec3 E2 = c - a;
        T = ray.m_start - a;
        vec3 P = glm::cross(ray.m_dir, E2);
        float d = glm::dot(P, E1);
        float u = glm::dot(P, T) / d;
        if (u < 0.0f)
            continue;
        vec3 Q = glm::cross(T, E1);
        float v = glm::dot(Q, ray.m_dir) / d;
        if (v < 0.0f)
            continue;
        if (u + v > 1.0f)
            continue;
        float t = glm::dot(Q, E2) / d;
        if (t < 0.0f)
            continue;
        if (IntersectedObj.distance > t || IntersectedObj.mesh_type == NONE)
        {
            IntersectedObj.mesh_type = UNDEFINED;
            IntersectedObj.distance = t;
            IntersectedObj.triangle_num = i;
        }
    }

    return IntersectedObj;
}

SIntrsectionObj CTracer::CheckRayCrossSphere(SSphere& sphere, SRay& ray) // пересечение со сферой
{
    SIntrsectionObj IntersectedObj;
    glm::vec3 m = ray.m_dir;
    glm::vec3 dist = ray.m_start - sphere.centre;
    glm::vec3 rad = sphere.r;
    float a = glm::dot(m, m);
    float b = glm::dot(m, dist);
    float c = glm::dot(dist, dist) - glm::dot(rad, rad);
    float D = b * b - a * c;
    if (D < 0) {
        return IntersectedObj;
    }
    float t = (- b - sqrt(D)) / a;
    if (t < 0) {
        t = t = (- b + sqrt(D)) / a;
    }
    IntersectedObj.mesh_type = UNDEFINED;
    IntersectedObj.distance = t;
    IntersectedObj.triangle_num = 0;
    return IntersectedObj;
}

float glossScal(glm::vec3 norm,  glm::vec3 L, glm::vec3 V, int coeff)
{
    float gloss_scalar = 2 * glm::dot(norm, L) * glm::dot(norm, V) - glm::dot(L, V);
    return glm::dot(norm, L) < 0 ? 0 : pow(gloss_scalar, coeff);

}

vec3 CTracer::ProcessRay(SRay& ray, SIntrsectionObj& IntersectedObj)
{
    float gloss_scalar = 0;
    float light_spot = 0;
    vec3 color = vec3(0, 0, 0);
    glm::vec3 point = ray.m_start + ray.m_dir * IntersectedObj.distance;
    glm::vec3 norm = vec3 (0, 0, 0);
    glm::vec3 L = glm::normalize(m_light.m_pos - point);
    glm::vec3 V = - ray.m_dir;
    if (IntersectedObj.mesh_type == SPHERE) {
        norm = glm::normalize(point - m_pScene->sphere_.centre);
        gloss_scalar = glossScal(norm, L, V, 8);;
        color = vec3(0.2, 0.2, 0.4);
    }
    if (IntersectedObj.mesh_type == BACK_WALL) {
        norm = vec3(0, 0, 1);
        gloss_scalar = glossScal(norm, L, V,4);
        color = vec3(0.35, 0.33, 0);
    }
    if (IntersectedObj.mesh_type == FLOOR) { //calculating photon map
        norm = vec3(0, -1, 0);
        gloss_scalar = glossScal(norm, L, V,4);
        int i = (point.x + 40.0f) / 80.0f * m_photonMap.ResX;
        int j = (point.z + 40.0f) / 80.0f * m_photonMap.ResY;
        int r = m_photonMap.Rad;
        float sum = 0;
        for (int ic = - r; ic < r; ic++) {
            for (int jc = - r; jc < r; jc++) {
                if (ic * ic + jc * jc > r * r 
                || ic + i < 0 
                || jc + j < 0
                || ic + i + 1 >= m_photonMap.ResX
                || jc + j + 1 >= m_photonMap.ResY) {
                    continue;
                }
                sum += m_photonMap.m_pixels[(i + ic) * m_photonMap.ResY + (j + jc)];
            }
        }
        light_spot = sum / (M_PI * r * r);
        color = vec3(0, 0.4, 0.9) / 3.0f;
    }
    if (IntersectedObj.mesh_type == RIGHT_WALL) {
        norm = vec3(-1, 0, 0);
        gloss_scalar = glossScal(norm, L, V,4);
        color = vec3(0.2, 0.3, 0.1);
    }
    if (IntersectedObj.mesh_type == CEILING) {
        norm = vec3(0, 1, 0);
        color = vec3(0.6, 0.3, 0.2);
    }
    if (IntersectedObj.mesh_type == LEFT_WALL) {
        norm = vec3(1, 0, 0);
        gloss_scalar = glossScal(norm, L, V,4);
        color = vec3(0.7, 0, 0);
    }
    if (IntersectedObj.mesh_type == BOX_LEFT_WALL) {
        norm = vec3(-1, 0, 0);
        gloss_scalar = glossScal(norm, L, V,4);
        color = vec3(0.33, 0.33, 0.5);
    }
    if (IntersectedObj.mesh_type == BOX_RIGHT_WALL) {
        norm = vec3(1, 0, 0);
        gloss_scalar = glossScal(norm, L, V,4);
        color = vec3(0.33, 0.33, 0.5);
    }
    if (IntersectedObj.mesh_type == BOX_CEILING) {
        norm = vec3(0, -1, 0);
        gloss_scalar = glossScal(norm, L, V,4);
        color = vec3(0.33, 0.33, 0.5);
    }
    if (IntersectedObj.mesh_type == BOX_FRONT_WALL) {
        norm = vec3(0, 0, 1);
        gloss_scalar = glossScal(norm, L, V,4);
        color = vec3(0.33, 0.33, 0.5);
    }
    if (IntersectedObj.mesh_type == LIGHT) {
        norm = glm::vec3(0, 0, 1);
        color = vec3(0.1, 0.1, 0.1);
    }
    SRay ShadeRay;
    ShadeRay.m_start = m_light.m_pos;
    ShadeRay.m_dir = glm::normalize(point - m_light.m_pos);
    // shading
    SIntrsectionObj sphereShade = CheckRayCrossSphere (m_pScene->sphere_, ShadeRay);
    if (sphereShade.mesh_type == UNDEFINED && IntersectedObj.mesh_type != SPHERE) {
        if (sphereShade.distance < glm::distance(m_light.m_pos, point) && sphereShade.distance > 0) {
            return color * 0.189f;
        }
    }
    SIntrsectionObj boxRightWallShadow = CheckRayCrossMesh (m_pScene->box_right_wall_, ShadeRay);
    if (boxRightWallShadow.mesh_type == UNDEFINED && IntersectedObj.mesh_type != BOX_RIGHT_WALL) {
        if (boxRightWallShadow.distance < glm::distance(m_light.m_pos, point) && boxRightWallShadow.distance > 0) {
            return color * 0.189f;
        }
    }
    SIntrsectionObj boxLeftWallShadow = CheckRayCrossMesh (m_pScene->box_left_wall_, ShadeRay); 
    if (boxLeftWallShadow.mesh_type == UNDEFINED && IntersectedObj.mesh_type != BOX_LEFT_WALL) {
        if (boxLeftWallShadow.distance < glm::distance(m_light.m_pos, point) && boxLeftWallShadow.distance > 0) {
            return color * 0.189f;
        }
    }
    SIntrsectionObj boxFrontWallShadow = CheckRayCrossMesh (m_pScene->box_front_wall_, ShadeRay); 
    if (boxFrontWallShadow.mesh_type == UNDEFINED && IntersectedObj.mesh_type != BOX_FRONT_WALL) {
        if (boxFrontWallShadow.distance < glm::distance(m_light.m_pos, point) && boxFrontWallShadow.distance > 0) {
            return color * 0.189f;
        }
    }
    SIntrsectionObj boxCeilingShadow = CheckRayCrossMesh (m_pScene->box_ceiling_, ShadeRay);
    if (boxCeilingShadow.mesh_type == UNDEFINED && IntersectedObj.mesh_type != BOX_CEILING) {
        if (boxCeilingShadow.distance < glm::distance(m_light.m_pos, point) && boxCeilingShadow.distance > 0) {
            return color * 0.189f;
        }
    }
    //ends here
    float dot = glm::dot(L, norm) < 0 ? 0 : glm::dot(L, norm);
    return color * (2.0f * dot + 0.2f + gloss_scalar + light_spot);
}


float CTracer::FresnelReflectance(float IncidentCosine, float RI) {
    float ci2 = IncidentCosine * IncidentCosine;
    float si2 = 1.0f - ci2;
    float si4 = si2 * si2;
    float a = ci2 + RI * RI - 1.0f;
    float sqa = 2.0f * sqrtf(a) * IncidentCosine;
    float b = ci2 + a;
    float c = ci2 * a + si4;
    float d = sqa * si2;
    return (b - sqa) / (b + sqa) * (1.0f - d / (c + d));
}

void CTracer::lighting(SPhotonMap &m_points, int n, float s) {
    m_points.m_pixels.resize(n + sqrt(n) + 1);
    m_points.ResX = sqrt(n);
    m_points.ResY = sqrt(n);
    m_points.Rad = sqrt(n) / 100;
    for (int k = 0; k < n; k++) {
        m_points.m_pixels[k] = 0;
    }
    for (int k = 0; k < n; k ++) {
        glm::vec2 point_s = vec2(s * (rand() / (1.0f * RAND_MAX) - 0.5),
        s * (rand() / (1.0f * RAND_MAX) - 0.5));
        float sins = (2.0f * rand() / RAND_MAX - 1);
        float D = (m_lens.n - 1) * (1 / m_lens.R1 + 1 / m_lens.R2);
        float f = 1 / D; 
        glm::vec3 m1 = glm::cross(vec3(1, 0, 0), m_light.m_dir);
        glm::vec3 m2 = glm::cross(m1, m_light.m_dir);
        
        float df = length(point_s) / sins;
        float d = (f + df) / (D * (f + df) - 1);
        
        
        glm::vec2 lensPoint = point_s * (1 + f / df);
        
        if (length (m1 * lensPoint.x + m2 * lensPoint.y) > 10) {
            continue;
        }
        
        SRay light_ray; // постоение луча
        light_ray.m_start = (m_light.m_pos + m_light.m_dir * f + m1 * lensPoint.x + m2 * lensPoint.y);
        
        if (d < 1e-3) {
            continue;
        } else if (d > 10000) {
            light_ray.m_dir = m_light.m_dir;
        } else {
            light_ray.m_dir = glm::normalize(m_light.m_pos + m_light.m_dir * (f + d) - light_ray.m_start);
        }


        SIntrsectionObj floorspot = CheckRayCrossMesh(m_pScene->floor_, light_ray);
        glm::vec3 normN = vec3 (0, 1, 0);
        glm::vec3 intersection_point = light_ray.m_start + light_ray.m_dir * floorspot.distance;
        int i = (intersection_point.x + 40.0f) / 80.0f * m_points.ResX;
        int j = (intersection_point.z + 40.0f) / 80.0f * m_points.ResY;


        float lightedArea = glm::dot(normN, light_ray.m_dir) < 0 ? 0 : glm::dot(normN, light_ray.m_dir);
        m_points.m_pixels [i * m_points.ResX + j] += 2 * lightedArea *
            FresnelReflectance((float)(abs(glm::dot(light_ray.m_dir, m_light.m_dir))), m_lens.n - 1);
        
    }
}

void CTracer::RenderImage(int xRes, int yRes, int AA)
{
    // Rendering
    m_camera.m_up = glm::vec3(0, 1, 0);
    m_camera.m_forward = glm::normalize(m_camera.m_forward);
    m_camera.m_right = glm::cross(m_camera.m_forward, m_camera.m_up);
    
    right_width = m_camera.m_right * tan(m_camera.m_viewAngle.x / 2);
    up_height = m_camera.m_up * tan(m_camera.m_viewAngle.y / 2);
    
    m_camera.m_resolution = uvec2(xRes, yRes);
    m_camera.m_pixels.resize(xRes * yRes);
    
	SRay ray;
	vec3 color;
    if (AA == 1){
        for(uint i = 0; i < yRes; i++) {
        for(uint j = 0; j < xRes; j++) {
            m_camera.m_pixels[i * xRes + j] = vec3(0);
			ray = MakeRay(uvec2(j, i), vec2(0.25));
			color = TraceRay(ray);
			m_camera.m_pixels[i * xRes + j] += color;
			ray = MakeRay(uvec2(j, i), vec2(0.75));
			color = TraceRay(ray);
			m_camera.m_pixels[i * xRes + j] += color;
			ray = MakeRay(uvec2(j, i), vec2(0.25, 0.75));
			color = TraceRay(ray);
			m_camera.m_pixels[i * xRes + j] += color;
			ray = MakeRay(uvec2(j, i), vec2(0.75, 0.25));
			color = TraceRay(ray);
			m_camera.m_pixels[i * xRes + j] += color;
			ray = MakeRay(uvec2(j, i), vec2(0.5));
			color = TraceRay(ray);
			m_camera.m_pixels[i * xRes + j] += color;
			m_camera.m_pixels[i * xRes + j] *= vec3(0.2);
        }
    }

    }
    else{
        for(uint i = 0; i < yRes; i++) {
		for(uint j = 0; j < xRes; j++) {
			m_camera.m_pixels[i * xRes + j] = vec3(0);
			ray = MakeRay(uvec2(j, i), vec2(0.5));
			color = TraceRay(ray);
			m_camera.m_pixels[i * xRes + j] += color;
		}
    }
    }
}

void CTracer::SaveImageToFile(std::string fileName)
{
    BMP image;

    int width = m_camera.m_resolution.x;
    int height = m_camera.m_resolution.y;

    image.SetSize(width, height);

    RGBApixel p;

    int textureDisplacement = 0;

    for  (int i = 0; i < height; i++) {
        for  (int j = 0; j < width; j++) {
            vec3 color = m_camera.m_pixels[textureDisplacement + j];

            p.Red   = clamp(color.r, 0.0f, 1.0f) * 255.0f;
            p.Green = clamp(color.g, 0.0f, 1.0f) * 255.0f;
            p.Blue  = clamp(color.b, 0.0f, 1.0f) * 255.0f;

            image.SetPixel(j, i, p);
        }

        textureDisplacement += width;
    }

    image.WriteToFile(fileName.c_str());
}
