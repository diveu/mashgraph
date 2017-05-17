#include "Tracer.h"
#include "stdio.h"
#include <stdio.h>
#include <time.h>


int main(int argc, char **argv)
{
    CTracer tracer;
    clock_t tStart = clock();
    int xRes = 1024;
    int yRes = 1024;
    tracer.m_camera.m_pos     = glm::vec3(0, 0, 0);
    tracer.m_camera.m_forward = glm::vec3(0, 1, 0);
    tracer.m_camera.m_viewAngle = glm::vec2(0.48, 0.48);
    tracer.m_light.m_pos = glm::vec3(-35, -25, -20);
    tracer.m_light.m_dir = glm::normalize(glm::vec3(0.3, 1, 0));
    tracer.m_lens.R1 = 10;
    tracer.m_lens.R2 = 10;
    tracer.m_lens.n = 1.6;
    int AA = 0;

    if(argc == 3) {
        FILE* file = fopen(argv[1], "r");

        if(file) {
            int xFile, yFile;
            float x = 0, y = 0, z = 0;
            if(fscanf(file, "%d %d", &xFile, &yFile) == 2) {
                xRes = xFile;
                yRes = yFile;
            } else
                printf("Invalid config format! Using default parameters.\r\n");
            if (fscanf(file, "%d", &AA) == 1){
                printf("AA=%d\n", AA);
            }
            else{
                printf("AA disabled\n");
            }
            if(fscanf(file, "%f %f %f", &x, &y, &z) == 3) {
                tracer.m_camera.m_pos = glm::vec3(x, y, z);
            } else
                printf("Invalid config format! Using default parameters.\r\n");
            
            if(fscanf(file, "%f %f %f", &x, &y, &z) == 3) {
                tracer.m_camera.m_forward = glm::vec3(x, y, z);
            } else
                printf("Invalid config format! Using default parameters.\r\n");
            
            if(fscanf(file, "%f %f", &x, &y) == 2) {
                tracer.m_camera.m_viewAngle = glm::vec2(x, y);
            } else
                printf("Invalid config format! Using default parameters.\r\n");
            
            if(fscanf(file, "%f %f %f", &x, &y, &z) == 3) {
                tracer.m_light.m_pos = glm::vec3(x, y, z);
            } else
                printf("Invalid config format! Using default parameters.\r\n");
            
            if(fscanf(file, "%f %f %f", &x, &y, &z) == 3) {
                tracer.m_light.m_dir = glm::normalize(glm::vec3(x, y, z));
            } else
                printf("Invalid config format! Using default parameters.\r\n");
            if(fscanf(file, "%f %f %f", &x, &y, &z) == 3) {
                tracer.m_lens.R1 = x;
                tracer.m_lens.R2 = y;
                tracer.m_lens.n = z;
            } else
                printf("Invalid config format! Using default parameters.\r\n");

            fclose(file);
        } else
            printf("Invalid config path! Using default parameters.\r\n");
    } else
        printf("No config! Using default parameters.\r\n");
    
    CScene scene(tracer.m_light.m_pos, glm::normalize(tracer.m_light.m_dir));
    tracer.m_pScene = &scene;
    tracer.lighting(tracer.m_photonMap, 5000000, 3.0f);
    tracer.RenderImage(xRes, yRes, AA);
    tracer.SaveImageToFile(argv[2]);
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
}
