#pragma once

#include "Types.h"
#include "glm.hpp"
#include "vector"

class CScene
{
public:
    // Set of meshes
    SMesh ceiling_;
    SMesh floor_;
    SMesh left_wall_;
    SMesh right_wall_;
    SMesh back_wall_;
    SMesh box_ceiling_;
    SMesh box_right_wall_;
    SMesh box_left_wall_;
    SMesh box_front_wall_;
    SSphere sphere_;
    SMesh light_;
public:
    CScene(glm::vec3 LightPos, glm::vec3 LightDir);
};
