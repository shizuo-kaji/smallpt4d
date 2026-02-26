#pragma once

#include <vector>

#include <Eigen/Dense>

struct SceneSphere {
    double radius;
    Eigen::Vector4d center;
    Eigen::Array3d emission;
    Eigen::Array3d colour;
    int refl;
};

std::vector<SceneSphere> build_default_scene();

