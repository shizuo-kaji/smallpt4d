#include "scene_definition.h"

#include <cmath>

namespace {

constexpr int kDiff = 0;
constexpr int kSpec = 1;
constexpr int kRefr = 2;

} // namespace

std::vector<SceneSphere> build_default_scene() {
    std::vector<SceneSphere> scene;
    scene.reserve(1000);

    scene.push_back(SceneSphere{1e5, Eigen::Vector4d(-1e5 - 175, 0, 0, 0), Eigen::Array3d::Zero(), Eigen::Array3d(.75, .25, .25), kSpec}); // Left wall
    scene.push_back(SceneSphere{1e5, Eigen::Vector4d(1e5 + 175, 0, 0, 0), Eigen::Array3d::Zero(), Eigen::Array3d(.25, .25, .75), kSpec}); // Right wall
    scene.push_back(SceneSphere{1e5, Eigen::Vector4d(0, 0, 1e5 - 175, 0), Eigen::Array3d::Zero(), Eigen::Array3d(.25, .75, .25), kSpec}); // Back wall
    scene.push_back(SceneSphere{1e5, Eigen::Vector4d(0, 0, -1e5 + 200, 0), Eigen::Array3d::Zero(), Eigen::Array3d(.75, .75, .75), kSpec}); // Front wall
    scene.push_back(SceneSphere{1e5, Eigen::Vector4d(0, -1e5 - 175, 0, 0), Eigen::Array3d(1, 1, 1), Eigen::Array3d(.35, .75, .55), kSpec}); // Floor light
    scene.push_back(SceneSphere{1e5, Eigen::Vector4d(0, 1e5 + 175, 0, 0), Eigen::Array3d(0, 0, 0), Eigen::Array3d(.35, .25, .55), kSpec}); // Ceiling
    scene.push_back(SceneSphere{1e5, Eigen::Vector4d(0, 0, 0, 175 + 1e5), Eigen::Array3d::Zero(), Eigen::Array3d(.35, .55, .75), kSpec}); // W+
    scene.push_back(SceneSphere{1e5, Eigen::Vector4d(0, 0, 0, -175 - 1e5), Eigen::Array3d::Zero(), Eigen::Array3d(.75, .55, .35), kSpec}); // W-

    const int m = 2;
    for (int i = -m; i <= m; ++i) {
        for (int j = -m; j <= m; ++j) {
            for (int k = -m; k <= m; ++k) {
                for (int l = -m; l <= m; ++l) {
                    const Eigen::Vector4d center(i * 75 / m, j * 75 / m, k * 75 / m, l * 75 / m);
                    const Eigen::Array3d tint(.50 + .2 * i / m, .50 + .2 * j / m, .50 + .2 * l / m);
                    switch (std::abs(i + j + k + l) % 3) {
                        case 0:
                            scene.push_back(SceneSphere{10, center, tint, tint, kDiff});
                            break;
                        case 1:
                            scene.push_back(SceneSphere{10, center, Eigen::Array3d::Zero(), tint, kSpec});
                            break;
                        default:
                            scene.push_back(SceneSphere{10, center, Eigen::Array3d::Zero(), tint, kRefr});
                            break;
                    }
                }
            }
        }
    }

    return scene;
}

