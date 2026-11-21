#pragma once

#include <string>
#include <vector>

struct OpenClFloat4 {
    double x;
    double y;
    double z;
    double w;
};

struct OpenClCamera {
    OpenClFloat4 origin;
    OpenClFloat4 corner;
    OpenClFloat4 u;
    OpenClFloat4 v;
};

struct OpenClSphere {
    OpenClFloat4 center;
    OpenClFloat4 emission;
    OpenClFloat4 colour;
    int refl;
    double radius;
};

bool render_opencl_image(
    const std::string& output_path,
    int width,
    int height,
    int sample_count,
    int min_bounce,
    double env_map,
    int frame_seed,
    const OpenClCamera& camera,
    const std::vector<OpenClSphere>& scene,
    unsigned progress_bar_width = 50
);
