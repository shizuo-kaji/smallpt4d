#include "opencl_renderer.h"

#include "png_writer.h"
#include "progress_bar.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#if defined(SMALLPT4D_WITH_OPENCL)
#if defined(__APPLE__)
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif
#endif

namespace {

inline double clamp(double x) {
    return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x);
}

inline uint8_t to_byte(double x) {
    return static_cast<uint8_t>(std::pow(clamp(x), 1.0 / 2.2) * 255.0 + 0.5);
}

struct RGBs {
    RGBs() = default;
    RGBs(uint8_t r_, uint8_t g_, uint8_t b_) : r(r_), g(g_), b(b_) {}
    uint8_t r;
    uint8_t g;
    uint8_t b;
};

static_assert(sizeof(RGBs) == 3, "RGB pixel must be tightly packed.");

#if defined(SMALLPT4D_WITH_OPENCL)

struct alignas(16) SphereCLF {
    cl_float4 center;
    cl_float4 emission;
    cl_float4 colour;
    cl_int refl;
    cl_float radius;
    cl_int pad0;
    cl_int pad1;
};

struct alignas(32) SphereCLD {
    cl_double4 center;
    cl_double4 emission;
    cl_double4 colour;
    cl_int refl;
    cl_double radius;
    cl_int pad0;
    cl_int pad1;
};

static_assert(sizeof(SphereCLF) == 64, "SphereCLF layout mismatch.");
static_assert(sizeof(SphereCLD) == 128, "SphereCLD layout mismatch.");

struct Float4 {
    float x;
    float y;
    float z;
    float w;
};

constexpr const char* kOpenClKernel = R"CLC(
#if defined(USE_FP64)
#if defined(cl_khr_fp64)
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
#error "This OpenCL device does not support fp64."
#endif
#endif

#if defined(USE_FP64)
typedef double real;
typedef double3 real3;
typedef double4 real4;
#else
typedef float real;
typedef float3 real3;
typedef float4 real4;
#endif

typedef struct {
    real4 center;
    real4 emission;
    real4 colour;
    int refl;
    real radius;
    int pad0;
    int pad1;
} SphereCL;

uint xorshift32(uint x) {
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return x;
}

uint mix_seed(uint x) {
    x ^= x >> 16;
    x *= 0x7feb352du;
    x ^= x >> 15;
    x *= 0x846ca68bu;
    x ^= x >> 16;
    return x;
}

real rand01(uint* state) {
    *state = xorshift32(*state);
    return (real)(*state & 0x00FFFFFFu) * ((real)1.0 / (real)16777216.0);
}

real4 random_unit4(uint* state) {
    while (1) {
        real4 v = (real4)(
            rand01(state) * (real)2.0 - (real)1.0,
            rand01(state) * (real)2.0 - (real)1.0,
            rand01(state) * (real)2.0 - (real)1.0,
            rand01(state) * (real)2.0 - (real)1.0
        );
        real n2 = dot(v, v);
        if (n2 > (real)1e-12 && n2 <= (real)1.0) {
            return normalize(v);
        }
    }
}

int hit_scene(
    __global const SphereCL* scene,
    int scene_size,
    real4 ro,
    real4 rd,
    real t_min,
    real t_max,
    real* t_hit
) {
    int hit_idx = -1;
    real closest = t_max;
    for (int i = 0; i < scene_size; ++i) {
        real4 oc = scene[i].center - ro;
        real b = dot(oc, rd);
        real discriminant = b * b - dot(oc, oc) + scene[i].radius * scene[i].radius;
        if (discriminant > (real)0.0) {
            real s = sqrt(discriminant);
            real t = b - s;
            if (t > t_min && t < closest) {
                closest = t;
                hit_idx = i;
            } else {
                t = b + s;
                if (t > t_min && t < closest) {
                    closest = t;
                    hit_idx = i;
                }
            }
        }
    }
    *t_hit = closest;
    return hit_idx;
}

__kernel void render_kernel(
    __global const SphereCL* scene,
    int scene_size,
    int width,
    int height,
    int sample_count,
    int min_bounce,
    real env_map,
    int max_depth,
    real4 cam_origin,
    real4 cam_corner,
    real4 cam_u,
    real4 cam_v,
    uint frame_seed,
    real ray_t_min,
    int row_offset,
    __global float4* output
) {
    int x = get_global_id(0);
    int local_y = get_global_id(1);
    int y = row_offset + local_y;
    if (x >= width || y >= height) {
        return;
    }

    int pixel_idx = y * width + x;
    uint state = mix_seed(frame_seed ^ ((uint)pixel_idx * 747796405u + 2891336453u));

    typedef struct {
        real4 ro;
        real4 rd;
        real3 throughput;
        int depth;
    } PathState;

    real3 accum = (real3)((real)0.0, (real)0.0, (real)0.0);
    for (int s = 0; s < sample_count; ++s) {
        real u = ((real)x + rand01(&state)) / (real)width;
        real v = ((real)y + rand01(&state)) / (real)height;

        PathState stack[8];
        int stack_size = 1;
        stack[0].ro = cam_origin;
        stack[0].rd = normalize(cam_corner + u * cam_u + v * cam_v - cam_origin);
        stack[0].throughput = (real3)((real)1.0, (real)1.0, (real)1.0);
        stack[0].depth = 0;

        while (stack_size > 0) {
            PathState path = stack[--stack_size];
            real4 ro = path.ro;
            real4 rd = path.rd;
            real3 throughput = path.throughput;
            int depth = path.depth;

            while (depth < max_depth) {
                real t = (real)1e20;
                int hit_idx = hit_scene(scene, scene_size, ro, rd, ray_t_min, (real)1e20, &t);
                if (hit_idx < 0) {
                    real alpha = (real)0.5 * (rd.y + (real)1.0);
                    real3 bg = env_map * ((real3)((real)1.0, (real)1.0, (real)1.0) * ((real)1.0 - alpha)
                                         + (real3)((real)0.5, (real)0.7, (real)1.0) * alpha);
                    accum += throughput * bg;
                    break;
                }

                SphereCL sp = scene[hit_idx];
                real4 hitpos = ro + t * rd;
                real4 normal = normalize(hitpos - sp.center);
                real4 nl = dot(normal, rd) < (real)0.0 ? normal : -normal;
                real3 f = sp.colour.xyz;
                real p = fmax(f.x, fmax(f.y, f.z));

                if ((depth + 1 > min_bounce) || p <= (real)0.0) {
                    if (p <= (real)0.0 || rand01(&state) >= p) {
                        accum += throughput * sp.emission.xyz;
                        break;
                    }
                    f /= p;
                }

                accum += throughput * sp.emission.xyz;

                if (sp.refl == 0) {
                    real4 d = random_unit4(&state);
                    if (dot(d, nl) < (real)0.0) {
                        d = -d;
                    }
                    ro = hitpos;
                    rd = d;
                    throughput *= f;
                    depth += 1;
                    continue;
                }

                if (sp.refl == 1) {
                    ro = hitpos;
                    rd = normalize(rd - (real)2.0 * dot(normal, rd) * normal);
                    throughput *= f;
                    depth += 1;
                    continue;
                }

                real4 refl = normalize(rd - (real)2.0 * dot(normal, rd) * normal);
                int into = dot(normal, nl) > (real)0.0;
                real nc = (real)1.0;
                real nt = (real)1.5;
                real nnt = into ? (nc / nt) : (nt / nc);
                real ddn = dot(rd, nl);
                real cos2t = (real)1.0 - nnt * nnt * ((real)1.0 - ddn * ddn);
                if (cos2t < (real)0.0) {
                    ro = hitpos;
                    rd = refl;
                    throughput *= f;
                    depth += 1;
                    continue;
                }

                real4 tdir = normalize(rd * nnt - normal * ((into ? (real)1.0 : (real)-1.0) * (ddn * nnt + sqrt(cos2t))));
                real a = nt - nc;
                real b = nt + nc;
                real R0 = (a * a) / (b * b);
                real c = (real)1.0 - (into ? -ddn : dot(tdir, normal));
                real Re = R0 + ((real)1.0 - R0) * c * c * c * c * c;
                real Tr = (real)1.0 - Re;

                if (depth <= 1) {
                    if (stack_size < 8) {
                        PathState trans_path;
                        trans_path.ro = hitpos;
                        trans_path.rd = tdir;
                        trans_path.throughput = throughput * f * Tr;
                        trans_path.depth = depth + 1;
                        stack[stack_size++] = trans_path;
                    }
                    ro = hitpos;
                    rd = refl;
                    throughput *= f * Re;
                    depth += 1;
                    continue;
                }

                real P = (real)0.25 + (real)0.5 * Re;
                if (rand01(&state) < P) {
                    ro = hitpos;
                    rd = refl;
                    throughput *= f * (Re / P);
                } else {
                    ro = hitpos;
                    rd = tdir;
                    throughput *= f * (Tr / ((real)1.0 - P));
                }
                depth += 1;
            }
        }
    }

    real inv_samples = (real)1.0 / (real)sample_count;
    real3 out_color = accum * inv_samples;
    output[pixel_idx] = (float4)((float)out_color.x, (float)out_color.y, (float)out_color.z, 1.0f);
}
)CLC";

bool check_cl(cl_int err, const char* step) {
    if (err == CL_SUCCESS) {
        return true;
    }
    std::cerr << "OpenCL error (" << err << ") at " << step << std::endl;
    return false;
}

struct OpenClState {
    cl_device_id device = nullptr;
    cl_context context = nullptr;
    cl_command_queue queue = nullptr;
    cl_program program = nullptr;
    cl_kernel kernel = nullptr;
    bool use_fp64 = false;
    bool initialized = false;
    bool available = false;

    ~OpenClState() {
        if (kernel) clReleaseKernel(kernel);
        if (program) clReleaseProgram(program);
        if (queue) clReleaseCommandQueue(queue);
        if (context) clReleaseContext(context);
    }

    bool init() {
        if (initialized) {
            return available;
        }
        initialized = true;

        cl_uint platform_count = 0;
        cl_int err = clGetPlatformIDs(0, nullptr, &platform_count);
        if (!check_cl(err, "clGetPlatformIDs(count)") || platform_count == 0) {
            return false;
        }

        std::vector<cl_platform_id> platforms(platform_count);
        err = clGetPlatformIDs(platform_count, platforms.data(), nullptr);
        if (!check_cl(err, "clGetPlatformIDs(list)")) {
            return false;
        }

        for (cl_platform_id platform : platforms) {
            cl_uint gpu_count = 0;
            if (clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, nullptr, &gpu_count) == CL_SUCCESS && gpu_count > 0) {
                std::vector<cl_device_id> devices(gpu_count);
                if (clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, gpu_count, devices.data(), nullptr) == CL_SUCCESS) {
                    device = devices[0];
                    break;
                }
            }
        }
        if (!device) {
            for (cl_platform_id platform : platforms) {
                cl_uint dev_count = 0;
                if (clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 0, nullptr, &dev_count) == CL_SUCCESS && dev_count > 0) {
                    std::vector<cl_device_id> devices(dev_count);
                    if (clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, dev_count, devices.data(), nullptr) == CL_SUCCESS) {
                        device = devices[0];
                        break;
                    }
                }
            }
        }
        if (!device) {
            std::cerr << "OpenCL device not found." << std::endl;
            return false;
        }

        size_t ext_size = 0;
        err = clGetDeviceInfo(device, CL_DEVICE_EXTENSIONS, 0, nullptr, &ext_size);
        if (!check_cl(err, "clGetDeviceInfo(CL_DEVICE_EXTENSIONS size)")) {
            return false;
        }
        std::vector<char> extensions(ext_size + 1, '\0');
        err = clGetDeviceInfo(device, CL_DEVICE_EXTENSIONS, ext_size, extensions.data(), nullptr);
        if (!check_cl(err, "clGetDeviceInfo(CL_DEVICE_EXTENSIONS value)")) {
            return false;
        }
        const std::string ext_str(extensions.data());
        use_fp64 = (ext_str.find("cl_khr_fp64") != std::string::npos) ||
                   (ext_str.find("cl_amd_fp64") != std::string::npos);
        if (!use_fp64) {
            std::cerr << "OpenCL device does not support fp64. Using float32 fallback." << std::endl;
        }

        context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &err);
        if (!check_cl(err, "clCreateContext")) {
            return false;
        }

        queue = clCreateCommandQueue(context, device, 0, &err);
        if (!check_cl(err, "clCreateCommandQueue")) {
            return false;
        }

        const char* source = kOpenClKernel;
        const size_t source_len = std::strlen(source);
        program = clCreateProgramWithSource(context, 1, &source, &source_len, &err);
        if (!check_cl(err, "clCreateProgramWithSource")) {
            return false;
        }

        const char* build_opts = use_fp64 ? "-DUSE_FP64=1" : "";
        err = clBuildProgram(program, 1, &device, build_opts, nullptr, nullptr);
        if (err != CL_SUCCESS) {
            size_t log_size = 0;
            clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &log_size);
            std::vector<char> log(log_size + 1, '\0');
            clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, log_size, log.data(), nullptr);
            std::cerr << "OpenCL program build failed:\n" << log.data() << std::endl;
            check_cl(err, "clBuildProgram");
            return false;
        }

        kernel = clCreateKernel(program, "render_kernel", &err);
        if (!check_cl(err, "clCreateKernel")) {
            return false;
        }

        available = true;
        return true;
    }
};

#endif

} // namespace

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
    unsigned progress_bar_width
) {
#if !defined(SMALLPT4D_WITH_OPENCL)
    (void)output_path;
    (void)width;
    (void)height;
    (void)sample_count;
    (void)min_bounce;
    (void)env_map;
    (void)frame_seed;
    (void)camera;
    (void)scene;
    (void)progress_bar_width;
    return false;
#else
    if (width <= 0 || height <= 0 || sample_count <= 0 || scene.empty()) {
        return false;
    }

    static OpenClState state;
    if (!state.init()) {
        return false;
    }

    const std::size_t pixel_count = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    std::vector<Float4> out_pixels(pixel_count);

    cl_int err = CL_SUCCESS;
    cl_mem cl_scene = nullptr;

    cl_mem cl_output = clCreateBuffer(
        state.context,
        CL_MEM_WRITE_ONLY,
        out_pixels.size() * sizeof(Float4),
        nullptr,
        &err
    );
    if (!check_cl(err, "clCreateBuffer(output)")) {
        return false;
    }

    const double fallback_scene_scale = 1e-4;
    const double scene_scale = state.use_fp64 ? 1.0 : fallback_scene_scale;
    const int scene_size = static_cast<int>(scene.size());
    const int max_depth = 64;
    const cl_uint cl_seed = static_cast<cl_uint>(frame_seed * 7919 + 12345);

    const int arg_scene = 0;
    const int arg_scene_size = 1;
    const int arg_width = 2;
    const int arg_height = 3;
    const int arg_sample_count = 4;
    const int arg_min_bounce = 5;
    const int arg_env_map = 6;
    const int arg_max_depth = 7;
    const int arg_cam_origin = 8;
    const int arg_cam_corner = 9;
    const int arg_cam_u = 10;
    const int arg_cam_v = 11;
    const int arg_seed = 12;
    const int arg_ray_t_min = 13;
    const int arg_row_offset = 14;
    const int arg_output = 15;

    err  = clSetKernelArg(state.kernel, arg_scene_size, sizeof(int), &scene_size);
    err |= clSetKernelArg(state.kernel, arg_width, sizeof(int), &width);
    err |= clSetKernelArg(state.kernel, arg_height, sizeof(int), &height);
    err |= clSetKernelArg(state.kernel, arg_sample_count, sizeof(int), &sample_count);
    err |= clSetKernelArg(state.kernel, arg_min_bounce, sizeof(int), &min_bounce);
    err |= clSetKernelArg(state.kernel, arg_max_depth, sizeof(int), &max_depth);
    err |= clSetKernelArg(state.kernel, arg_seed, sizeof(cl_uint), &cl_seed);
    err |= clSetKernelArg(state.kernel, arg_output, sizeof(cl_mem), &cl_output);
    if (!check_cl(err, "clSetKernelArg(static)")) {
        clReleaseMemObject(cl_output);
        return false;
    }

    if (state.use_fp64) {
        std::vector<SphereCLD> scene_buffer;
        scene_buffer.reserve(scene.size());
        for (const OpenClSphere& sp : scene) {
            SphereCLD out{};
            out.center = cl_double4{
                sp.center.x * scene_scale, sp.center.y * scene_scale,
                sp.center.z * scene_scale, sp.center.w * scene_scale
            };
            out.emission = cl_double4{sp.emission.x, sp.emission.y, sp.emission.z, sp.emission.w};
            out.colour = cl_double4{sp.colour.x, sp.colour.y, sp.colour.z, sp.colour.w};
            out.refl = static_cast<cl_int>(sp.refl);
            out.radius = sp.radius * scene_scale;
            scene_buffer.push_back(out);
        }

        cl_scene = clCreateBuffer(
            state.context,
            CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
            scene_buffer.size() * sizeof(SphereCLD),
            scene_buffer.data(),
            &err
        );
        if (!check_cl(err, "clCreateBuffer(scene fp64)")) {
            clReleaseMemObject(cl_output);
            return false;
        }

        const cl_double4 cam_origin = cl_double4{
            camera.origin.x * scene_scale, camera.origin.y * scene_scale,
            camera.origin.z * scene_scale, camera.origin.w * scene_scale
        };
        const cl_double4 cam_corner = cl_double4{
            camera.corner.x * scene_scale, camera.corner.y * scene_scale,
            camera.corner.z * scene_scale, camera.corner.w * scene_scale
        };
        const cl_double4 cam_u = cl_double4{
            camera.u.x * scene_scale, camera.u.y * scene_scale,
            camera.u.z * scene_scale, camera.u.w * scene_scale
        };
        const cl_double4 cam_v = cl_double4{
            camera.v.x * scene_scale, camera.v.y * scene_scale,
            camera.v.z * scene_scale, camera.v.w * scene_scale
        };
        const double env_map64 = env_map;
        const double ray_t_min64 = 1e-4;

        err  = clSetKernelArg(state.kernel, arg_scene, sizeof(cl_mem), &cl_scene);
        err |= clSetKernelArg(state.kernel, arg_env_map, sizeof(double), &env_map64);
        err |= clSetKernelArg(state.kernel, arg_cam_origin, sizeof(cl_double4), &cam_origin);
        err |= clSetKernelArg(state.kernel, arg_cam_corner, sizeof(cl_double4), &cam_corner);
        err |= clSetKernelArg(state.kernel, arg_cam_u, sizeof(cl_double4), &cam_u);
        err |= clSetKernelArg(state.kernel, arg_cam_v, sizeof(cl_double4), &cam_v);
        err |= clSetKernelArg(state.kernel, arg_ray_t_min, sizeof(double), &ray_t_min64);
        if (!check_cl(err, "clSetKernelArg(fp64)")) {
            clReleaseMemObject(cl_output);
            clReleaseMemObject(cl_scene);
            return false;
        }
    } else {
        std::vector<SphereCLF> scene_buffer;
        scene_buffer.reserve(scene.size());
        for (const OpenClSphere& sp : scene) {
            SphereCLF out{};
            out.center = cl_float4{
                static_cast<cl_float>(sp.center.x * scene_scale),
                static_cast<cl_float>(sp.center.y * scene_scale),
                static_cast<cl_float>(sp.center.z * scene_scale),
                static_cast<cl_float>(sp.center.w * scene_scale)
            };
            out.emission = cl_float4{
                static_cast<cl_float>(sp.emission.x),
                static_cast<cl_float>(sp.emission.y),
                static_cast<cl_float>(sp.emission.z),
                static_cast<cl_float>(sp.emission.w)
            };
            out.colour = cl_float4{
                static_cast<cl_float>(sp.colour.x),
                static_cast<cl_float>(sp.colour.y),
                static_cast<cl_float>(sp.colour.z),
                static_cast<cl_float>(sp.colour.w)
            };
            out.refl = static_cast<cl_int>(sp.refl);
            out.radius = static_cast<cl_float>(sp.radius * scene_scale);
            out.pad0 = 0;
            out.pad1 = 0;
            scene_buffer.push_back(out);
        }

        cl_scene = clCreateBuffer(
            state.context,
            CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
            scene_buffer.size() * sizeof(SphereCLF),
            scene_buffer.data(),
            &err
        );
        if (!check_cl(err, "clCreateBuffer(scene fp32)")) {
            clReleaseMemObject(cl_output);
            return false;
        }

        const cl_float4 cam_origin = cl_float4{
            static_cast<cl_float>(camera.origin.x * scene_scale),
            static_cast<cl_float>(camera.origin.y * scene_scale),
            static_cast<cl_float>(camera.origin.z * scene_scale),
            static_cast<cl_float>(camera.origin.w * scene_scale)
        };
        const cl_float4 cam_corner = cl_float4{
            static_cast<cl_float>(camera.corner.x * scene_scale),
            static_cast<cl_float>(camera.corner.y * scene_scale),
            static_cast<cl_float>(camera.corner.z * scene_scale),
            static_cast<cl_float>(camera.corner.w * scene_scale)
        };
        const cl_float4 cam_u = cl_float4{
            static_cast<cl_float>(camera.u.x * scene_scale),
            static_cast<cl_float>(camera.u.y * scene_scale),
            static_cast<cl_float>(camera.u.z * scene_scale),
            static_cast<cl_float>(camera.u.w * scene_scale)
        };
        const cl_float4 cam_v = cl_float4{
            static_cast<cl_float>(camera.v.x * scene_scale),
            static_cast<cl_float>(camera.v.y * scene_scale),
            static_cast<cl_float>(camera.v.z * scene_scale),
            static_cast<cl_float>(camera.v.w * scene_scale)
        };
        const float env_map32 = static_cast<float>(env_map);
        const float ray_t_min32 = 1e-6f;

        err  = clSetKernelArg(state.kernel, arg_scene, sizeof(cl_mem), &cl_scene);
        err |= clSetKernelArg(state.kernel, arg_env_map, sizeof(float), &env_map32);
        err |= clSetKernelArg(state.kernel, arg_cam_origin, sizeof(cl_float4), &cam_origin);
        err |= clSetKernelArg(state.kernel, arg_cam_corner, sizeof(cl_float4), &cam_corner);
        err |= clSetKernelArg(state.kernel, arg_cam_u, sizeof(cl_float4), &cam_u);
        err |= clSetKernelArg(state.kernel, arg_cam_v, sizeof(cl_float4), &cam_v);
        err |= clSetKernelArg(state.kernel, arg_ray_t_min, sizeof(float), &ray_t_min32);
        if (!check_cl(err, "clSetKernelArg(fp32)")) {
            clReleaseMemObject(cl_output);
            clReleaseMemObject(cl_scene);
            return false;
        }
    }

    const int progress_step = std::max(1, height / 10);
    const unsigned progress_ticks =
        static_cast<unsigned>((height + progress_step - 1) / progress_step);
    ProgressBar pb(progress_bar_width, progress_ticks);

    for (int row_offset = 0; row_offset < height; row_offset += progress_step) {
        const int rows_this_pass = std::min(progress_step, height - row_offset);
        err = clSetKernelArg(state.kernel, arg_row_offset, sizeof(int), &row_offset);
        if (!check_cl(err, "clSetKernelArg(row_offset)")) {
            clReleaseMemObject(cl_output);
            clReleaseMemObject(cl_scene);
            return false;
        }

        const size_t global_size[2] = {
            static_cast<size_t>(width),
            static_cast<size_t>(rows_this_pass)
        };
        err = clEnqueueNDRangeKernel(
            state.queue,
            state.kernel,
            2,
            nullptr,
            global_size,
            nullptr,
            0,
            nullptr,
            nullptr
        );
        if (!check_cl(err, "clEnqueueNDRangeKernel")) {
            clReleaseMemObject(cl_output);
            clReleaseMemObject(cl_scene);
            return false;
        }

        err = clFinish(state.queue);
        if (!check_cl(err, "clFinish")) {
            clReleaseMemObject(cl_output);
            clReleaseMemObject(cl_scene);
            return false;
        }

        pb.tick();
        pb.print();
    }
    pb.end();

    err = clEnqueueReadBuffer(
        state.queue,
        cl_output,
        CL_TRUE,
        0,
        out_pixels.size() * sizeof(Float4),
        out_pixels.data(),
        0,
        nullptr,
        nullptr
    );
    clReleaseMemObject(cl_output);
    clReleaseMemObject(cl_scene);
    if (!check_cl(err, "clEnqueueReadBuffer")) {
        return false;
    }

    std::vector<RGBs> bitmap(pixel_count);
    for (std::size_t i = 0; i < pixel_count; ++i) {
        bitmap[i] = RGBs(
            to_byte(out_pixels[i].x),
            to_byte(out_pixels[i].y),
            to_byte(out_pixels[i].z)
        );
    }

    if (!simple_png::write_png_rgb8(
            output_path,
            width,
            height,
            reinterpret_cast<const uint8_t*>(bitmap.data()),
            3 * width)) {
        std::cerr << "Failed to write image: " << output_path << std::endl;
        return false;
    }

    return true;
#endif
}
