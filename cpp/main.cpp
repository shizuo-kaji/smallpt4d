/*
 smallpt4D: simple ray tracer for 4D scenes by Shizuo Kaji
 Based on: smallpt: Global Illumination in 99 lines of C++ by Kevin Beason
 Depends on: Eigen 3 (http://eigen.tuxfamily.org)
*/


#include "stdafx.h"
#include "opencl_renderer.h"
#include "png_writer.h"
#include "progress_bar.h"
#include "scene_definition.h"

#if defined(_WIN32)|| defined(_WIN64)
	#include "getopt_win.h"
#else
	#include <getopt.h>
#endif
#include <limits>

using namespace Eigen;

// configuration parameters
struct Configuration{
	// oculus' resolution is 1080x1200 for each eye
    int width = 1080;
    int height = 1200;
	int frames = 200;
	int min_bounce = 5;
    int number_of_samples = 100;   // default 100
    int emitter_sampling_off = 1;   // set 0 to sample emitters
    std::string output="output.png";
    double env_map = 1;
    bool use_gpu = false;
};

// random number generator
class Sampling{
public:
    explicit Sampling(int s): gen(s), zero_one(0.0, 1.0){}
    double random_double() { return zero_one(gen); }
    double random_signed() { return 2.0 * random_double() - 1.0; }
private:
    std::mt19937 gen;
    std::uniform_real_distribution<> zero_one;
};

inline Vector4d random_unit_vector(Sampling& sampler) {
    while (true) {
        const Vector4d candidate(
            sampler.random_signed(),
            sampler.random_signed(),
            sampler.random_signed(),
            sampler.random_signed()
        );
        const double norm2 = candidate.squaredNorm();
        if (norm2 > 1e-12 && norm2 <= 1.0) {
            return candidate.normalized();
        }
    }
}


struct Ray {
public:
    Ray(const Vector4d& o_, const Vector4d& d_, double t_max=std::numeric_limits<double>::infinity()) :
        o(o_), d(d_.normalized()), t_min(1e-4), t_max(t_max) {}
    Vector4d o, d;  // origin & direction
    double t_min, t_max;  // max/min distance for hit test
};

class Camera {
public:
    Camera(){}
    Camera(Vector4d origin, // origin of the camera
		Vector4d rel, // centre of the 2D screen relative to origin
		Vector4d u, Vector4d v): // coordinate vectors of the 2D screen
        origin(origin), upper_left_corner(origin+rel-u/2.0-v/2.0), u(u), v(v) { }
    
    Ray get_ray(double x, double y) const { // ray emitted from the origin towards (u,v) on the screen
        return Ray(origin, upper_left_corner + x * u + y * v - origin);
    }
    Vector4d origin;
    Vector4d upper_left_corner;
    Vector4d u,v;
};



enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()

class Sphere {
public:
    Sphere(double radius, const Vector4d& center, const Array3d& emission, const Array3d& colour, Refl_t refl):
         center(center), emission(emission), colour(colour), refl(refl), radius(radius) {}
    
    double intersect(const Ray& ray) const{
        Vector4d oc = center - ray.o;
        double b = oc.dot(ray.d);
        double discriminant = b*b - oc.squaredNorm() + radius * radius;
        if(discriminant > 0) // ray intersects with the sphere
        {
            double t = b-std::sqrt(discriminant); // distance between the ray and the centre of the sphere
            if(t > ray.t_min)  return t;
            t = b+std::sqrt(discriminant);
            if(t > ray.t_min)  return t;
        }
        return 0;
    }
    Vector4d center;
    Array3d emission;
    Array3d colour;
    Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
    double radius;
};


// find hit object closest to the screen by iterating over all objects; very inefficient!
inline const Sphere* hit(const std::vector<Sphere>& scene, const Ray& r, double& t) {
    const Sphere* hitobj = nullptr;
    t = r.t_max;
    for (const Sphere& sphere : scene) {
        const double tt = sphere.intersect(r);
        if (tt > 0 && tt < t) {
            t = tt;
            hitobj = &sphere;
        }
    }
    return hitobj;
}

// Gram-Schmidt orthogonalisation
inline Vector4d gs(const Vector4d& v, const Vector4d& orth) {
	return(v - v.dot(orth)*orth).normalized();
}

class PathTracer {
public:
    Array3d radiance(
        const std::vector<Sphere>& scene,
        const Ray& r,
        int depth,
        Sampling& sampler,
        const Configuration& configuration,
        int E = 1
    ) const {
        double t = 0;
        const Sphere* hitobj = hit(scene, r, t);
        // return the background if nothing is hit
        if(!hitobj){
            const double alpha = 0.5 * (r.d[1] + 1.0);
            return configuration.env_map*(Array3d(1.0,1.0,1.0) * (1-alpha) + Array3d(0.5, 0.7, 1.0) * alpha);
        }
        
        // something is hit
        const Vector4d hitpos = r.o + t * r.d;
        const Vector4d normal = (hitpos - hitobj->center).normalized();
        const Vector4d nl = normal.dot(r.d) < 0 ? normal : -normal;
        Array3d f = (*hitobj).colour;
        const double p = f.maxCoeff();
        
        ++depth;
        if (depth > configuration.min_bounce || !p){
            if (sampler.random_double()<p) f=f/p; else return E*(*hitobj).emission;
        }
        if (hitobj->refl == DIFF){
            Vector4d d = random_unit_vector(sampler);
            if(d.dot(nl) < 0) {
                d = -d;
            }
            const Ray newray(hitpos, d); // diffuse reflection
            Array3d e(0, 0, 0); // emitter samling colour
            if (configuration.emitter_sampling_off == 0) {
                // emitter sampling
                for (const Sphere& light : scene) {
                    if (light.emission.maxCoeff() <= 0) {
                        continue; // skip non-lights
                    }
                    const Vector4d tt = light.center - hitpos;
                    const double cos_a_max = std::sqrt(1 - light.radius * light.radius / tt.squaredNorm());
                    const Vector4d su = random_unit_vector(sampler);  // TODO: fix this biased sampling
                    const Vector4d l = (tt + su * light.radius * sampler.random_double()).normalized();
                    const Ray lr(hitpos, l);
                    double light_t = 0;
                    const Sphere* closest = hit(scene, lr, light_t);
                    if (closest == &light && l.dot(nl) > 0) {
                        e += f * (light.emission * l.dot(nl) * 2 * (1 - cos_a_max));
                    }
                }
            }
            // bounce
            return e+E*(*hitobj).emission + f * radiance(scene, newray, depth, sampler, configuration, configuration.emitter_sampling_off);
        }else if (hitobj->refl == SPEC) {
            const Ray newray(hitpos,r.d-2*normal.dot(r.d)*normal);
            return (*hitobj).emission + f * radiance(scene, newray, depth, sampler, configuration);
        }else if (hitobj->refl == REFR) {
            const Ray reflRay(hitpos, r.d-2*normal.dot(r.d)*normal);
            const bool into = normal.dot(nl)>0;                // Ray from outside going in?
            const double nc=1;
            const double nt=1.5;
            const double nnt=into ? nc/nt : nt/nc;
            const double ddn=r.d.dot(nl);
            const double cos2t=1-nnt*nnt*(1-ddn*ddn);
            if (cos2t<0){    // Total internal reflection
                return (*hitobj).emission + f * radiance(scene,reflRay,depth,sampler,configuration);
            }
            const Vector4d tdir = (r.d*nnt - normal*((into?1:-1)*(ddn*nnt+std::sqrt(cos2t)))).normalized();
            const double a=nt-nc;
            const double b=nt+nc;
            const double R0=a*a/(b*b);
            const double c = 1-(into?-ddn:tdir.dot(normal));
            const double Re=R0+(1-R0)*c*c*c*c*c;
            const double Tr=1-Re;
            const double P=.25+.5*Re;
            const double RP=Re/P;
            const double TP=Tr/(1-P);
            Array3d e = (*hitobj).emission;
            const Ray tRay(hitpos,tdir);
            if(depth>2){
                if(sampler.random_double()<P){
                    e += f * radiance(scene,reflRay,depth,sampler,configuration)*RP;
                }else{
                    e += f * radiance(scene,tRay,depth,sampler,configuration)*TP;
                }
            }else{
                e += f * (radiance(scene,reflRay,depth,sampler,configuration) * Re
                         +radiance(scene,tRay,depth,sampler,configuration)*Tr);
            }
            return e;
        }
        return Array3d::Zero();
    }
};
// utility
inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }
inline uint8_t toInt(double x){ return static_cast<uint8_t>(std::pow(clamp(x),1/2.2)*255+.5); }

// for writing png
struct RGBs{
    RGBs(){}
    RGBs(uint8_t r, uint8_t g, uint8_t b):r(r),g(g),b(b){}
    uint8_t r,g,b;
};

static_assert(sizeof(RGBs) == 3, "RGB pixel must be tightly packed.");

inline std::string make_frame_filename(char eye, int frame_index) {
    std::ostringstream ss;
    ss << "output_" << eye << "_"
       << std::setfill('0') << std::setw(3) << std::right
       << frame_index << ".png";
    return ss.str();
}

constexpr double kOpenClSceneScale = 1.0;

inline OpenClFloat4 to_opencl_float4(const Vector4d& v) {
    return OpenClFloat4{
        static_cast<double>(v[0] * kOpenClSceneScale),
        static_cast<double>(v[1] * kOpenClSceneScale),
        static_cast<double>(v[2] * kOpenClSceneScale),
        static_cast<double>(v[3] * kOpenClSceneScale)
    };
}

inline OpenClFloat4 to_opencl_rgb(const Array3d& c) {
    return OpenClFloat4{
        static_cast<double>(c[0]),
        static_cast<double>(c[1]),
        static_cast<double>(c[2]),
        0.0
    };
}

OpenClCamera to_opencl_camera(const Camera& cam) {
    OpenClCamera out{};
    out.origin = to_opencl_float4(cam.origin);
    out.corner = to_opencl_float4(cam.upper_left_corner);
    out.u = to_opencl_float4(cam.u);
    out.v = to_opencl_float4(cam.v);
    return out;
}

std::vector<OpenClSphere> to_opencl_scene(const std::vector<Sphere>& scene) {
    std::vector<OpenClSphere> out;
    out.reserve(scene.size());
    for (const Sphere& sp : scene) {
        out.push_back(OpenClSphere{
            to_opencl_float4(sp.center),
            to_opencl_rgb(sp.emission),
            to_opencl_rgb(sp.colour),
            static_cast<int>(sp.refl),
            static_cast<double>(sp.radius * kOpenClSceneScale)
        });
    }
    return out;
}

// rendering
void render_cpu(const Configuration& configuration, const Camera& cam, const std::vector<Sphere>& scene){
    const int width = configuration.width;
    const int height = configuration.height;
    const int sample_count = configuration.number_of_samples;
    const unsigned progress_step = static_cast<unsigned>(std::max(1, height / 10));
    const unsigned progress_ticks = static_cast<unsigned>((height + static_cast<int>(progress_step) - 1) / static_cast<int>(progress_step));

    ProgressBar pb(50, progress_ticks);
    std::vector<RGBs> bitmap(static_cast<std::size_t>(width) * static_cast<std::size_t>(height));
    PathTracer integrator;
    const double inv_width = 1.0 / static_cast<double>(width);
    const double inv_height = 1.0 / static_cast<double>(height);
    const double sample_weight = 1.0 / static_cast<double>(sample_count);
    
#pragma omp parallel for schedule(dynamic, 1)    // OpenMP
    for(int y=0; y<configuration.height; y++){
        Sampling sampler((y + 1) * (y + 1) * (y + 1));
        for (int x=0; x<width; ++x){
            Array3d col(0,0,0); // initial colour of pixel
            for (int s=0; s<sample_count; ++s)
            {
                const double u = (x + sampler.random_double()) * inv_width;
                const double v = (y + sampler.random_double()) * inv_height;
                const Ray ray = cam.get_ray(u, v);
                col += integrator.radiance(scene, ray, 0, sampler, configuration) * sample_weight;
            }

            const std::size_t index =
                static_cast<std::size_t>(y) * static_cast<std::size_t>(width)
                + static_cast<std::size_t>(x);
            bitmap[index] = RGBs(toInt(col[0]), toInt(col[1]), toInt(col[2]));
        }
		// progress bar
		if ((y % static_cast<int>(progress_step)) == 0) {
			pb.tick();
			pb.print();
		}
    }
    pb.end();
    
    if (!simple_png::write_png_rgb8(
            configuration.output,
            width,
            height,
            reinterpret_cast<const uint8_t*>(bitmap.data()),
            3 * width)) {
        std::cerr << "Failed to write image: " << configuration.output << std::endl;
    }
}

// entry point
int main(int argc, char *argv[]){
    
    int opt;
	Configuration configuration;
    while((opt=getopt(argc,argv,"m:s:w:h:f:g"))!=-1){
        switch(opt){
            case 'm':  // minimum number of bouncing
                configuration.min_bounce = atoi(optarg);
                break;
            case 's':  // number of samples
                configuration.number_of_samples = atoi(optarg);
                break;
            case 'w':  // width of output
                configuration.width = atoi(optarg);
                break;
            case 'h':  // height
                configuration.height = atoi(optarg);
                break;
				case 'f':  // number of frames
					configuration.frames = atoi(optarg);
					break;
            case 'g':  // use GPU backend when available
                configuration.use_gpu = true;
                break;
			case ':':
                std::cout << opt << " needs arguments"<< std::endl;
                break;
            case '?':
                std::cout << "unknown option: " << opt << std::endl;
                break;
        }
    }

    if (configuration.width <= 0 || configuration.height <= 0 || configuration.frames <= 0 ||
        configuration.number_of_samples <= 0 || configuration.min_bounce < 0) {
        std::cerr << "Invalid render configuration. "
                  << "Expected positive width/height/frames/samples and min_bounce >= 0." << std::endl;
        return 1;
    }
    

    const std::vector<SceneSphere> scene_def = build_default_scene();
    std::vector<Sphere> scene;
    scene.reserve(scene_def.size());
    for (const SceneSphere& sp : scene_def) {
        scene.emplace_back(
            sp.radius,
            sp.center,
            sp.emission,
            sp.colour,
            static_cast<Refl_t>(sp.refl)
        );
    }
    

	    // 2D screen axis
	    Vector4d v = Vector4d(0, -1, 0.01, 0.02);

    const std::vector<OpenClSphere> opencl_scene =
        configuration.use_gpu ? to_opencl_scene(scene) : std::vector<OpenClSphere>();

		// render images for both eyes
	    bool warned_gpu_fallback = false;
	    for(int i=0;i<configuration.frames;i++){
        const double theta = M_PI * static_cast<double>(i) / static_cast<double>(configuration.frames);
		const Vector4d origin(100*std::sin(theta), 100-200.0*i/configuration.frames, 190,  100*std::cos(theta));
		Vector4d eye_vect(4*std::cos(theta/3),4*std::sin(theta/3),0,0);
        eye_vect = gs(eye_vect,origin.normalized());
//        std::cout << origin << std::endl << eye_vect << std::endl;
    	// left eye
		const Vector4d left_rel((-(origin+eye_vect).normalized()));
		const Vector4d left_u = gs(eye_vect,left_rel);
		const Vector4d left_v = gs(gs(v,left_rel),left_u);
        Camera cam = Camera(origin+eye_vect, left_rel, left_u,left_v);
        configuration.output = make_frame_filename('L', i);
	        std::cout << configuration.output << std::endl;
	        bool rendered = false;
	        if (configuration.use_gpu) {
            rendered = render_opencl_image(
                configuration.output,
                configuration.width,
                configuration.height,
                configuration.number_of_samples,
                configuration.min_bounce,
                static_cast<double>(configuration.env_map),
                2 * i,
                to_opencl_camera(cam),
                opencl_scene,
                50
            );
	            if (!rendered && !warned_gpu_fallback) {
#if defined(SMALLPT4D_WITH_OPENCL)
	                std::cerr << "GPU render failed. Falling back to CPU." << std::endl;
#else
	                std::cerr << "GPU backend requested but OpenCL is not enabled in this build. Using CPU." << std::endl;
#endif
	                warned_gpu_fallback = true;
	            }
	        }
	        if (!rendered) {
	            render_cpu(configuration, cam, scene);
	        }
    	// right eye
		const Vector4d right_rel((-(origin-eye_vect).normalized()));
		const Vector4d right_u = gs(eye_vect,right_rel);
		const Vector4d right_v = gs(gs(v,right_rel),right_u);
        cam = Camera(origin-eye_vect, right_rel, right_u,right_v);
        configuration.output = make_frame_filename('R', i);
	        std::cout << configuration.output << std::endl;
	        rendered = false;
	        if (configuration.use_gpu) {
            rendered = render_opencl_image(
                configuration.output,
                configuration.width,
                configuration.height,
                configuration.number_of_samples,
                configuration.min_bounce,
                static_cast<double>(configuration.env_map),
                2 * i + 1,
                to_opencl_camera(cam),
                opencl_scene,
                50
            );
	            if (!rendered && !warned_gpu_fallback) {
#if defined(SMALLPT4D_WITH_OPENCL)
	                std::cerr << "GPU render failed. Falling back to CPU." << std::endl;
#else
                std::cerr << "GPU backend requested but OpenCL is not enabled in this build. Using CPU." << std::endl;
#endif
                warned_gpu_fallback = true;
            }
	        }
	        if (!rendered) {
	            render_cpu(configuration, cam, scene);
	        }
    }

    return 0;
}
