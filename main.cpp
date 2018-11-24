/*
 smallpt4D: simple ray tracer for 4D scenes by Shizuo Kaji
 Based on: smallpt: Global Illumination in 99 lines of C++ by Kevin Beason
 Depends on: Eigen 3 (http://eigen.tuxfamily.org)
             stb_image_write.h by Sean Barrett (included)
*/


#include "stdafx.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#if defined(_WIN32)|| defined(_WIN64)
	#include "getopt_win.h"
#else
	#include <getopt.h>
#endif
#include <mutex>

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
};

// random number generator
class Sampling{
public:
    Sampling(int s): gen(s), zero_one(0,1){}
    double random_double() { return zero_one(gen); }
private:
    std::mt19937 gen;
    std::uniform_real_distribution<> zero_one;
};


struct Ray {
public:
    Ray(Vector4d o_, Vector4d d_, double t_max=std::numeric_limits<double>::infinity()) :
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
    Sphere(double radius, Vector4d center, Array3d emission, Array3d colour, Refl_t refl):
         radius(radius), center(center), emission(emission), colour(colour), refl(refl) {}
    
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
inline Sphere const *hit(const std::vector<Sphere> &scene, Ray &r, double &t) {
    Sphere const *hitobj = 0;
    t = r.t_max;
    for (std::vector<Sphere>::const_iterator i = scene.begin(); i != scene.end() ; ++i) {
        double tt = (*i).intersect(r);
        if (tt > 0 && tt < t) {
            t = tt;
            hitobj = &(*i);
        }
    }
    return hitobj;
}

// Gram-Schmidt orthogonalisation
inline Vector4d gs(Vector4d v, Vector4d orth) {
	return(v - v.dot(orth)*orth).normalized();
}

class PathTracer {
public:
    PathTracer() {}
    Array3d radiance(const std::vector<Sphere> &scene, Ray &r, int depth, Sampling &sampler, const Configuration& configuration, int E=1) const{
        double t=0;
        Sphere const *hitobj =  hit(scene, r, t);
        // return the background if nothing is hit
        if(!hitobj){
            double alpha = 0.5 * (r.d[1] + 1.0);
            return configuration.env_map*(Array3d(1.0,1.0,1.0) * (1-alpha) + Array3d(0.5, 0.7, 1.0) * alpha);
        }
        
        // something is hit
        Vector4d hitpos=r.o + t * r.d;
        Vector4d normal=(hitpos-(*hitobj).center).normalized();
        Vector4d nl=normal.dot(r.d)<0 ? normal : -normal;
        Array3d f = (*hitobj).colour;
        double p = f.maxCoeff();
        
        if (++depth>configuration.min_bounce || !p){
            if (sampler.random_double()<p) f=f/p; else return E*(*hitobj).emission;
        }
        if ((*hitobj).refl == DIFF){
            Vector4d d = Vector4d::Random().normalized();
            if(d.dot(nl)<0) d = -d;
            Ray newray(hitpos, d); // diffuse reflection
            Array3d e(0,0,0); // emitter samling colour
            if(configuration.emitter_sampling_off == 0){
                // emitter sampling
                for (std::vector<Sphere>::const_iterator light = scene.begin(); light != scene.end() ; ++light) {
                    if ((*light).emission.maxCoeff()<=0) continue; // skip non-lights
					Vector4d tt = (*light).center - hitpos;
					double cos_a_max = sqrt(1 - (*light).radius * (*light).radius / tt.squaredNorm());
					Vector4d su = Vector4d::Random().normalized();  // TODO: fix this biased sampling
					Vector4d l = (tt + su * (*light).radius*sampler.random_double()).normalized();
					Ray lr = Ray(hitpos, l);
                    double d=0;
                    Sphere const *closest = hit(scene,lr,d);
                    if(closest == &(*light) && l.dot(nl) > 0){
                        e += f * ((*light).emission * l.dot(nl)*2*(1-cos_a_max));
                    }
                }
            }
            // bounce
            return e+E*(*hitobj).emission + f * radiance(scene, newray, depth, sampler, configuration, configuration.emitter_sampling_off);
        }else if ((*hitobj).refl == SPEC) {
            Ray newray(hitpos,r.d-2*normal.dot(r.d)*normal);
            return (*hitobj).emission + f * radiance(scene, newray, depth, sampler, configuration);
        }else if ((*hitobj).refl == REFR) {
            Ray reflRay(hitpos, r.d-2*normal.dot(r.d)*normal);
            bool into = normal.dot(nl)>0;                // Ray from outside going in?
            double nc=1, nt=1.5, nnt=into ? nc/nt : nt/nc, ddn=r.d.dot(nl);
            double cos2t=1-nnt*nnt*(1-ddn*ddn);
            if (cos2t<0){    // Total internal reflection
                return (*hitobj).emission + f * radiance(scene,reflRay,depth,sampler,configuration);
            }
            Vector4d tdir = (r.d*nnt - normal*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).normalized();
            double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(normal));
            double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
            Array3d e = (*hitobj).emission;
            Ray tRay(hitpos,tdir);
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



// print progress bar
class ProgressBar{
public:
    ProgressBar(int size, int max_ticks): size(size),max_ticks(max_ticks),current(0){}
    
    void end() const { std::cout<<std::endl; }
    void print() {
        std::lock_guard<std::mutex> guard(tick_mutex);
        int pos = size * (float)current/max_ticks;
        std::cout << "[";
        for (int i = 0; i < size; i++) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << ".";
        }
        std::cout << "] " << (int((float)current / max_ticks * 100.0)) << " %" << std::endl;
        std::cout.flush();
    }
    
    void tick() {
        std::lock_guard<std::mutex> guard(tick_mutex);
        if(current < max_ticks) current++;
    }
    
private:
    unsigned size;
    unsigned max_ticks;
    unsigned current;
    std::mutex tick_mutex;
};

// utility
inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }
inline uint8_t toInt(double x){ return uint8_t(pow(clamp(x),1/2.2)*255+.5); }

// for writing png
struct RGBs{
    RGBs(){}
    RGBs(uint8_t r, uint8_t g, uint8_t b):r(r),g(g),b(b){}
    uint8_t r,g,b;
};

// rendering
void render(Configuration& configuration, Camera& cam, std::vector<Sphere>& scene){
    ProgressBar pb(50, configuration.height/10);
    MatrixXd R(configuration.width, configuration.height);
    MatrixXd G(configuration.width, configuration.height);
    MatrixXd B(configuration.width, configuration.height);
    std::unique_ptr<PathTracer> integrator(new PathTracer);
    
#pragma omp parallel for schedule(dynamic, 1)    // OpenMP
    for(int y=0; y<configuration.height; y++){
        Sampling sampler(y*y*y);
        for (unsigned short x=0; x<configuration.width; x++){
            Array3d col(0,0,0); // initial colour of pixel
            for (int s=0; s<configuration.number_of_samples; s++)
            {
                double u = (x + sampler.random_double()) / configuration.width;
                double v = (y + sampler.random_double()) / configuration.height;
                Ray ray = cam.get_ray(u, v);
                col = col + integrator->radiance(scene,ray,0,sampler,configuration)*(1./configuration.number_of_samples);
            }
            R(x,y) = clamp(col[0]);
            G(x,y) = clamp(col[1]);
            B(x,y) = clamp(col[2]);
            
        }
		// progress bar
		if (y % 10 == 0) {
			pb.tick();
			pb.print();
		}
    }
    std::cout<<"\n";
    
    // write image to file
    std::vector<RGBs> bitmap(configuration.width * configuration.height);
    for(int y = 0 ; y < configuration.height ; ++y)
    {
        for (int x=0; x<configuration.width; ++x)
            bitmap[y * configuration.width + x] = RGBs(toInt(R(x,y)), toInt(G(x,y)), toInt(B(x,y)));
    }
    stbi_write_png(configuration.output.c_str(), configuration.width, configuration.height, 3, bitmap.data(), 3 * configuration.width);
}

// entry point
int main(int argc, char *argv[]){
    
    int opt;
	Configuration configuration;
    while((opt=getopt(argc,argv,"m:s:w:h:f:"))!=-1){
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
			case ':':
                std::cout << opt << " needs arguments"<< std::endl;
                break;
            case '?':
                std::cout << "unknown option: " << opt << std::endl;
                break;
        }
    }
    

////////////// modify below /////////////    
// scene construction
    std::vector<Sphere> scene;
    scene.clear();
    scene.reserve(1000);
    scene.push_back(Sphere(1e5, Vector4d(-1e5-175,0,0,0), Array3d::Zero(),Array3d(.75,.25,.25),SPEC)); //Left wall
    scene.push_back(Sphere(1e5, Vector4d(1e5+175,0,0,0),  Array3d::Zero(),Array3d(.25,.25,.75),SPEC));//Right wall
    scene.push_back(Sphere(1e5, Vector4d(0,0, 1e5-175,0),      Array3d::Zero(),Array3d(.25,.75,.25),SPEC));//Back wall
    scene.push_back(Sphere(1e5, Vector4d(0,0,-1e5+200,0),  Array3d::Zero(),Array3d(.75,.75,.75),  SPEC));//Frnt wall
    scene.push_back(Sphere(1e5, Vector4d(0, -1e5-175,0,0),    Array3d(1,1,1),Array3d(.35,.75,.55),SPEC));//floor light
//    scene.push_back(Sphere(1e5, Vector4d(0,-1e5-175,0,0),     Array3d(0,0,0),Array3d(.35,.75,.55),SPEC));//floor
    scene.push_back(Sphere(1e5, Vector4d(0, 1e5+175,0,0),    Array3d(0,0,0),Array3d(.35,.25,.55),SPEC));//ceiling
    scene.push_back(Sphere(1e5, Vector4d(0,0,0,175+1e5),     Array3d::Zero(),Array3d(.35,.55,.75),SPEC));// w+
    scene.push_back(Sphere(1e5, Vector4d(0,0,0,-175-1e5),     Array3d::Zero(),Array3d(.75,.55,.35),SPEC));// w-
        // objects
	int m = 2;
    for(int i=-m; i<=m; i++){
        for(int j=-m; j<=m; j++){
            for(int k=-m; k<=m; k++){
                for(int l=-m; l<=m; l++){
					switch(abs(i+j+k+l) % 3){
						case 0:
							scene.push_back(Sphere(10, Vector4d(i*75/m, j* 75 / m, k* 75 / m, l* 75 / m), Array3d(.50+.2*i/m, .50+.2*j/m, .50+.2*l/m),Array3d(.50+.2*i/m, .50+.2*j/m, .50+.2*l/m), DIFF));
							break;
						case 1:
							scene.push_back(Sphere(10, Vector4d(i * 75 / m, j * 75 / m, k * 75 / m, l * 75 / m), Array3d::Zero(), Array3d(.50 + .2*i / m, .50 + .2*j / m, .50 + .2*l / m), SPEC));
							break;
						case 2:
							scene.push_back(Sphere(10, Vector4d(i * 75 / m, j * 75 / m, k * 75 / m, l * 75 / m), Array3d::Zero(), Array3d(.50 + .2*i / m, .50 + .2*j / m, .50 + .2*l / m), REFR));
							break;
					}
                }
            }
        }
    }
//        Sphere(16.5,Vector4d(50,20.0,80,0),       Array3d::Zero(),Array3d(1,1,1)*.999, REFR),// glass ball
//    scene.push_back(Sphere(16.5,Vector4d(23,-10,90,5),       Array3d::Zero(),Array3d(1,1,1)*.999, SPEC));// mirror ball
//        Sphere(16.5,Vector4d(80,16.5,100,10),       Array3d::Zero(),Array3d(.8,.25,0), DIFF), //red ball
//        Sphere(600, Vector4d(50,681.6-.27,81.6,0), Array3d(1,1,1)*12,  Array3d::Zero(), DIFF) //ceiling light
    

    // 2D screen axis
    Vector4d v = Vector4d(0, -1, 0.01, 0.02);

	// render images for both eyes
    for(int i=0;i<configuration.frames;i++){
        double theta = 0.5*M_PI*i*2.0/configuration.frames;
		Vector4d origin(100*sin(theta), 100-200.0*i/configuration.frames, 190,  100*cos(theta));
		Vector4d eye_vect(4*cos(theta/3),4*sin(theta/3),0,0);
        eye_vect = gs(eye_vect,origin.normalized());
//        std::cout << origin << std::endl << eye_vect << std::endl;
    	// left eye
		Vector4d left_rel((-(origin+eye_vect).normalized()));
		Vector4d left_u = gs(eye_vect,left_rel);
		Vector4d left_v = gs(gs(v,left_rel),left_u);
        Camera cam = Camera(origin+eye_vect, left_rel, left_u,left_v);
        std::stringstream ss;
        ss << "output_L_" << std::setfill('0') << std::setw(3) << std::right << std::to_string(i) << ".png";
        ss >> configuration.output;
        std::cout << configuration.output << std::endl;
        render(configuration, cam, scene);
    	// right eye
		Vector4d right_rel((-(origin-eye_vect).normalized()));
		Vector4d right_u = gs(eye_vect,right_rel);
		Vector4d right_v = gs(gs(v,right_rel),right_u);
        cam = Camera(origin-eye_vect, right_rel, right_u,right_v);
		std::stringstream ssr;
		ssr << "output_R_" << std::setfill('0') << std::setw(3) << std::right << std::to_string(i) << ".png";
        ssr >> configuration.output;
        std::cout << configuration.output << std::endl;
        render(configuration, cam, scene);
    }
}
