#include <math.h> 
#include <cmath>
#include <vector>
#include <string>
#include "vector.cpp"
using namespace std;

class Intersection {
    public:
    bool exists;            //does the ray intersect the sphere
    double t1;              //distance to first intersection
    double t2;              //distance to second intersection
    double t;               //distance to relevant intersection (t1 or t2)
    Vector P;               //intersection point (relevant intersection)
    Vector N;               //unit normal at P
};

class Ray {
    public:
    Vector O;       //origin vector
    Vector u;       //direction unit vector 

    Ray(Vector O, Vector u) {
        this->O = O;
        this->u = u;
    }
};

class Sphere {
    public:
    Vector C;           //center C
    double R;           //radius R
    Vector albedo;      //color RGB 0-1
    string type;
    Intersection intersect(const Ray& r);

    Sphere(Vector C, double R, Vector albedo, string type) {
        this->C = C;
        this->R = R;
        this->albedo = albedo;
        this->type = type;
    }
};

class Scene {
    public:
    vector<Sphere> spheres;
    int closest_intersect(const Ray& r);
    Intersection intersection(const Ray& r);
    Vector get_color(const Ray& ray , int ray_depth, Scene scene, Sphere light, double intensity, bool last_bounce_diffuse);

    Scene(vector<Sphere> spheres) {
        this->spheres = spheres;
    }
};