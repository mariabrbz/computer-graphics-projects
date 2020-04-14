#include <math.h> 
#include <vector>
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

    Sphere(Vector C, double R, Vector albedo) {
        this->C = C;
        this->R = R;
        this->albedo = albedo;
    }

    Intersection intersect(Ray r) {                                     //returns the intersection of sphere with r
        Intersection intersection;                                  
        Vector V = r.O - C;
        Vector minus_V = C - r.O;
        double delta = dot(r.u, V)*dot(r.u, V) - dot(V,  V) + R*R; 
        intersection.exists = (delta >= 0);

        if (intersection.exists) {
            intersection.t1 = dot(r.u, minus_V) - sqrt(delta);
            intersection.t2 = dot(r.u, minus_V) + sqrt(delta);

            if (intersection.t2 < 0) { intersection.exists = false;}
            else {
                if (intersection.t1 >= 0) { intersection.t = intersection.t1;}
                else { intersection.t = intersection.t2;}
            }
        }

        intersection.P = r.O + intersection.t*r.u;
        intersection.N = (intersection.P - C) / sqrt(dot(intersection.P - C, intersection.P - C));
        return intersection;
    }
};

class Scene {
    public:
    vector<Sphere> spheres;

    Scene(vector<Sphere> spheres) {
        this->spheres = spheres;
    }

    int closest_intersect(Ray r) {                   //returns the index of the closest sphere that intersects the ray r
        int min = pow(10, 7);
        int closest = -1;
        for (int i = 0; i < spheres.size(); i++) {
            if (spheres[i].intersect(r).exists && spheres[i].intersect(r).t < min) {
                min = spheres[i].intersect(r).t;
                closest = i;
            }
        }
        return closest;
    }

    Intersection intersection(Ray r) {
        Intersection un;
        int i = this->closest_intersect(r);
        if (i == -1) {
            un.exists = false;
            return un;
        }
        un.exists = true;
        Sphere closest_sphere = spheres[i];
        un = closest_sphere.intersect(r);
        return un;
    }
};