#define PI 3.14159265
#include "classes.cpp"
#include <cstdlib>
#include <iostream>

double epsilon = 0.001;   //for noise reduction later on

//Diffuse objects light function

Vector intensity(Scene scene, Vector sigma, Vector P, Vector N, double I, Vector S) {
    double V_p;
    P += epsilon*N;
    double d = sqrt(dot(S - P, S - P));
    Vector omega = (S - P) / d;
    Ray r = Ray(P, omega);

    if (!scene.intersection(r).exists) { V_p = 1;}
    else {  
        if (scene.intersection(r).t > d) {
            V_p = 1;
        }
        else { V_p = 0;}
    }
    return (I / (4*PI*PI*d*d)) * V_p * max(dot(N, omega), 0.) * sigma;
}

//Sphere functions

Intersection Sphere::intersect(const Ray& r) {                                     //returns the intersection of sphere with r
    Intersection intersection;                                  
    Vector V = r.O - C;
    Vector minus_V = C - r.O;
    double delta = dot(r.u, V)*dot(r.u, V) - dot(V,  V) + R*R; 
    intersection.exists = (delta >= 0);

    if (intersection.exists) {
        double sq_root = sqrt(delta);
        double dot_product = dot(r.u, minus_V);
        intersection.t1 = dot_product - sq_root;
        intersection.t2 = dot_product + sq_root;

        if (intersection.t2 < 0) { intersection.exists = false;}
        else {
            if (intersection.t1 >= 0) { intersection.t = intersection.t1;}
            else { intersection.t = intersection.t2;}
        }

        intersection.P = r.O + intersection.t*r.u;
        intersection.N = (intersection.P - C) / sqrt(dot(intersection.P - C, intersection.P - C));  
    }
    return intersection;
}

//Scene functions

int Scene::closest_intersect(const Ray& r) {                   //returns the index of the closest sphere that intersects the ray r
    int min = 100000;
    int closest = -1;
    for (int i = 0; i < spheres.size(); i++) {
        if (spheres[i].intersect(r).exists && spheres[i].intersect(r).t < min) {
            min = spheres[i].intersect(r).t;
            closest = i;
        }
    }
    return closest;
}

Intersection Scene::intersection(const Ray& r) {
    Intersection un;
    int i = this->closest_intersect(r);
    if (i == -1) {
        un.exists = false;
        return un;
    }
    Sphere closest_sphere = spheres[i];
    un = closest_sphere.intersect(r);
    return un;
}

Vector Scene::get_color(const Ray& ray , int ray_depth, Vector light) {
    if (ray_depth < 0) {
        return Vector(0., 0., 0.);
    }

    int sphere_id = closest_intersect(ray);
    Intersection inter = spheres[sphere_id].intersect(ray);

    if (inter.exists) {
        Vector N = inter.N;
        Vector P = inter.P + epsilon*N;

        //reflection
        if (spheres[sphere_id].mirror) {                                   
            Ray reflected = Ray(P, ray.u - (2 * dot(ray.u, N)) * N);
            return get_color(reflected, ray_depth - 1, light);
        }

        //refraction
        else if (spheres[sphere_id].transparent) {       
            double n1 = 1;
            double n2 = 1.5;     
            P = P - 2*epsilon*N;   

            double k0 = (n1 - n2) * (n1 - n2)/((n1 + n2) * (n1 + n2));
            double R = k0 + (1 - k0) * pow((1 - abs(dot(N, ray.u)) ), 5);
            double T = 1 - R;
            double ran = double(rand()) / double(RAND_MAX);

            if (ran < R) {
                Ray reflected = Ray(P, ray.u - (2 * dot(ray.u, N)) * N);
                return get_color(reflected, ray_depth - 1, light);
            }
            else {
                if (dot(ray.u, N) > 0) {
                    P = P + 2*epsilon*N;
                    N = Vector(0., 0., 0.) - N;
                    n1 = 1.5;
                    n2 = 1;
                }   
                Vector wt = (n1 / n2)*(ray.u - dot(ray.u, N)*N);   //tangential component of direction
                Vector wn = Vector(0., 0., 0.) - (sqrt(1 - (n1/n2)*(n1/n2)*(1 - dot(ray.u, N)*dot(ray.u, N))) * N);   //normal component of direction
                Vector w = wn + wt;
                Ray refracted = Ray(P, w);
                return get_color(refracted, ray_depth - 1, light);
            }
        }

        //diffuse objects
        else {
            Vector color = intensity(*this, spheres[sphere_id].albedo, P, N, 100000, light);    
            return color;
        }
    }
}