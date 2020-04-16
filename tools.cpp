#define PI 3.14159265
#include "classes.cpp"

//Diffuse objects light function

Vector intensity(Scene scene, Vector sigma, Vector P, Vector N, double I, Vector S) {
    double V_p;
    double epsilon = 0.001;
    P += epsilon*N;
    double d = sqrt(dot(S - P, S - P));
    Vector omega = (S - P) / d;
    Ray r = Ray(S, Vector(0, 0, 0) - omega);

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
        Vector P = inter.P + 0.001*N;       //offset for output refinement

        if (spheres[sphere_id].mirror) {
            Ray reflected = Ray(P, ray.u - (2 * dot(ray.u, N)) * N);
            return get_color(reflected, ray_depth - 1, light);
        }

        else {
            Vector color = intensity(*this, spheres[sphere_id].albedo, P, N, 100000, light);
            return color;
        }
    }
}

