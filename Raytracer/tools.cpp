#define PI 3.14159265
#include "mesh.cpp"
#include <cstdlib>
#include <iostream>
#include <random>

static std::default_random_engine engine(10); 
static std::uniform_real_distribution<double> uniform(0, 1);
double epsilon = 0.00001;   //for noise reduction later on

//useful functions
Vector random_cos(const Vector &N) {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    double rad = sqrt(1 - r2);
    double x = cos(2 * PI * r1) * rad;
    double y = sin(2 * PI * r1) * rad;
    double z = sqrt(r2);

    double min_comp = min(min(abs(N[0]), abs(N[1])), abs(N[2]));
    Vector t1;
    for (int i = 0; i < 3; i++) {
        if (abs(N[i]) == min_comp) {
            switch (i) {
            case 0:
                t1 = Vector(0, -N[2], N[1]);
                break;
            case 1:
                t1 = Vector(-N[2], 0, N[0]);
                break;
            case 2:
                t1 = Vector(-N[1], N[0], 0);
                break;    
            }
            break;
        }
    }
    t1 = t1 / sqrt(max(dot(t1, t1), 0.));
    Vector t2 = cross(N, t1);
    return x * t1 + y * t2 + z * N;
}

Vector boxMuller() {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    double x = sqrt(-2 * log(r1)) * cos(2 * PI * r2);
    double y = sqrt(-2 * log(r1)) * sin(2 * PI * r2);
    return Vector(x, y, 0.);
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

        intersection.P = r.O + intersection.t * r.u;
        intersection.N = (intersection.P - C) / sqrt(max(dot(intersection.P - C, intersection.P - C), 0.));  
    }
    return intersection;
}

//Mesh functions
Intersection TriangleMesh::intersect(const Ray& r) {
    Intersection intersection;
    double d = 100000000;
    for (int i = 0; i < indices.size(); i++) {
        TriangleIndices triangle = indices[i];
        Vector A = vertices[triangle.vtxi];
        Vector B = vertices[triangle.vtxj];
        Vector C = vertices[triangle.vtxk];
        Vector e1 = B - A;
        Vector e2 = C - A;
        Vector N = cross(e1, e2);
        double dotun = dot(r.u, N);
        Vector cross_aminou = cross(A - r.O, r.u);

        double b = dot(e2, cross_aminou) / dotun;
        double c = - dot(e1, cross_aminou) / dotun;
        double a = 1. - b - c;

        intersection.t = dot(A - r.O, N) / dotun;

        if (intersection.t >= 0. && 0. <= a && a <= 1. && 0. <= b && b <= 1. && 0. <= c && c <= 1.) {
            intersection.exists = true;
            if (intersection.t < d) {
                d = intersection.t;
                intersection.P = r.O + intersection.t * r.u;
                intersection.N = a * normals[triangle.ni] + b * normals[triangle.nj] + c * normals[triangle.nk];
                intersection.N = intersection.N / norm(intersection.N);
            }
        }
    }
    if (d == 100000000) {
        intersection.exists = false;
    }
    return intersection;
}

//Scene functions

int Scene::closest_intersect(const Ray& r) {                   //returns the index of the closest object that intersects the ray r
    int min = 100000;
    int closest = -1;
    for (int i = 0; i < objects.size(); i++) {
        if (objects[i]->intersect(r).exists && objects[i]->intersect(r).t < min) {
            min = objects[i]->intersect(r).t;
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

    un = objects[i]->intersect(r);
    return un;
}

Vector Scene::get_color(const Ray& ray, int ray_depth, Sphere light, double intensity, bool last_bounce_diffuse) {
    if (ray_depth < 0) {
        return Vector(0., 0., 0.);
    }

    int object_id = closest_intersect(ray);
    Intersection inter = this->intersection(ray);

    if (inter.exists) {
        Vector N = inter.N;
        Vector P = inter.P + epsilon * N;

        //light source
        if (objects[object_id]->type == "light") {
            if (last_bounce_diffuse) {
                return Vector(0., 0., 0.);
            }
            return (intensity / (4 * PI * PI * light.R * light.R)) * Vector(1., 1., 1.);
        }

        //reflection
        else if (objects[object_id]->type == "mirror") {                                   
            Ray reflected = Ray(P, ray.u - (2 * dot(ray.u, N)) * N);
            return get_color(reflected, ray_depth - 1, light, intensity, false);
        }

        //refraction
        else if (objects[object_id]->type == "transparent") {       
            double n1 = 1;
            double n2 = 1.5;     
            P = P - 2 * epsilon * N;   

            double k0 = (n1 - n2) * (n1 - n2)/((n1 + n2) * (n1 + n2));
            double R = k0 + (1 - k0) * pow((1 - abs(dot(N, ray.u))), 5);
            double T = 1 - R;
            double ran = double(rand()) / double(RAND_MAX);

            if (ran < R) {
                Ray reflected = Ray(P, ray.u - (2 * dot(ray.u, N)) * N);
                return get_color(reflected, ray_depth - 1, light, intensity, false);
            }
            else {
                if (dot(ray.u, N) > 0) {
                    P = P + 2 * epsilon * N;
                    N =  -N;
                    n1 = 1.5;
                    n2 = 1;
                }   
                Vector wt = (n1 / n2)*(ray.u - dot(ray.u, N)*N);   //tangential component of direction
                double radic = 1 - (n1/n2)*(n1/n2)*(1 - dot(ray.u, N)*dot(ray.u, N));
                if (radic < 0) {
                    Ray reflected = Ray(P, ray.u - (2 * dot(ray.u, N)) * N);
                    return get_color(reflected, ray_depth - 1, light, intensity, false);
                }
                Vector wn = -(sqrt(radic) * N);   //normal component of direction
                Vector w = wn + wt;
                Ray refracted = Ray(P, w);
                return get_color(refracted, ray_depth - 1, light, intensity, false);
            }
        }

        //diffuse objects
        else {
            //direct lighting
            Vector D = (P - light.C) / norm(P - light.C);
            Vector xprime = light.R * random_cos(D) + light.C;
            Vector Nprime = (xprime - light.C) / norm(xprime - light.C);  
            double d = norm(xprime - P);
            Vector omega_i = (xprime - P) / d;
            double pdf = max(dot(Nprime, D), 0.) / (PI * light.R * light.R);
            Ray r = Ray(P, omega_i);

            double visible;
            if (!this->intersection(r).exists) {visible = 1;}
            else {  
                if (this->intersection(r).t + epsilon > d) {
                    visible = 1;
                }
                else { visible = 0;}
            }

            Vector color(0., 0., 0.);
            if (pdf != 0.) {
                color = (intensity / (4 * PI * PI * PI * light.R * light.R)) * visible * max(dot(N, omega_i), 0.) * max(dot(Nprime, -omega_i), 0.) / (sq_norm(xprime - P) * pdf) * objects[object_id]->albedo;
            }

            //indirect lighting
            Ray random_ray = Ray(P, random_cos(N));
            color += get_color(random_ray, ray_depth - 1, light, intensity, true) * objects[object_id]->albedo;
            return color / 2;
        }
    }
}