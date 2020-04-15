#define STB_IMAGE_WRITE_IMPLEMENTATION
#define PI 3.14159265
#include "stb_image_write.h"
#include <math.h> 
#include <vector>
#include "tools.cpp"
#include <iostream>
using namespace std;

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

int main() {

    //generating the spheres
    Sphere red_sphere(Vector(0, 1000, 0), 940, Vector(1, 0, 0));
    Sphere green_sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0));
    Sphere blue_sphere(Vector(0, -1000, 0), 990, Vector(0, 0, 1));
    Sphere pink_sphere(Vector(0, 0, 1000), 940, Vector(1, 0, 1));
    Sphere cyan_sphere(Vector(1000, 0, 0), 940, Vector(0, 1, 1));
    Sphere yellow_sphere(Vector(-1000, 0, 0), 940, Vector(1, 1, 0));
    Sphere object(Vector(0, 0, 0), 10, Vector(1, 1, 1));
    Vector light_source = Vector(-10, 20, 40);

    //creating the scene
    static Sphere A[] = {red_sphere, blue_sphere, green_sphere, pink_sphere, cyan_sphere, yellow_sphere, object};
    vector<Sphere> scene_components(A, A + sizeof(A) / sizeof(A[0]));
    Scene scene(scene_components);

    //camera and pixel grid
    Vector Q = Vector(0, 0, 55);                //camera center
    double W = 1000;                             //grid width
    double H = 1000;                             //grid height
    double fov = PI/1.5;                          //alpha, field of view
    vector<unsigned char> img(W*H*3);                  //image vector

    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector V;
            V[0] = Q[0] + 0.5 + j - (W / 2);
            V[1] = Q[1] - i - 0.5 + (H / 2);
            V[2] = Q[2] - (W / (2 * tan(fov / 2))); 
            Vector n = (V - Q) / sqrt(dot(V - Q, V - Q));                       //normalized ray direction

            Intersection x = scene.intersection(Ray(Q, n));
            int index = scene.closest_intersect(Ray(Q, n));
            Vector color = intensity(scene, scene.spheres[index].albedo, x.P, x.N, 100000, light_source);

            double power = 1./2.2; 
            img[(i*W+j)*3+0] = min(255.,max(0., pow(color[0], power)*255));
            img[(i*W+j)*3+1] = min(255.,max(0., pow(color[1], power)*255));
            img[(i*W+j)*3+2] = min(255.,max(0., pow(color[2], power)*255));
        }
    }
    stbi_write_png("test.png", W, H, 3, &img[0], 0);
    return 0;
}

