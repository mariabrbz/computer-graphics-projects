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
    Ray r = Ray(P, omega);
    if (!scene.intersection(r).exists) { V_p = 1;}
    else {  
        if (scene.intersection(r).t > d) {
            V_p = 1;
        }
        else { V_p = 0;}
    }
    Vector res = (I / (4*PI*PI*d*d)) * V_p * max(dot(N, omega), 0.) * sigma;
    return res;
}

int main(int argc, char **argv) {

    //generating the spheres
    Sphere red_sphere(Vector(0, 1000, 0), 940, Vector(1, 0, 0));
    Sphere green_sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0));
    Sphere blue_sphere(Vector(0, -1000, 0), 990, Vector(0, 0, 1));
    Sphere pink_sphere(Vector(0, 0, 1000), 940, Vector(1, 0, 1));
    Sphere object(Vector(0, 0, 0), 10, Vector(1, 1, 1));
    Vector light_source = Vector(-10, 20, 40);

    //creating the scene
    vector<Sphere> scene_components;            //this probably needs to be optimized
    scene_components.push_back(red_sphere);
    scene_components.push_back(green_sphere); 
    scene_components.push_back(blue_sphere);
    scene_components.push_back(pink_sphere); 
    scene_components.push_back(object);

    Scene scene(scene_components);

    //camera and pixel grid
    Vector Q = Vector(0, 0, 55);              //camera center
    double W = 600;                            //grid width
    double H = 600;                            //grid height
    double fov = PI/3;                          //alpha, field of view
    vector< vector<Vector> > grid;            //grid of pixels  
    vector<unsigned char> img;

    for (int i = 0; i < H; i++) {
        vector<Vector> row;
        for (int j = 0; j < W; j++) {
            Vector V;
            V[0] = Q[0] + 0.5 + j - (W / 2);
            V[1] = Q[1] - i - 0.5 + (H / 2);
            V[2] = Q[2] - (W / (2 * tan(fov / 2))); 
            Vector n = (V - Q) / sqrt(dot(V - Q, V - Q));                       //normalized ray direction

            /*Intersection x = scene.intersection(Ray(Q, n));
            int index = scene.closest_intersect(Ray(Q, n));
            Vector color = intensity(scene, scene.spheres[index].albedo, x.P, x.N, 38000, light_source);
            img.push_back(color[0]*255);
            img.push_back(color[1]*255);
            img.push_back(color[2]*255);*/

            if (object.intersect(Ray(Q, n)).exists) {
                img.push_back(255);
                img.push_back(255);
                img.push_back(255);
            }
            else { 
                img.push_back(0);
                img.push_back(0);
                img.push_back(0);
            }
            row.push_back(n);
        }
        grid.push_back(row);
    }
    stbi_write_png("first_image.png", W, H, 3, &img[0], 0);
}

