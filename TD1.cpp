#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "stdlib.h"
#include "tools.cpp"

int main() {
    double n1 = 1;
    double n2 = 1.5;

    //generating the spheres
    Sphere red_sphere(Vector(0, 1000, 0), 940, Vector(1, 0, 0), false, false);
    Sphere green_sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0), false, false);
    Sphere blue_sphere(Vector(0, -1000, 0), 990, Vector(0, 0, 1), false, false);
    Sphere pink_sphere(Vector(0, 0, 1000), 940, Vector(1, 0, 1), false, false);
    Sphere cyan_sphere(Vector(1000, 0, 0), 940, Vector(0, 1, 1), false, false);
    Sphere yellow_sphere(Vector(-1000, 0, 0), 940, Vector(1, 1, 0), false, false);
    Sphere object(Vector(0, 0, 0), 10, Vector(1, 1, 1), true, false);
    Sphere object2(Vector(20.5, 0, 0), 10, Vector(1, 1, 1), false, true);
    Sphere object3(Vector(-20.5, 0, 0), 10, Vector(1, 1, 1), false, true);
    Vector light_source = Vector(-10, 20, 40);

    //creating the scene
    static Sphere A[] = {red_sphere, blue_sphere, green_sphere, pink_sphere, cyan_sphere, yellow_sphere, object, object2, object3};
    vector<Sphere> scene_components(A, A + sizeof(A) / sizeof(A[0]));
    Scene scene(scene_components);

    //camera and pixel grid
    Vector Q = Vector(0, 0, 55);                //camera center
    double W = 600;                             //grid width
    double H = 512;                             //grid height
    double fov = PI/3;                          //alpha, field of view
    vector<unsigned char> img(W*H*3);           //image vector

    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector V;
            V[0] = Q[0] + 0.5 + j - (W / 2);
            V[1] = Q[1] - i - 0.5 + (H / 2);
            V[2] = Q[2] - (W / (2 * tan(fov / 2))); 
            Vector n = (V - Q) / sqrt(dot(V - Q, V - Q));                       //normalized ray direction

            Vector color = scene.get_color(Ray(Q, n), 2, light_source);

            double power = 1./2.2; 
            img[(i*W+j)*3+0] = min(255.,max(0., pow(color[0], power)*255));
            img[(i*W+j)*3+1] = min(255.,max(0., pow(color[1], power)*255));
            img[(i*W+j)*3+2] = min(255.,max(0., pow(color[2], power)*255));
        }
    }
    stbi_write_png("test.png", W, H, 3, &img[0], 0);
    return 0;
}