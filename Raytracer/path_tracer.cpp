#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "stdlib.h"
#include "tools.cpp"
using namespace std;

int main() {

    // generating the spheres
    Sphere* red_sphere = new Sphere(Vector(0, 1000, 0), 940, Vector(1, 0, 0), "diffuse");
    Sphere* green_sphere = new Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0), "diffuse");
    Sphere* blue_sphere = new Sphere(Vector(0, -1000, 0), 990, Vector(0, 0, 1), "diffuse");
    Sphere* pink_sphere = new Sphere(Vector(0, 0, 1000), 940, Vector(1, 0, 1), "diffuse");
    Sphere* cyan_sphere = new Sphere(Vector(1000, 0, 0), 940, Vector(0, 1, 1), "diffuse");
    Sphere* yellow_sphere = new Sphere(Vector(-1000, 0, 0), 940, Vector(1, 1, 0), "diffuse");
    Sphere* object = new Sphere(Vector(0, 0, 0), 10, Vector(1, 1, 1), "diffuse");
    Sphere* object2 = new Sphere(Vector(20.5, 0, 0), 10, Vector(1, 1, 1), "transparent");
    Sphere* object3 = new Sphere(Vector(-20.5, 0, 0), 10, Vector(1, 1, 1), "mirror");
    Sphere* light_source = new Sphere(Vector(-10, 20, 40), 2, Vector(1, 1, 1), "light");

    //generating the mesh
    TriangleMesh* mesh = new TriangleMesh(Vector(1, 1, 1), "mirror");
    mesh->readOBJ("./cadnav.com_model./Models_F0202A090./cat.obj");
    /*for (int i = 0; i < mesh->vertices.size(); i++) {
        mesh->vertices[i] = 0.6 * mesh->vertices[i] + Vector(0, -10, 0);
    }*/

    //creating the scene
    static Geometry *A[] = {red_sphere, blue_sphere, green_sphere, pink_sphere, cyan_sphere, yellow_sphere, mesh, light_source};
    vector<Geometry *> scene_components(A, A + sizeof(A) / sizeof(A[0]));
    Scene scene(scene_components);

    //camera and pixel grid
    Vector Q = Vector(0, 0, 55);               //camera center
    double W = 600;                            //grid width
    double H = 512;                            //grid height
    double fov = PI/2;                         //alpha, field of view
    int limit = 1;                            //amount of rays
    vector<unsigned char> img(W*H*3);          //image vector

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(H); i++) {

        #pragma omp parallel for
        for (int j = 0; j < static_cast<int>(W); j++) {
            Vector V;
            V[2] = Q[2] - (W / (2 * tan(fov / 2))); 
            Vector sum(0., 0., 0.);

            #pragma omp parallel for
            for (int k = 0; k < limit; k++) {
                Vector M = boxMuller();
                V[0] = Q[0] + 0.5 + j - (W / 2) + M[0];
                V[1] = Q[1] - i - 0.5 + (H / 2) + M[1];
                Vector n = (V - Q) / norm(V - Q);             
                                                                
                //Dof follows
                /* double D = (W / (2 * tan(fov / 2))) * norm(V - Q) / (norm(V - Q) - (W / (2 * tan(fov / 2))));
                Vector P = Q + (D / abs(n[2])) * n;
                Vector Qprime(0., 0., Q[2]);
                double emily = 0.65;     //r_max
                double r = (double(rand()) / double(RAND_MAX)) * emily;
                double arthur = (double(rand()) / double(RAND_MAX)) * 2 * PI;    //teta but with a better name
                Qprime[0] = Q[0] + cos(arthur) * r;
                Qprime[1] = Q[1] + sin(arthur) * r;*/

                Vector color = scene.get_color(Ray(Q, n), 5, *light_source, 200000, false);
                sum += color;

            }
            sum = sum / limit;
            double power = 1. / 2.2; 
            img[(i * W + j) * 3 + 0] = min(255, max(0, int(pow(sum[0], power) * 255)));
            img[(i * W + j) * 3 + 1] = min(255, max(0, int(pow(sum[1], power) * 255)));
            img[(i * W + j) * 3 + 2] = min(255, max(0, int(pow(sum[2], power) * 255)));
        }
    }

    stbi_write_png("cat_mirror3.png", W, H, 3, &img[0], 0);
    return 0;
}