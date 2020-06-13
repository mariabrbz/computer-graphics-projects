#define PI 3.14159265
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb-master/stb_image.h"
#include "stb-master/stb_image_write.h"
#include "../Raytracer/vector.cpp"
#include <utility>
#include <algorithm>
#include <vector>
#include <random>
#include <iostream>
using namespace std;

static std::default_random_engine engine(10); 
static std::uniform_real_distribution<double> uniform(0, 1);

Vector random_direction() {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    double rad = sqrt(r2 * (1 - r2));
    double x = cos(2 * PI * r1) * rad;
    double y = sin(2 * PI * r1) * rad;
    double z = 1 - 2 * r2;
    Vector v(x, y, z);
    return v;
}

void color_match(unsigned char* I, unsigned char* M, int n) {
    vector<pair<double, int> > proj_I(n), proj_M(n);

    for (int k = 0; k < 2; k++) {
        Vector v = random_direction();
        int j = 0;
        for (int i = 0; i < 3 * n; i += 3) {
            proj_I[j].first = dot(Vector(I[i], I[i+1], I[i+2]), v);
            proj_I[j].second = i;
            proj_M[j].first = dot(Vector(M[i], M[i+1], M[i+2]), v);
            proj_M[j].second = i;
            j++;
        }

        sort(proj_I.begin(), proj_I.end());
        sort(proj_M.begin(), proj_M.end());

        for (int i = 0; i < n; i++) {
            Vector a = (proj_M[i].first - proj_I[i].first) * v;
            int b = proj_I[i].second;
            I[b] = min(255, max(0, int(I[b] + a[0])));
            I[b+1] = min(255, max(0, int(I[b+1] + a[1])));
            I[b+2] = min(255, max(0, int(I[b+2] + a[2])));
        }
    }
}

int main() {
    int width, height, channels;
    unsigned char *I = stbi_load("lego.png", &width, &height, &channels, 3);
    cout<<I[0]<<endl;
    unsigned char *M = stbi_load("model2.png", &width, &height, &channels, 3);
    int n = 40000;
    color_match(I, M, n);
    stbi_write_png("test2.png", 200, 200, 3, &I[0], 600);
    return 0;
}