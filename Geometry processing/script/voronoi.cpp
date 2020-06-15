#include "svg_polygon.cpp"
#include <random>
using namespace std;

static default_random_engine engine(10); 
static uniform_real_distribution<double> uniform(0, 1);
static uniform_real_distribution<double> uniform_modified(0, log(100)/100);

vector<Vector> random_pointcloud(int n) {
    vector<Vector> point_cloud(n);
    for (int i = 0; i < n; i++) {
        point_cloud[i] = Vector(uniform(engine), uniform(engine), 0.);
    }
    return point_cloud;
}

Vector intersect(Vector A, Vector B, Vector u, Vector v) {
    Vector M = (u + v) / 2;
    double t = dot(M - A, u - v) / dot(B - A, u - v);
    Vector P(0., 0., 0.);
    
    if (t >= 0 && t <= 1) {
        P = A + t * (B - A);
    }

    return P;
}

Vector intersect_weighted(Vector A, Vector B, Vector u, Vector v, double wu, double wv) {
    Vector M = (u + v) / 2;
    M += ((wu - wv) / (2 * sq_norm(u - v))) * (v - u);
    double t = dot(M - A, u - v) / dot(B - A, u - v);
    Vector P(0., 0., 0.);
    
    if (t >= 0 && t <= 1) {
        P = A + t * (B - A);
    }

    return P;
}

bool inside(Vector P, Vector u, Vector v) {
    Vector M = (u + v) / 2;
    return (dot(P - M, v - u) < 0);
}

bool inside_weighted(Vector P, Vector u, Vector v, double wu, double wv) {
    Vector M = (u + v) / 2;
    M += ((wu - wv) / (2 * sq_norm(u - v))) * (v - u);
    return (dot(P - M, v - u) < 0);
}

Polygon sutherland_hodgman(Polygon subject, Polygon clip) {
    Polygon* out;
    for (int i = 0; i < clip.vertices.size() - 1; i++) {
        vector<Vector> clip_edge(2);
        clip_edge[0] = clip.vertices[i];
        clip_edge[1] = (i == clip.vertices.size() - 1) ? clip.vertices[0] : clip.vertices[i + 1];

        out = new Polygon();
        for(int j = 0; j < subject.vertices.size(); j++) {
            Vector current_vertex = subject.vertices[j];
            Vector previous_vertex = subject.vertices[(j > 0) ? (j - 1) : subject.vertices.size() - 1];
            Vector intersection = intersect(previous_vertex, current_vertex, clip_edge[0], clip_edge[1]);

            if (inside(current_vertex, clip_edge[0], clip_edge[1])) {
                if(!inside(previous_vertex, clip_edge[0], clip_edge[1])) {
                    out->vertices.push_back(intersection);
                }
                out->vertices.push_back(current_vertex);
            }

            else if (inside(previous_vertex, clip_edge[0], clip_edge[1])) {
                out->vertices.push_back(intersection);
            }
        }

        subject = *out;
    }
    return *out;
}

vector<Polygon> voronoi(Polygon subject, Polygon clip) {
    vector<Polygon> res;

    #pragma omp parallel for
    for (int i = 0; i < subject.vertices.size(); i++) {
        Polygon current = clip;

        #pragma omp parallel for
        for (int j = 0; j < subject.vertices.size(); j++) {
            if (i != j) {
                Polygon out;

                #pragma omp parallel for
                for (int k = 0; k < current.vertices.size(); k++) {
                    Vector current_vertex = current.vertices[k];
                    Vector previous_vertex = current.vertices[(k > 0) ? (k - 1) : current.vertices.size() - 1];
                    Vector intersection = intersect(previous_vertex, current_vertex, subject.vertices[i], subject.vertices[j]);

                    if (inside(current_vertex, subject.vertices[i], subject.vertices[j])) {
                        if(!inside(previous_vertex, subject.vertices[i], subject.vertices[j])) {
                            out.vertices.push_back(intersection);
                        }
                        out.vertices.push_back(current_vertex);
                    }

                    else if (inside(previous_vertex, subject.vertices[i], subject.vertices[j])) {
                        out.vertices.push_back(intersection);
                    }
                }
            current = out;
            }
        }
        res.push_back(current);
    }
    return res;
}

vector<Polygon> voronoi_weighted(Polygon subject, Polygon clip, vector<double> weights) {
    vector<Polygon> res;

    #pragma omp parallel for
    for (int i = 0; i < subject.vertices.size(); i++) {
        Polygon current = clip;

        #pragma omp parallel for
        for (int j = 0; j < subject.vertices.size(); j++) {
            if (i != j) {
                Polygon out;

                #pragma omp parallel for
                for (int k = 0; k < current.vertices.size(); k++) {
                    Vector current_vertex = current.vertices[k];
                    Vector previous_vertex = current.vertices[(k > 0) ? (k - 1) : current.vertices.size() - 1];
                    Vector intersection = intersect_weighted(previous_vertex, current_vertex, subject.vertices[i], subject.vertices[j], weights[i], weights[j]);

                    if (inside_weighted(current_vertex, subject.vertices[i], subject.vertices[j], weights[i], weights[j])) {
                        if(!inside_weighted(previous_vertex, subject.vertices[i], subject.vertices[j], weights[i], weights[j])) {
                            out.vertices.push_back(intersection);
                        }
                        out.vertices.push_back(current_vertex);
                    }

                    else if (inside_weighted(previous_vertex, subject.vertices[i], subject.vertices[j], weights[i], weights[j])) {
                        out.vertices.push_back(intersection);
                    }
                }
            current = out;
            }
        }
        res.push_back(current);
    }
    return res;
}

int main() {
    int n = 100;
    vector<Vector> point_cloud = random_pointcloud(n);
    vector<double> weights(n);
    for (int i = 0; i < n; i++) {
        weights[i] = uniform_modified(engine);
    }
    vector<Vector> rectangle{Vector(0, 0, 0), Vector(0, 1, 0), Vector(1, 1, 0), Vector(1, 0, 0)};
    Polygon subject(point_cloud);
    Polygon clip(rectangle);

    save_svg(voronoi_weighted(subject, clip, weights), subject, "100points_weighted.svg") ;
    return 0;
}