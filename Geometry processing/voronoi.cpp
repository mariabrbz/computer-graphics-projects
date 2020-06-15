#include "svg_polygon.cpp"

Vector intersect(Vector u, Vector v, vector<Vector> edge) {
    Vector A = edge[0];
    Vector B = edge[1];
    Vector M = Vector((u[0] + v[0]) / 2, (u[1] + v[1]) / 2, 0.);
    int t = dot(M - A, u - v) / dot(B - A, u - v);
    Vector P(0., 0., 0.);
    
    if (t >= 0 && t <= 1) {
        P = A + t * (B - A);
    }

    return P;
}

bool inside(Vector P, Vector u, Vector v) {
    Vector M = Vector((u[0] + v[0]) / 2, (u[1] + v[1]) / 2, 0.);
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
            Vector intersection = intersect(previous_vertex, current_vertex, clip_edge);

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

vector<Polygon> voronoi( ) {
    
}

int main() {
    
    return 0;
}